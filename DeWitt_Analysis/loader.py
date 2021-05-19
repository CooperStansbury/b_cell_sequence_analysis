import numpy as np
import pandas as pd
import sys
import os
from scipy.spatial.distance import pdist
import Levenshtein
import itertools
import multiprocessing as mp


def build_distance_matrix(sample):
    """A function to build a pairwise distance matrix from sequences

    args:
        : sample (pd.DataFrame): data frame with `vSegment` column

    returns:
        : dists (np.array): a distance matrix
    """
    sequences = sample['nucleotide'].tolist()
    N = len(sequences)
    A = np.zeros((N,N), dtype=np.float64)

    for i in range(0, N):
        for j in range(0, N):
            levy_dist = leven_wrapper(sequences[i],sequences[j])
            A[i,j] = levy_dist

    np.fill_diagonal(A, 1)
    return A


def leven_wrapper(s1, s2):
    """A function to wrap levenshietn distance to handle pandas input
    
    args:
        : s1 (str): string 1
        : s2 (str): string 2
        
    returns:
        : distance (float): distance between two strings
    """
    s1 = str(s1)
    s2 = str(s2)
    return Levenshtein.distance(s1, s2)


class SequenceSampler():
    """A class to manage building samples from a directory of tsv files
    """
    
    def __init__(self, root_dir):
        
        self.root_dir = root_dir
        self.data_dirs = ['D1-M','D1-Na','D1-Nb',
                          'D2-N','D2-M',
                          'D3-M','D3-N']
        
        self.samples = None
        self.data_dir = None
        self.data = None
        self.matrices = None

        
    def _slice_V_segement(self, row):
        """A function to create a new column in a dataframe with the v segement
        
        args:
            : row (row of pd.DataFrame): input dataframe row
            
        returns:
            : df (pd.DataFrame): input dataframe with new column
        """
        seq = row['nucleotide']
        start = row['vIndex']
        end = start + 130
        vSegement = seq[start:end]
        return vSegement
        

    def _read_file(self, filepath, n_rows):
        """A function to read an individial tsv file
        
        args:
            : filepath (str): the file path
            : n_rows (int or None): if int, will read only the first n rows
            
        returns:
            : df (pd.DataFrame): a dataframe
        """
        
        dtypes = {
            'nucleotide': 'object',
            'aminoAcid': 'object',
            'copy': 'int',
            'frequency': 'float',
            'cdr3Length': 'int',
            'vFamilyName': 'object',
            'vGeneName': 'object',
            'vGeneAllele': 'float64',
            'vGeneNameTies': 'object',
            'dFamilyName': 'object',
            'dGeneName': 'object',
            'dGeneAllele': 'float64',
            'jFamilyName': 'object',
            'jGeneName': 'object',
            'jGeneAllele': 'float64',
            'd5Deletion': 'int64',
            'd3Deletion': 'int64',
            'jDeletion': 'int64',
            'n2Insertion': 'int64',
            'n1Insertion': 'int64',
            'vIndex': 'int64',
            'n2Index': 'int64',
            'dIndex': 'int64',
            'n1Index': 'int64',
            'jIndex': 'int64',
            'vAlignLength': 'int64',
            'vAlignSubstitutionCount': 'int64',
            'dAlignLength': 'int64',
            'dAlignSubstitutionCount': 'float64',
            'jAlignLength': 'int64',
            'jAlignSubstitutionCount': 'int64'
        }
        
        df = pd.read_csv(filepath, sep="\t", 
                         usecols=dtypes.keys(), 
                         dtype=dtypes, 
                         nrows=n_rows)
#         df['vSegment'] = df.apply(lambda row: self._slice_V_segement(row), axis=1)
        return df
    
    
    def load_dir(self, datadir, n_rows=None):
        """A function to load all files in a directory 
        
        args:
            : datadir (str): the directory of the files
            : n_rows (int or None): if int, will read only the first n rows
        """
        self.data_dir = datadir
        dirpath = f"{self.root_dir}{datadir}/"
        df_list = []
        for f in os.listdir(dirpath):
            filepath = f"{dirpath}{f}"
            df = self._read_file(filepath, n_rows)
            df['filepath'] = f
            df_list.append(df)
            
        self.data = pd.concat(df_list, ignore_index=True)
    
    
    def get_samples(self, sample_size=1000, n_samples=100, 
                    return_list=False):
        """a function to create a set of n_samples of size sample_size
        
        NOTE: depends on having a data loaded via load_dir

        args:
            : sample_size (int): the number of records in each sample
            : n_samples (int): the number of samples to create
            : return_list (bool): if true, returns the list

        returns:
            : samples (list of pd.DataFrame): a list of dataframes
        """
        samples = []
        for _ in range(n_samples):
            samples.append(self.data.sample(sample_size))
        self.samples = samples
        
        
    def samples_to_distance_matrices(self):
        """ A function to convert samples to distance matrices"""
    
        matrices = []
        pool = mp.Pool(36)
        self.matrices = pool.map(build_distance_matrix, self.samples)
            
        
            
            
            
            