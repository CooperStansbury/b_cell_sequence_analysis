import numpy as np
import pandas 
import scipy 
from scipy.sparse import csgraph
from scipy.stats import entropy
import networkx as nx



class MatUtil():
    """A matrix utility class """
    
    def __init__(self, M):
        """ input maxtrix """
        self.M = M
        
        
    def to_adj(self, tolerance=np.mean):
        """A function to covert distance matrix to adjacency matrix 

        args:
            : tolerance (str or func): how to covert the matrix to an
            adjacency. If float, all distances greater than or equal to 
            tolerance are coverted to 1. Any function that computes a floating
            point value from a matrix is acceptable otherwise.
        """ 
        if not callable(tolerance): 
            A = np.where(self.M >= tolerance, 1, 0)
        else:
            thresh = tolerance(self.M)
            A = np.where(self.M >= thresh, 1, 0)
        self.M = A.astype(int)
        
        
    def get_fiedler_number(self):
        """A function to return the fiedler number for the
        normed laplacian"""
        L = csgraph.laplacian(self.M, normed=True)
        w, _ = np.linalg.eigh(L, UPLO='L')
        return w[1] # second smallest item
        
    def get_von_Neumann_Entropy(self):
        """A function to compute the entropy of a normed 
        laplacian.
        
        NOTE: this follows from: the following:
        
            E(G) = - sum(F(lambda))
            
            where:
                G: symmetric graph
                F(x): continuous function, x * log(x) if x > 0, else 0
                lambda: eigenvalues of G
        """
        L = csgraph.laplacian(self.M, normed=True)
        w, _ = np.linalg.eigh(L, UPLO='L') 
        # w = w/w.sum()
        # S = -np.sum([x*np.log(x) for x in w if x > 0])
        w = w[w>0]
        return entropy(w)
        
    
        
        
        
        
        
        
        
        
        
        
        
        



class Networker():
    """A class to manage a few common metrics based on a distance matrices """
        
    def __init__(self, M):
        """ input maxtrix """
        
        self.matrix = M
        self.graph = None
    
    
    def to_adj(self, tolerance=np.mean):
        """A function to covert distance matrix to adjacency matrix 
        
        args:
            : tolerance (str or func): how to covert the matrix to an
            adjacency. If float, all distances greater than or equal to 
            tolerance are coverted to 1. Any function that computes a floating
            point value from a matrix is acceptable otherwise.
        """ 
        if not callable(tolerance): 
            A = np.where(self.matrix >= tolerance, 1, 0)
        else:
            thresh = tolerance(self.matrix)
            A = np.where(self.matrix >= thresh, 1, 0)
        self.matrix = A.astype(int)
            
    
    def to_graph(self, kind='adjacency'):
        """A function to covert to graph from array 
        
        args:
            : kind (str): either `adjacency` or `distance`
        """
        if kind == 'adjacency':      
            G = nx.from_numpy_matrix(self.matrix, parallel_edges=True)
        else:
            G = nx.from_numpy_matrix(self.matrix, parallel_edges=False)
            
        G.pos = nx.spring_layout(G)
        self.graph = G
            
            
    def _redux(self, node_list, func):
        """A utility function to reduce a set of node paramters by an 
        aggregation function
        
        args:
            : node_list (dict): list of nodes/node params
        
        returns:
            agg (float): a floating point number
        """
        return func(list(node_list.values()))
        
            
    def compute_centralities(self, redux_func=np.max):
        """A function to return a result dict for different cenrtality 
        metrics 
        
        NOTE: depends on having the graph object defined via to_graph.
        
        args:
            : redux_func (callable): an aggregation function
        
        returns:
            : centralities (dict): dictionary of results for a given matrix
        """
        
        if self.graph is None:
            raise ValueError("self.graph cannot be `None`. Call self.to_graph().")
        
        C = {}
        
        C['degree'] = self._redux(nx.degree_centrality(self.graph), redux_func)
        C['eigenvector'] = self._redux(nx.eigenvector_centrality(self.graph), redux_func)
        C['closeness'] = self._redux(nx.closeness_centrality(self.graph), redux_func) 
        C['information'] = self._redux(nx.information_centrality(self.graph), redux_func) 
        C['betweenness'] = self._redux(nx.betweenness_centrality(self.graph), redux_func) 
        C['harmonic'] = self._redux(nx.harmonic_centrality(self.graph), redux_func)     
        return C
        
        
        
        
            
        
        
    
    