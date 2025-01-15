# author: Yueqi Cao
# date: 15/01/2025
# contact: y.cao21@imperial.ac.uk
#
# This script defines a customized graph class which extends the 
# functionality of a NetworkX graph and includes additional methods
# to compute the tropical Abel--Jacobi transform.

import networkx as nx
import numpy as np

class MetricGraph(nx.Graph):
    '''
    the metric graph class is a subclass of NetworkX.Graph
    '''

    def __init__(self, *args, **kwargs):
        '''
        Initialize the CustomGraph by calling the base class constructor.
        '''
        super().__init__(*args, **kwargs)

    def tree_orientation(self, ST):
        '''
        Orient the metric graph with respect to a given spanning tree ST.

        The edges in the spanning tree are directed away from the root, and the
        remaining edges are oriented consistently based on the spanning tree structure.

        Return a directed version of the graph with canonical edge orientations.
        '''

        if not nx.is_connected(self):
            raise ValueError("The graph must be connected to compute a spanning tree.")
        
        root = list(ST.nodes())[0]

        # Create a directed graph
        D = nx.DiGraph()
        D.add_nodes_from(self.nodes())
        
        # Direct the spanning tree edges away from the root
        for edge in nx.bfs_edges(ST, source=root):
            u, v = edge
            D.add_edge(u, v)
        
        # Direct the remaining edges in G\ST
        tree_edges = set(ST.edges())
        for u, v in self.edges():
            if (u, v) not in tree_edges and (v, u) not in tree_edges:
                # Use the relative depth in the BFS tree to decide the orientation
                if nx.shortest_path_length(ST, source=root, target=u) < nx.shortest_path_length(ST, source=root, target=v):
                    D.add_edge(u, v)
                else:
                    D.add_edge(v, u)
        
        return D

    def refine(self, ratio=2, method="random"):
        '''
        Refine the combinatorial model by edge subdivision
        '''
        pass

    def simplify(self):
        '''
        Simplify the combinatorial model by edge contraction 
        '''
        pass
