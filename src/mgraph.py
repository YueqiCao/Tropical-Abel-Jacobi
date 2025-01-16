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
        Initialize the MetricGraph by calling the base class constructor.
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
        for u, v in nx.bfs_edges(ST, source=root):
            weight = self[u][v].get("weight")
            D.add_edge(u, v, weight=weight)
        
        # Direct the remaining edges in G\ST
        tree_edges = set(ST.edges())
        for u, v in self.edges():
            if (u, v) not in tree_edges and (v, u) not in tree_edges:
                # Use the relative depth in the BFS tree to decide the orientation
                if nx.shortest_path_length(ST, source=root, target=u) < nx.shortest_path_length(ST, source=root, target=v):
                    D.add_edge(u, v, weight = self[u][v].get("weight"))
                else:
                    D.add_edge(v, u, weight = self[u][v].get("weight"))
        
        return D
    
    def homology_basis(self, ST):
        '''
        compute the homology basis of the metric graph with respect to a given spanning tree
        each edge not in the spanning tree corresponds to a unique homology cycle
        '''

        # orient the graph with respect to the spanning tree
        D = self.tree_orientation(ST)

        # find edges not in the spanning tree
        tree_edges = set(ST.edges())
        non_tree_edges = [
            edge for edge in D.edges() 
            if (edge[0], edge[1]) not in tree_edges and (edge[1], edge[0]) not in tree_edges
        ]
        
        homology_basis = []

        for edge in non_tree_edges:

            # create the cycle graph
            cycle = nx.DiGraph()

            # add the non-tree edge
            u, v = edge
            weight = D[u][v]["weight"]
            cycle.add_edge(u, v, weight = weight, sign = 1)

            # find the path in the spanning tree between u and v
            path = nx.shortest_path(ST, source = u, target = v)

            # add edges along the path in the spanning tree
            # "sign" indicates whether the directed edge in D 
            # is consistent with the cycle orientation
            for i in range(len(path)-1):
                e = (path[i], path[i+1]) 
                if e in D.edges():
                    sign = -1
                    x, y = e[0], e[1]
                else:
                    sign = 1
                    x, y = e[1], e[0]
                cycle.add_edge(x, y, weight = self[x][y]["weight"], sign = sign)

            homology_basis.append(cycle)

        return homology_basis

    def refine(self, ratio=2, method="random"):
        '''
        Refine the combinatorial model by edge subdivision
        '''
        pass

    def simplify(self):
        '''
        Simplify the combinatorial model by edge contraction 
        Multigraph case not implemented yet
        '''
        pass
