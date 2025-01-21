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
    The metric graph class is a subclass of NetworkX.Graph.
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
            # add edges in the spanning tree
            D.add_edge(u, v, weight=weight)
        
        # Direct the remaining edges in G\ST
        tree_edges = set(ST.edges())
        for u, v in self.edges():
            if (u, v) not in tree_edges and (v, u) not in tree_edges:
                # Use the relative depth in the BFS tree to decide the orientation
                # add edges not in the spanning tree
                if nx.shortest_path_length(ST, source=root, target=u) < nx.shortest_path_length(ST, source=root, target=v):
                    D.add_edge(u, v, weight = self[u][v].get("weight"))
                else:
                    D.add_edge(v, u, weight = self[u][v].get("weight"))
        
        return D
    
    def homology_basis(self, ST):
        '''
        Compute the homology basis of the metric graph with respect to a given spanning tree.
        Each edge not in the spanning tree corresponds to a unique homology cycle.
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
    
    def cycle_edge(self, ST):
        '''
        Compute the cycle-edge incidence matrix with respect to a given spanning tree.
        The (i,j)th element in the matrix indicates the coefficient of edge j in the
        ith homology cycle. 
        '''
        
        # orient the graph with respect to the spanning tree
        D = self.tree_orientation(ST)

        # compute the homology basis
        homology_basis = self.homology_basis(ST)

        # compute the cycle-edge incidence matrix
        g = len(homology_basis)
        edge_list = list(D.edges())
        m = len(edge_list)
        C = np.zeros((g,m))
        for i in range(g):
            cycle = homology_basis[i].edges(data=True)
            for edge in cycle:
                index = edge_list.index((edge[0], edge[1]))
                C[i, index] = edge[2]["sign"]

        return C

    def reduced_cycle_edge(self, ST):
        '''
        Compute the reduced cycle-edge incidence matrix with respect to a given spanning tree.
        The reduced version only retains the edges in the spanning tree
        '''

        # compute the homology basis
        homology_basis = self.homology_basis(ST)

        # compute the cycle-edge incidence matrix
        g = len(homology_basis)
        edge_list = list(ST.edges())
        n_minus = len(edge_list)
        C = np.zeros((g,n_minus))
        for i in range(g):
            cycle = homology_basis[i].edges(data=True)
            for edge in list(cycle)[1:]: # by construction the first edge is not in ST
                try:
                    index = edge_list.index((edge[0], edge[1]))
                except ValueError:
                    index = edge_list.index((edge[1], edge[0]))
                C[i, index] = edge[2]["sign"]

        return C

    def path_edge(self, ST, base_point):
        '''
        Compute the path-edge incidence matrix with respect to a given spanning tree and a base point.
        The (i,j)th element in the matrix indicates the coefficient of edge j in the path from the base
        point to the ith vertex. 
        '''
        
        n = len(self.nodes())
        m = len(self.edges())
        Y = np.zeros((n,m))

        # orient the graph with respect to the spanning tree
        D = self.tree_orientation(ST)

        # compute the path-edge incidence matrix
        for i, u in enumerate(list(self.nodes())):
            path = nx.shortest_path(ST, source=base_point, target=u)
            if len(path)>0:
                for j in range(len(path)-1):
                    try:
                        index = list(D.edges()).index((path[j],path[j+1]))
                        Y[i, index] = 1
                    except ValueError:
                        index = list(D.edges()).index((path[j+1],path[j]))
                        Y[i, index] = -1

        return Y

    def reduced_path_edge(self, ST, base_point):
        '''
        Compute the reduced path-edge incidence matrix with respect to a given spanning tree
        and a base point. The reduced version only retains the edges in the spanning tree.
        '''

        n = len(self.nodes())
        Y = np.zeros((n,n-1))

        # orient the graph with respect to the spanning tree
        D = self.tree_orientation(ST)
        # get oriented edges for the spanning tree
        ori_ST = []
        for edge in ST.edges():
            if edge in D.edges():
                ori_ST.append(edge)
            else:
                ori_ST.append((edge[1],edge[0]))

        # compute the path-edge incidence matrix
        for i, u in enumerate(list(self.nodes())):
            path = nx.shortest_path(ST, source=base_point, target=u)
            if len(path)>0:
                for j in range(len(path)-1):
                    try:
                        index = ori_ST.index((path[j],path[j+1]))
                        Y[i, index] = 1
                    except ValueError:
                        index = ori_ST.index((path[j+1],path[j]))
                        Y[i, index] = -1

        return Y
    
    def edge_length(self, ST):
        '''
        Compute the edge length matrix
        '''

        # to ensure the edge length is sorted in a consistent way with 
        # previous constructions of C and Y
        D = self.tree_orientation(ST)

        weight_list = []
        for edge in D.edges(data=True):
            weight_list.append(edge[2]["weight"])
        
        return np.diag(np.array(weight_list))
        

    def reduced_edge_length(self, ST):
        '''
        Compute the reduced edge lenth matrix with respect to a spanning tree.
        '''
        
        weight_list = []
        for edge in ST.edges(data=True):
            weight_list.append(edge[2]["weight"])

        return np.diag(np.array(weight_list))

    def trop_transform(self, ST, base_point):
        '''
        Compute the tropical Abel--Jacobi transform with respect to a spanning tree
        and a base point.
        '''
        
        C_st = self.reduced_cycle_edge(ST)
        L_st = self.reduced_edge_length(ST)
        Y_st = self.reduced_path_edge(ST, base_point)

        return C_st @ L_st @ Y_st.T
    
    def trop_polarization(self, ST):
        '''
        Compute the tropical polarization matrix, which is also the matrix of lattice basis 
        '''

        # construct edge length matric for edges not in the spanning tree
        L_g_list = []
        homology_basis = self.homology_basis(ST)
        for cycle in homology_basis:
            edges = list(cycle.edges(data=True))
            L_g_list.append(edges[0][2]["weight"])

        L_g = np.diag(np.array(L_g_list))
        
        # compute the tropical polarization matrix
        C_st = self.reduced_cycle_edge(ST)
        L_st = self.reduced_edge_length(ST)
        Q = C_st @ L_st @ C_st.T + L_g 

        return Q
    
    def interpolate(self, ST, V, ratio, method="equidistant"):
        '''
        Sample points from the tropical Abel--Jacobi transform by interpolation
        ''' 

        # orient the graph with respect to the spanning tree
        D = self.tree_orientation(ST)

        # get oriented edges for the spanning tree
        ori_ST = []
        # exclude bridges
        bridges = list(nx.bridges(self))
        for edge in ST.edges():
            # if the edge is a bridge, skip
            if edge in bridges:
                continue
            # otherwise, add it
            if edge in D.edges():
                ori_ST.append(edge)
            else:
                ori_ST.append((edge[1],edge[0]))

        # new matrix for the interpolated points        
        V_new = V

        # interpolate points in edges in the spanning tree
        for edge in ori_ST:
            u, v = edge
            j_ini = list(self.nodes()).index(u)
            x = V[:,j_ini].reshape(-1,1)
            j_ter = list(self.nodes()).index(v)
            y = V[:,j_ter].reshape(-1,1)
            if method == "equidistant":
                for i in range(1, ratio+1):
                    V_new = np.hstack((V_new, (1-i/(ratio+1))*x + i/(ratio+1)*y))
            elif method == "random":
                for i in range(1, ratio+1):
                    theta = np.random.rand()
                    V_new = np.hstack((V_new, theta*x + (1-theta)*y))
        
        # interpolate points in edges not in the spanning tree
        homology_basis = self.homology_basis(ST)
        g = len(homology_basis)
        for i in range(g):
            cycle = homology_basis[i]
            edge = list(cycle.edges(data=True))[0]
            u = edge[0]
            w = np.zeros((g,1))
            w[i] = edge[2]["weight"]
            j_ini = list(self.nodes()).index(u)
            x = V[:,j_ini].reshape(-1,1)
            if method == "equidistant":
                for j in range(1, ratio+1):
                    V_new = np.hstack((V_new, x + j/(ratio+1)*w))
            elif method == "random":
                for j in range(1, ratio+1):
                    theta = np.random.rand()
                    V_new = np.hstack((V_new, x + theta*w))
        
        return V_new

    def edge_subdivision(self, ratio=1):
        '''
        Refine the combinatorial model by edge subdivision
        '''

        # Create a list to hold new edges and nodes
        new_nodes = []
        new_edges = []
        remove_edges = []
        
        # Process each edge in the graph
        for edge in self.edges(data=True):  
            u, v = edge[0], edge[1]
            weight = edge[2]["weight"]  
            
            # Compute the weight of each sub-edge
            sub_edge_weight = weight / (ratio + 1)
            
            # Create new vertices along the edge
            interior_vertices = [f"{u}-{v}-{i}" for i in range(1, ratio + 1)]
            
            # Add the new vertices and edges
            previous_vertex = u
            for new_vertex in interior_vertices:
                new_nodes.append(new_vertex)  # Add the new vertex
                new_edges.append((previous_vertex, new_vertex, sub_edge_weight))  # Add edge to new vertex
                previous_vertex = new_vertex
            
            # Add the final edge to v
            new_edges.append((previous_vertex, v, sub_edge_weight))
            
            # Remove the original edge
            remove_edges.append(edge)
        
        # Add new vertices and edges to the graph
        self.add_nodes_from(new_nodes)
        self.add_weighted_edges_from(new_edges)
        self.remove_edges_from(remove_edges)

    def remove_bridges(self):
        '''
        Contract bridges from the graph
        '''

        bridges = list(nx.bridges(self))

        while bridges != []:
            for u, v in bridges:
                self.remove_edge(u,v)
                # Contract the bridge edge by merging u and v into a single node
                # Map all neighbours of v to u, preserving edge weights
                for neighbour in list(self.neighbors(v)):
                    self.add_edge(u, neighbour, weight=self[v][neighbour]["weight"])
                    self.remove_edge(v, neighbour)
                    if (v, neighbour) in bridges:
                        bridges.remove((v, neighbour))
                self.remove_node(v)
            bridges = list(nx.bridges(self))

    def remove_interior_nodes(self):
        '''
        Simplify the combinatorial model by remove interior nodes
        '''

        # list of vertices to process
        remove_nodes = [node for node in self.nodes() if self.degree[node] == 2]

        for node in remove_nodes:
            # skip nodes with self-loops (special case)
            if self.has_edge(node, node):
                continue
            
            # get the neighbours of the node
            neighbours = list(self.neighbors(node))
            if len(neighbours) != 2:
                continue
            
            # extract the two neighbours and the weights of the edges
            u, v = neighbours
            weight1 = self[u][node]["weight"]
            weight2 = self[node][v]["weight"]
            
            if self.has_edge(u, v):
                continue
            else:
                # Remove the node and its edges
                self.remove_node(node)
                
                # Add a new edge between the two neighbours, summing the weights
                new_weight = weight1 + weight2
                self.add_edge(u, v, weight=new_weight)