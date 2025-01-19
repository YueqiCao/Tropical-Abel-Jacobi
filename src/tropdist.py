# author: Yueqi Cao
# date: 15/01/2025
# contact: y.cao21@imperial.ac.uk
#
# This script contains functions to compute tropical distances 
# on the tropical Jacobian. 

import numpy as np
import mgraph as mg
import networkx as nx
import fpylll as fp

def generate_random_graph(n, max_edges=None):
    """
    Generate a random connected graph with n nodes.
    max_edges: Maximum number of edges in the graph. If None, defaults to n.
    """
    
    # initialize the graph
    G = mg.MetricGraph()
    G.add_nodes_from(range(n))

    # create a line graph
    for i in range(n-1):
        # round the weight for better visualization
        weight = round(np.random.rand(), 2)  
        G.add_edge(i, i+1, weight=weight)
    
    # Randomly connect edges
    # Add additional random edges (up to max_edges if specified)
    max_edges = max_edges or n
    while G.number_of_edges() < max_edges:
        node1, node2 = np.random.choice(range(n), size=2, replace=False)
        if not G.has_edge(node1, node2):
            weight = round(np.random.rand(),2)
            G.add_edge(node1, node2, weight=weight)
    
    return G

def LLL_reduced_lattice(lattice, scale=100):
    '''
    Compute the reduced lattice basis using LLL algorithm.
    
    Parameters:
    lattice: np.array, the lattice basis (NB: numpy is not compatible with fpylll)
    scale: int, the scaling factor for the lattice basis (fpylll requires integer matrix)
    '''
    # fpylll requires IntegerMatrix type
    # in fpylll, the lattice vectors are row vectors
    lattice_adjusted = (scale * lattice.T).astype(int)
    A = fp.IntegerMatrix.from_matrix(lattice_adjusted)
    A_LLL = fp.LLL.reduction(A)
    lattice_reduced = np.array(list(A_LLL)).T/scale # transpose back
    return A_LLL, lattice_reduced

def draw_lattice_2D(lattice, ax):
    '''
    Draw the fundamental domain of a flat torus given by a matrix of lattice basis.
    '''
    u = lattice[:,0].reshape(1,-1)
    v = lattice[:,1].reshape(1,-1)

    l1 = np.vstack((-0.2*u, 1.2*u))
    l2 = np.vstack((-0.2*u+v, 1.2*u+v))
    l3 = np.vstack((-0.2*v, 1.2*v))
    l4 = np.vstack((-0.2*v+u, 1.2*v+u))
    ax.plot(l1[:,0], l1[:,1], 'k--',
            l2[:,0], l2[:,1], 'k--',
            l3[:,0], l3[:,1], 'k--', 
            l4[:,0], l4[:,1], 'k--')

def translate_to_fundamental_domain(points, lattice):
    '''
    Translate vectors to the fundamental domain for better visualization
    '''
    # compute the coordinates of points in the basis of the lattice
    X = np.linalg.inv(lattice) @ points
    # compute the translations
    N = np.floor(X)
    # translate the points to the fundamental domain
    points_new = points - lattice @ N
    return points_new