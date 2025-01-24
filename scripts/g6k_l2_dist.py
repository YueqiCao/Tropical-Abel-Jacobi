# author: Yueqi Cao
# date: 22/01/2025
# contact: y.cao21@imperial.ac.uk
#
# This script is used to test the complexity of computing the full, exact tropical polarization matrix
# using the sieving algorithm in G6K

import sys
from pathlib import Path
src_path = Path("../src") 
sys.path.append(str(src_path))
sys.path.append('~/g6k')

from tropdist import *
import numpy as np
import pandas as pd
import networkx as nx
import time
import fpylll as fp
from g6k import Siever
import logging
logging.basicConfig(level=logging.ERROR)

# set random seed for reproducibility
np.random.seed(0)

# function to compute the tropical polarization distance matrix
def g6k_L2_distance(V_transformed, Q_sqrt, scale=1e5, method="gauss"):

    # extend the lattice matrix
    Q_ext = np.hstack((Q_sqrt, np.zeros((Q_sqrt.shape[0], 1))))
    Q_scaled = np.round(Q_ext * scale).astype("int")

    # compute the distance matrix
    g, n = V_transformed.shape
    D = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1, n):
            # embed the lattice to a higher dimension
            target = np.round((V_transformed[:,i] - V_transformed[:,j])*scale).astype("int")
            target_ext = np.append(target, 1)
            Q_new = np.vstack((Q_scaled, target_ext.reshape(1,-1)))
            # compute the closest vector using Kannan's embedding        
            A = fp.IntegerMatrix.from_matrix(Q_new)
            #A = fp.LLL.reduction(A)
            # Run the sieve to solve CVP
            siever = Siever(A)
            siever.initialize_local(0, 0, g+1)
            siever(alg=method)
            best_lifts = siever.best_lifts()
            if best_lifts:
                i, norm_pre, coeffs = best_lifts[0] # shortest vector sieved
                norm = int(round(norm_pre))
            else:
                norm = min([sum(v**2 for v in vector) for vector in A])
            # "norm" is the squared norm of the shortest vector
            D[i,j] = np.sqrt(norm - 1)/scale
            D[j,i] = D[i,j]
    
    return D

# set hyperparameters
num_exp = 6 # number of experiments in one loop
num_loops = 10 # number of loops
methods = ["gauss", "nv"] # other advanced sieve methods for large scale "hk3", "bgj1", "bdgl", "bdgl1", "bdgl2", "bdgl3"

# fixed genus
g, initial_nodes, step_nodes = 15, 20, 10

# loop over methods
fixg = []
for method in methods:
    # loop over different number of nodes
    for i in range(num_loops):
        # loop over different experiments
        for j in range(num_exp):
            # generate a random graph
            num_nodes = initial_nodes + i * step_nodes
            MG = generate_random_graph(num_nodes, max_edges=num_nodes+g-1)
            # compute the tropical Abel--Jacobi transform
            MST = nx.minimum_spanning_tree(MG)
            base_point = list(MST.nodes())[0]
            V0 = MG.trop_transform(MST, base_point)
            # sample more points from the transform and drop the repeated points
            V1 = np.unique(MG.interpolate(MST, V0, 5), axis=1)
            # fix the number of points to be computed
            V = V1[:, :num_nodes]
            Q = MG.trop_polarization(MST)
            V_L2, Q_sqrt = to_L2(V, Q)
            # record time
            start_time = time.time()
            dist_g6k = g6k_L2_distance(V_L2, Q_sqrt, method=method)
            end_time = time.time()
            fixg.append({"Graph Nodes": num_nodes, "Time": end_time-start_time, "Method": method})
            
# save DataFrame
df = pd.DataFrame(fixg)
df.to_csv("../data/outputs/g6k_fixg.csv", index=False)

# fixed number of nodes
num_nodes, initial_g, step_g = 50, 5, 4

#loop over methods
fixn = []
for method in methods:
    # loop over different genus
    for i in range(num_loops):
        # loop over different experiments
        for j in range(num_exp):
            # generate a random graph
            g = initial_g + i * step_g
            MG = generate_random_graph(num_nodes, max_edges=num_nodes+g-1)
            # compute the tropical Abel--Jacobi transform
            MST = nx.minimum_spanning_tree(MG)
            base_point = list(MST.nodes())[0]
            V = MG.trop_transform(MST, base_point)
            Q = MG.trop_polarization(MST)
            V_L2, Q_sqrt = to_L2(V, Q)
            # record time
            start_time = time.time()
            dist_g6k = g6k_L2_distance(V_L2, Q_sqrt, method=method)
            end_time = time.time()
            fixn.append({"Graph Genus": g, "Time": end_time-start_time, "Method": method})

# save DataFrame
df = pd.DataFrame(fixn)
df.to_csv("../data/outputs/g6k_fixn.csv", index=False)