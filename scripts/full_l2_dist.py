# author: Yueqi Cao
# date: 22/01/2025
# contact: y.cao21@imperial.ac.uk
#
# This script is used to test the complexity of computing the full, exact tropical polarization matrix
# using the enumeration algorithm in fpylll

import sys
from pathlib import Path
src_path = Path("../src") 
sys.path.append(str(src_path))

from tropdist import *
import numpy as np
import networkx as nx
import time

# set random seed for reproducibility
np.random.seed(0)

# set hyperparameters
num_exp = 10 # number of experiments in one loop
num_loops = 15 # number of loops

# fixed genus
g, initial_nodes, step_nodes = 15, 20, 10

time_fixg_fplll = np.zeros((num_loops, num_exp))

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
        # compute the full tropical polarization matrix
        start_time = time.time()
        dist_exact = exact_L2_distance(V_L2, Q_sqrt)
        end_time = time.time()
        time_fixg_fplll[i, j] = end_time - start_time

# fixed number of nodes
num_nodes, initial_g, step_g = 100, 5, 4
time_fixn_fplll = np.zeros((num_loops, num_exp))

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
        # compute the full tropical polarization matrix
        start_time = time.time()
        dist_exact = exact_L2_distance(V_L2, Q_sqrt)
        end_time = time.time()
        time_fixn_fplll[i, j] = end_time - start_time

# save the time records
np.savez("../data/outputs/full_l2_dist.npz", fixg_fplll=time_fixg_fplll, fixn_fplll=time_fixn_fplll)