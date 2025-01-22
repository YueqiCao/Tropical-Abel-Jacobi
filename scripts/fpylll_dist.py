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

time_recodes_fixg = np.zeros((num_loops, num_exp))

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
        V = MG.trop_transform(MST, base_point)
        Q = MG.trop_polarization(MST)
        V_L2, Q_sqrt = to_L2(V, Q)
        # compute the full tropical polarization matrix
        start_time = time.time()
        dist_exact = exact_L2_distance(V_L2, Q_sqrt)
        end_time = time.time()
        time_recodes_fixg[i, j] = end_time - start_time

# fixed number of nodes
num_nodes, initial_g, step_g = 100, 5, 4
time_recodes_fixn = np.zeros((num_loops, num_exp))

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
        time_recodes_fixn[i, j] = end_time - start_time

# save the time records
np.savez("../data/outputs/fpylll_dist.npz", fixg=time_recodes_fixg, fixn=time_recodes_fixn)