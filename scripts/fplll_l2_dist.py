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
import pandas as pd
import networkx as nx
import time

# set random seed for reproducibility
np.random.seed(0)

# set hyperparameters
num_exp = 6 # number of experiments in one loop
num_loops = 10 # number of loops

# fixed genus
g, initial_nodes, step_nodes = 15, 20, 10

# loop over different number of nodes
fixg = []
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
        fixg.append({"Graph Nodes": num_nodes, "Time": end_time-start_time, "Method": "fplll"})
            
# save DataFrame
df = pd.DataFrame(fixg)
df.to_csv("../data/outputs/fplll_fixg.csv", index=False)

# fixed number of nodes
num_nodes, initial_g, step_g = 50, 5, 4

# loop over different genus
fixn = []
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
        fixn.append({"Graph Genus": g, "Time": end_time-start_time, "Method": "fplll"})

# save DataFrame
df = pd.DataFrame(fixn)
df.to_csv("../data/outputs/fplll_fixn.csv", index=False)