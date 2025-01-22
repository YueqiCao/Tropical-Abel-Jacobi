# author: Yueqi Cao
# date: 22/01/2025
# contact: y.cao21@imperial.ac.uk
#
# pgf plot for fpylll_dist.py

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import sys
from pathlib import Path
src_path = Path("../data") 
sys.path.append(str(src_path))

# load the results
data = np.load("../outputs/full_l2_dist.npz")

# parameters
# set hyperparameters
num_exp = 10 # number of experiments in one loop
num_loops = 15 # number of loops

# fixed genus
g, initial_nodes, step_nodes = 15, 20, 10
node_list = [n for n in range(initial_nodes, initial_nodes + num_loops * step_nodes, step_nodes)]

time_fixg_fplll = data["fixg_fplll"] # more data to come

# Convert data to a long-format DataFrame for Seaborn
pd_data = []
for num_nodes, time_values in zip(node_list, time_fixg_fplll):
    for time in time_values:
        pd_data.append({"Graph Nodes": np.log(num_nodes), "Time": np.log(time), "lib": "fplll"})

# more data to come

# covert to pd
df = pd.DataFrame(pd_data)

# Plot with Seaborn
plt.figure(figsize=(6, 6))
sns.lineplot(
    data=df,
    x="Graph Nodes",
    y="Time",
    hue = "lib",
    style = "lib",
    errorbar=("ci", 95), 
    marker="o",
    linewidth=2,
    color="royalblue"
)

# Customise the plot
plt.title("Time complexity of computing the full tropical polarization matrix", fontsize=14)
plt.xlabel("log(Graph Nodes)", fontsize=12)
plt.ylabel("log(Time)", fontsize=12)
plt.legend(fontsize=12, title_fontsize=14)
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

plt.savefig("../figures/full_dist_fixg.pdf")

#####################################################################################

# fixed number of nodes
num_nodes, initial_g, step_g = 100, 5, 4
g_list = [g for g in range(initial_g, initial_g + num_loops * step_g, step_g)]
time_fixn_fplll = data["fixn_fplll"]


# Convert data to a long-format DataFrame for Seaborn
pd_data = []
for g, time_values in zip(g_list, time_fixn_fplll):
    for time in time_values:
        pd_data.append({"Graph Genus": np.log(g), "Time": np.log(time), "lib": "fplll"})
df = pd.DataFrame(pd_data)

# Plot with Seaborn
plt.figure(figsize=(6, 6))
sns.lineplot(
    data=df,
    x="Graph Genus",
    y="Time",
    hue = "lib",
    style = "lib",
    errorbar=("ci", 95), 
    marker="o",
    linewidth=2,
    color="royalblue"
)

# Customise the plot
plt.title("Time complexity of computing the full tropical polarization matrix", fontsize=14)
plt.xlabel("log(Genus)", fontsize=12)
plt.ylabel("log(Time)", fontsize=12)
plt.legend(fontsize=12, title_fontsize=14)
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

plt.savefig("../figures/full_dist_fixn.pdf")

