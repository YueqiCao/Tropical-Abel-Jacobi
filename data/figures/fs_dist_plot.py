# author: Yueqi Cao
# date: 22/01/2025
# contact: y.cao21@imperial.ac.uk
#
# plot for Foster-Zhang distance

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import sys
from pathlib import Path
src_path = Path("../data") 
sys.path.append(str(src_path))

# load the results
df = pd.read_csv("../outputs/fs_fixg.csv")
df["lognodes"] = np.log(df["Graph Nodes"])
df["logtime"] = np.log(df["Time"])
df["Solver"] = df["Solver"].replace("cbc", "CBC")
df["Solver"] = df["Solver"].replace("ipopt", "IPOPT")
df["Solver"] = df["Solver"].replace("glpk", "GLPK")
df["Solver"] = df["Solver"].replace("scip", "SCIP")

# Plot with Seaborn
plt.figure(figsize=(6, 6))
sns.lineplot(
    data=df,
    x="lognodes",
    y="logtime",
    hue = "Solver",
    style = "Solver",
    errorbar=("ci", 95), 
    marker="o",
    linewidth=2,
    color="royalblue"
)

# Customise the plot
#plt.title("Time complexity of computing the full tropical polarization matrix", fontsize=14)
plt.xlabel("log(#Nodes)", fontsize=12)
plt.ylabel("log(Time)", fontsize=12)
plt.legend(fontsize=12, title_fontsize=14)
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

plt.savefig("../figures/fs_dist_fixg.pdf")

#####################################################################################

df = pd.read_csv("../outputs/fs_fixn.csv")
df["loggenus"] = np.log(df["Graph Genus"])
df["logtime"] = np.log(df["Time"])
df["Solver"] = df["Solver"].replace("cbc", "CBC")
df["Solver"] = df["Solver"].replace("ipopt", "IPOPT")
df["Solver"] = df["Solver"].replace("glpk", "GLPK")
df["Solver"] = df["Solver"].replace("scip", "SCIP")

# Plot with Seaborn
plt.figure(figsize=(6, 6))
sns.lineplot(
    data=df,
    x="loggenus",
    y="logtime",
    hue = "Solver",
    style = "Solver",
    errorbar=("ci", 95), 
    marker="o",
    linewidth=2,
    color="royalblue"
)

# Customise the plot
#plt.title("Time complexity of computing the full tropical polarization matrix", fontsize=14)
plt.xlabel("log(Genus)", fontsize=12)
plt.ylabel("log(Time)", fontsize=12)
plt.legend(fontsize=12, title_fontsize=14)
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

plt.savefig("../figures/fs_dist_fixn.pdf")