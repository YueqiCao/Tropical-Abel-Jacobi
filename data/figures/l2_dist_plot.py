# author: Yueqi Cao
# date: 22/01/2025
# contact: y.cao21@imperial.ac.uk
#
# plot for l2 distance

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import sys
from pathlib import Path
src_path = Path("../data") 
sys.path.append(str(src_path))

# load the results
df1 = pd.read_csv("../outputs/fplll_fixg.csv")
df2 = pd.read_csv("../outputs/g6k_fixg.csv")
df = pd.concat([df1, df2], ignore_index=True)
df["lognodes"] = np.log(df["Graph Nodes"])
df["logtime"] = np.log(df["Time"])
df["Method"] = df["Method"].replace("fplll", "FPLLL enumeration")
df["Method"] = df["Method"].replace("nv", "Nguyen-Vidick sieve")
df["Method"] = df["Method"].replace("gauss", "Gauss sieve")

# Plot with Seaborn
plt.figure(figsize=(6, 6))
sns.lineplot(
    data=df,
    x="lognodes",
    y="logtime",
    hue = "Method",
    style = "Method",
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

plt.savefig("../figures/l2_dist_fixg.pdf")

#####################################################################################

df1 = pd.read_csv("../outputs/fplll_fixn.csv")
df2 = pd.read_csv("../outputs/g6k_fixn.csv")
df = pd.concat([df1, df2], ignore_index=True)
df["loggenus"] = np.log(df["Graph Genus"])
df["logtime"] = np.log(df["Time"])
df["Method"] = df["Method"].replace("fplll", "FPLLL enumeration")
df["Method"] = df["Method"].replace("nv", "Nguyen-Vidick sieve")
df["Method"] = df["Method"].replace("gauss", "Gauss sieve")

# Plot with Seaborn
plt.figure(figsize=(6, 6))
sns.lineplot(
    data=df,
    x="loggenus",
    y="logtime",
    hue = "Method",
    style = "Method",
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

plt.savefig("../figures/l2_dist_fixn.pdf")