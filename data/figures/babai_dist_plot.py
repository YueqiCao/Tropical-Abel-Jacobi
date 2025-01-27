# author: Yueqi Cao
# date: 22/01/2025
# contact: y.cao21@imperial.ac.uk
#
# plot for babai's method

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import sys
from pathlib import Path
src_path = Path("../data") 
sys.path.append(str(src_path))

# load the computational time for fixed genus
df = pd.read_csv("../outputs/babai_fixg.csv")
df["lognodes"] = np.log(df["Graph Nodes"])
df["logtime"] = np.log(df["Time"])
df["Method"] = df["Method"].replace("rounding", "Babai's rounding algorithm")
df["Method"] = df["Method"].replace("nearest", "Babai's nearest plane algorithm")

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
plt.xlabel("log(#Nodes)", fontsize=12)
plt.ylabel("log(Time)", fontsize=12)
plt.legend(fontsize=12, title_fontsize=14)
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

plt.savefig("../figures/babai_fixg.pdf")

#####################################################################################

# load the approximation error for fixed genus
df = pd.read_csv("../outputs/babai_errorg.csv")
df["lognodes"] = np.log(df["Graph Nodes"])
df["logmse"] = np.log(df["MSE"])
df["Method"] = df["Method"].replace("rounding", "Babai's rounding algorithm")
df["Method"] = df["Method"].replace("nearest", "Babai's nearest plane algorithm")

# Plot with Seaborn
plt.figure(figsize=(6, 6))
sns.lineplot(
    data=df,
    x="Graph Nodes",
    y="MSE",
    hue = "Method",
    style = "Method",
    errorbar=("ci", 95), 
    marker="o",
    linewidth=2,
    color="royalblue"
)

# Customise the plot
plt.xlabel("#Nodes", fontsize=12)
plt.ylabel("MSE", fontsize=12)
plt.legend(fontsize=12, title_fontsize=14)
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

plt.savefig("../figures/babai_errorg.pdf")


#####################################################################################

# load the computational time for fixed number of nodes
df = pd.read_csv("../outputs/babai_fixn.csv")
df["loggenus"] = np.log(df["Graph Genus"])
df["logtime"] = np.log(df["Time"])
df["Method"] = df["Method"].replace("rounding", "Babai's rounding algorithm")
df["Method"] = df["Method"].replace("nearest", "Babai's nearest plane algorithm")

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
plt.xlabel("log(Genus)", fontsize=12)
plt.ylabel("log(Time)", fontsize=12)
plt.legend(fontsize=12, title_fontsize=14)
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

plt.savefig("../figures/babai_fixn.pdf")

#####################################################################################

# load the approximation error for fixed number of nodes
df = pd.read_csv("../outputs/babai_errorn.csv")
df["loggenus"] = np.log(df["Graph Genus"])
df["logmse"] = np.log(df["MSE"])
df["Method"] = df["Method"].replace("rounding", "Babai's rounding algorithm")
df["Method"] = df["Method"].replace("nearest", "Babai's nearest plane algorithm")

# Plot with Seaborn
plt.figure(figsize=(6, 6))
sns.lineplot(
    data=df,
    x="Graph Genus",
    y="MSE",
    hue = "Method",
    style = "Method",
    errorbar=("ci", 95), 
    marker="o",
    linewidth=2,
    color="royalblue"
)

# Customise the plot
plt.xlabel("Genus", fontsize=12)
plt.ylabel("MSE", fontsize=12)
plt.legend(fontsize=12, loc='upper right', title_fontsize=14)
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

plt.savefig("../figures/babai_errorn.pdf")