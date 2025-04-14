# TropAJ: Tropical Abel--Jacobi Transform of Metric Graphs

this repository contains code for computation of tropical Abel--Jacobi transform of metric graphs.

You need to use WSL to install ``fpylll``.

To compute the Foster--Zhang distance matrix. There is no need to install ``fpylll`` (no lattice reduction)
but need to install ``pyomo`` and any solver (``coincbc``, ``gurobi``, ``cplex`` etc.).

## Usage

The repository is organized as follows:

- ```src``` contains modules for tropical transform and tropical distances:
    - ```mgraph.py```: metric graph class, with functions for operations on metric graphs and computations of the tropical Abel--Jacobi transform
    - ```tropdist```: functions for computations of tropical distances
- ```notebooks``` contains jupyter notebook tutorials:
    - ```tropical transform```: tutorial of tropical Abel--Jacobi transform
    - ```tropical distances```: tutorial of tropical distances on the tropical Jacobian
- ```scripts``` contains source codes for experiments in the original paper:
    - ```fplll_l2_dist.py```: computation using ```fplll```
    - ```g6k_l2_dist.py```: computation using ```g6k```
    - ```fs_dist.py```: computation using MIP solvers
    - ```babai_dist_time.py``` and ```babai_dist_error.py```: computation using Babai's algorithms
- ```data``` contains numerical outputs of experiments in the original paper


## Academic Use
You can use the following BibTex entry:
```
arXiv
```


## Requirements

fplll

## Contact

If you are interested in this work, or if you find any bug in the code, please contact y.cao21@imperial.ac.uk