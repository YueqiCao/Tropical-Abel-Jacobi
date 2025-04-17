# TropAJ: Tropical Abel--Jacobi Transform of Metric Graphs

``TropAJ`` is a Python repository for computing the tropical Abel--Jacobi transform of a metric graph and tropical distances on its tropical Jacobian. This repository also contains source codes for experiments in the paper

> [Computing the Tropical Abel--Jacobi Transform and Tropical Distances for Metric Graphs](https://arxiv.org/abs/2504.11619)

## Tutorial

Some tutorials about the tropical Abel--Jacobi transform and this repository are included in ``notebooks``.   

## Package Dependence

The repository depends on several Python libraries from different scientific fields. For your purpose it is not necessary to install all of them.

- The computation of tropical Abel--Jacobi transform depends on [NetworkX](https://networkx.org/); 
- The computation of tropical polarization distance depends on [fplll](https://github.com/fplll/fplll) and [G6K](https://github.com/fplll/g6k) (not available on Windows system);
- The computation of Foster--Zhang distance depends on [Pyomo](https://www.pyomo.org/), and any MIP solver (in the paper, [Cbc](https://github.com/coin-or/Cbc), [Ipopt](https://coin-or.github.io/Ipopt/), [GLPK](https://www.gnu.org/software/glpk/), [SCIP](https://www.scipopt.org/) are tested).

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
@misc{cao2025computingtropicalabeljacobitransform,
      title={Computing the Tropical Abel--Jacobi Transform and Tropical Distances for Metric Graphs}, 
      author={Yueqi Cao and Anthea Monod},
      year={2025},
      eprint={2504.11619},
      archivePrefix={arXiv},
      primaryClass={math.AG},
      url={https://arxiv.org/abs/2504.11619}, 
}
```

## Contact

If you are interested in this work, or if you find any bug in the code, please contact y.cao21@imperial.ac.uk