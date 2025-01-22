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

from tropdist import *
import numpy as np
import networkx as nx
import time
from g6k import Siever

# set random seed for reproducibility
np.random.seed(0)

from fpylll import IntegerMatrix, LLL, FPLLL

FPLLL.set_random_seed(0x1337)
A = IntegerMatrix.random(50, "qary", k=25, bits=20)
A = LLL.reduction(A)

g6k = Siever(A)
g6k.initialize_local(0, 0, 50)
g6k(alg="gauss")

i, norm, coeffs = g6k.best_lifts()[0]
l = int(round(norm))
print(l)
