# author: Yueqi Cao
# date: 15/01/2025
# contact: y.cao21@imperial.ac.uk
#
# This script contains functions to compute tropical distances 
# on the tropical Jacobian. 

import numpy as np
import mgraph as mg
import networkx as nx
import fpylll as fp
import pyomo.environ as pyo

def generate_random_graph(n, max_edges=None):
    """
    Generate a random connected graph with n nodes.
    max_edges: Maximum number of edges in the graph. If None, defaults to n.
    """
    
    # initialize the graph
    G = mg.MetricGraph()
    G.add_nodes_from(range(n))

    # create a line graph
    for i in range(n-1):
        # round the weight for better visualization
        weight = round(np.random.rand(), 2)  
        G.add_edge(i, i+1, weight=weight)
    
    # Randomly connect edges
    # Add additional random edges (up to max_edges if specified)
    max_edges = max_edges or n
    while G.number_of_edges() < max_edges:
        node1, node2 = np.random.choice(range(n), size=2, replace=False)
        if not G.has_edge(node1, node2):
            weight = round(np.random.rand(),2)
            G.add_edge(node1, node2, weight=weight)
    
    return G

def LLL_reduced_lattice(lattice, scale=1e6):
    '''
    Compute the reduced lattice basis using LLL algorithm.
    
    Parameters:
    lattice: np.array, the lattice basis (NB: numpy is not compatible with fpylll)
    scale: int, the scaling factor for the lattice basis (fpylll requires integer matrix)
    '''
    # fpylll requires IntegerMatrix type
    # in fpylll, the lattice vectors are row vectors
    lattice_adjusted = np.round(scale * lattice.T).astype(int)
    A = fp.IntegerMatrix.from_matrix(lattice_adjusted)
    A_LLL = fp.LLL.reduction(A)
    lattice_reduced = np.array(list(A_LLL)).T/scale # transpose back
    return A_LLL, lattice_reduced

def to_L2(V, Q):
    '''
    Transform the tropical Jacobian to a standard flat torus with the L2 metric.
    Parameters:
    V: np.array, the tropical transform of a metric graph
    Q: np.array, the tropical polarization matrix
    '''
    # compute the eigenvalue decomposition of the tropical polarization matrix
    # note that np.linalg.eigh only works for symmetric matrices
    eigval, eigvec = np.linalg.eigh(Q)
    Q_sqrt = eigvec @ np.diag(np.sqrt(eigval)) @ eigvec.T
    Q_sqrt_inv = eigvec @ np.diag(1.0/np.sqrt(eigval)) @ eigvec.T

    # the transformation is given by Q^{-1/2} V
    V_transformed = Q_sqrt_inv @ V

    return V_transformed, Q_sqrt

def exact_L2_distance(V_transformed, Q_sqrt, scale=1e6):
    '''
    compute the exact distance matrix in a flat torus with the L2 metric.
    Parameters:
    V_transformed: vectors in a flat torus with the L2 metric
    Q_sqrt: the square root of the tropical polarization matrix
    scale: int, the scaling factor for the lattice basis in fpylll
    '''

    # compute the LLL-reduced lattice basis
    Q_scaled = LLL_reduced_lattice(Q_sqrt, scale)[0]

    # compute the distance matrix
    n = V_transformed.shape[1]
    D = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1, n):
            # Convert target vector to integer coordinates. float64 is not supported!
            target_scaled = np.round((V_transformed[:,i] - V_transformed[:,j]) * scale).astype("int")
            cv = fp.CVP.closest_vector(Q_scaled, target_scaled) # cv is a tuple!
            D[i,j] = np.linalg.norm(target_scaled-np.array(cv))/scale
            D[j,i] = D[i,j]
    
    return D

def babai_rounding(V_transformed, Q_sqrt_LLL):
    '''
    compute the approximate distance matrix using Babai's rounding algorithm.
    Parameters:
    V_transformed: vectors in a flat torus with the L2 metric
    Q_sqrt_LLL: the LLL-reduced square root of the tropical polarization matrix
    '''

    # compute the inverse of Q_sqrt_LLL
    Q_inv = np.linalg.inv(Q_sqrt_LLL)

    # compute the approximate distance matrix
    n = V_transformed.shape[1]
    D_approx = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1, n):
            target = V_transformed[:,i] - V_transformed[:,j]
            int_coefficients = np.round(Q_inv @ target)
            D_approx[i,j] = np.linalg.norm(target-Q_sqrt_LLL@int_coefficients)
            D_approx[j,i] = D_approx[i,j]
    
    return D_approx

def babai_nearest_plane_fpylll(V_transformed, Q_sqrt, scale=1e6):
    '''
    compute the approximate distance matrix using Babai's nearest plane algorithm in fpylll.
    '''
    
    # compute the LLL-reduced lattice basis
    Q_scaled = LLL_reduced_lattice(Q_sqrt, scale)[0]
    # compute the Gram-Schmidt orthogonalization of the lattice basis
    M = fp.GSO.Mat(Q_scaled)
    _ = M.update_gso()

    # compute the approximate distance matrix
    n = V_transformed.shape[1]
    D_approx = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1, n):
            # Convert target vector to integer coordinates
            target_scaled = np.round((V_transformed[:,i] - V_transformed[:,j]) * scale).astype("int")
            int_coefficients = M.babai(target_scaled.tolist()) # the coefficients are relative to Q_scaled!
            cv = Q_scaled.multiply_left(int_coefficients)
            D_approx[i,j] = np.linalg.norm(target_scaled-np.array(cv))/scale
            D_approx[j,i] = D_approx[i,j]
    
    return D_approx


def babai_nearest_plane_original(V_transformed, Q_sqrt_LLL):
    '''
    compute the approximate distance matrix using Babai's nearest plane algorithm.
    Parameters:
    V_transformed: vectors in a flat torus with the L2 metric
    Q_sqrt_LLL: the LLL-reduced square root of the tropical polarization matrix
    '''

    # compute the Gram-Schmidt orthogonalization of the lattice basis
    g = Q_sqrt_LLL.shape[1]
    Q_GSO = Q_sqrt_LLL.copy() 
    for i in range(g): 
        sum = np.zeros(g)
        # subtract projections onto previous orthogonal vectors
        for j in range(i):
            proj = np.dot(Q_GSO[:,i], Q_GSO[:,j]) / np.dot(Q_GSO[:,j], Q_GSO[:,j]) * Q_GSO[:,j]
            sum += proj
        Q_GSO[:,i] -= sum

    # compute the approximate distance matrix
    n = V_transformed.shape[1]
    D_approx = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1, n):
            target = V_transformed[:,i] - V_transformed[:,j]
            w = target
            cv = np.zeros(g)
            # compute the closest vector using Babai's nearest plane algorithm 
            for k in range(g-1,-1,-1):
                proj = np.dot(w, Q_GSO[:,k]) / np.dot(Q_GSO[:,k], Q_GSO[:,k])
                coeff = np.round(proj)
                cv += coeff * Q_sqrt_LLL[:,k] 
                w = w - (proj - coeff) * Q_GSO[:,k] - coeff * Q_sqrt_LLL[:,k]
            # compute the approximate distance
            D_approx[i,j] = np.linalg.norm(target-cv)
            D_approx[j,i] = D_approx[i,j]
    
    return D_approx

def babai_nearest_plane_QR(V_transformed, Q_sqrt):
    '''
    compute the approximate distance matrix using Babai's nearest plane algorithm and QR decomposition.
    '''

    # compute the QR decomposition of Q_sqrt
    Q, R = np.linalg.qr(Q_sqrt)
    # rotate the vectors
    V_rotated = Q.T @ V_transformed

    # compute the LLL-reduction of R
    R_LLL = LLL_reduced_lattice(R)[1]
    # check if R_LLL is a valid upper triangular matrix
    assert np.allclose(np.triu(R_LLL), R_LLL)

    # compute the approximate distance matrix
    g, n = V_transformed.shape
    D_approx = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1, n):
            target = V_rotated[:,i] - V_rotated[:,j]
            # compute the integral coefficients from bottom to top
            coeff = np.zeros(g) 
            sum = 0
            for k in range(g-1,-1,-1):
                for l in range(g-1,k,-1):
                    sum += coeff[l] * R_LLL[k,l]
                coeff[k] = np.round((target[k] - sum)/R_LLL[k,k])
            # compute the approximate distance
            D_approx[i,j] = np.linalg.norm(target-R_LLL@coeff)
            D_approx[j,i] = D_approx[i,j]

    return D_approx

def FS_distance(C, V, Q, solver="cbc", mipgap=1e-2):
    '''
    Compute the Foster--Zhang distance matrix using pyomo and coincbc solver.
    Parameters:
    C: np.array, the (full) cycle-edge incidence matrix
    V: np.array, the tropical transform of a metric graph
    Q: np.array, the tropical polarization matrix
    '''

    # transform vectors to the Albanese torus
    V_transformed = np.linalg.inv(Q) @ V

    # formulate the MILP problem
    g, n = V_transformed.shape
    m = C.shape[1]
    dist_mat = np.zeros((n,n))

    for i in range(n):
        for j in range(i+1,n):
            target = V_transformed[:,i] - V_transformed[:,j]
            b = C.T @ target 
            model = pyo.ConcreteModel()
            # variable x is the integral coefficients
            model.x = pyo.Var(range(g), domain=pyo.Integers)
            # variable y is the distance, set bound for faster computation
            model.y = pyo.Var(bounds=(0,100))
            # objective is to minimize the distance
            model.OBJ = pyo.Objective(expr = model.y)
            model.con = pyo.ConstraintList()
            for k in range(m):
                expr1 = -model.y
                expr2 = -model.y
                for r in range(g):
                    expr1 += C[r,k]*model.x[r]
                    expr2 += -C[r,k]*model.x[r]
                model.con.add(expr1<=b[k])
                model.con.add(expr2<=-b[k])
            solver_instance = pyo.SolverFactory(solver)
            if solver == "cbc":
                solver_instance.options['ratioGap'] = mipgap
            if solver == "ipopt":
                solver_instance.options['acceptable_tol'] = mipgap
            if solver == "glpk":
                solver_instance.options['mipgap'] = mipgap
            if solver == "scip":
                solver_instance.options['limits/gap'] = mipgap
            result = solver_instance.solve(model,tee=False)
            dist_mat[i,j] = model.y.value
            dist_mat[j,i] = dist_mat[i,j]

    return dist_mat

def truncate_mat(mat, threshold):
    '''
    Truncate the distance matrix to drop possibly wrong values.
    '''

    mat_truncated = np.where(mat > threshold, np.inf, mat)

    return mat_truncated

def draw_lattice_2D(lattice, ax):
    '''
    Draw the fundamental domain of a flat torus given by a matrix of lattice basis.
    '''
    u = lattice[:,0].reshape(1,-1)
    v = lattice[:,1].reshape(1,-1)

    l1 = np.vstack((-0.2*u, 1.2*u))
    l2 = np.vstack((-0.2*u+v, 1.2*u+v))
    l3 = np.vstack((-0.2*v, 1.2*v))
    l4 = np.vstack((-0.2*v+u, 1.2*v+u))
    ax.plot(l1[:,0], l1[:,1], 'k--',
            l2[:,0], l2[:,1], 'k--',
            l3[:,0], l3[:,1], 'k--', 
            l4[:,0], l4[:,1], 'k--')

def translate_to_fundamental_domain(points, lattice):
    '''
    Translate vectors to the fundamental domain for better visualization
    '''
    # compute the coordinates of points in the basis of the lattice
    X = np.linalg.inv(lattice) @ points
    # compute the translations
    N = np.floor(X)
    # translate the points to the fundamental domain
    points_new = points - lattice @ N
    return points_new