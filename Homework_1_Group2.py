'''Homework 1, Computational Photonics, SS 2020:  FD mode solver.
'''
import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import eigs
import matplotlib.pyplot as plt


def guided_modes_1DTE(prm, k0, h):
    """Computes the effective permittivity of a TE polarized guided eigenmode.
    All dimensions are in µm.
    Note that modes are filtered to match the requirement that
    their effective permittivity is larger than the substrate (cladding).

    Parameters
    ----------
    prm : 1d-array
        Dielectric permittivity in the x-direction
    k0 : float
        Free space wavenumber
    h : float
        Spatial discretization

    Returns
    -------
    eff_eps : 1d-array
        Effective permittivity vector of calculated modes
    guided : 2d-array
        Field distributions of the guided eigenmodes
    """

    # 1D mode operator
    adj_elements = np.ones(len(prm)-1)/ h**2
    diag_elements = -2/h**2 + k0**2 * prm
    L = (np.diag(diag_elements, 0) + np.diag(adj_elements, -1) + np.diag(adj_elements, 1))/ k0**2

    eff_eps, field = np.linalg.eig(L)
    # f = plt.figure()
    # plt.plot(x,prm-e_substrate, label = "permittivity")
    # plt.plot(x, guided_modes[:,0,i], label = "field")
    # plt.plot(eff_eps)
    # print(e_substrate)
    # print(max(prm))
    # plt.show()
    for i in range(len(eff_eps)):
        if (eff_eps[i] <= e_substrate or eff_eps[i] >= max(prm)):
             eff_eps[i] = 0
             field[:,i] = 0
    # g = plt.figure()
    # plt.plot(eff_eps)
    # plt.show()
    guided = field[:,np.where(field[1]!=0)]
    return eff_eps[np.where(eff_eps!=0)], guided
    pass



def guided_modes_2D(prm, k0, h, numb):
    """Computes the effective permittivity of a quasi-TE polarized guided
    eigenmode. All dimensions are in µm.

    Parameters
    ----------
    prm  : 2d-array
        Dielectric permittivity in the xy-plane
    k0 : float
        Free space wavenumber
    h : float
        Spatial discretization
    numb : int
        Number of eigenmodes to be calculated

    Returns
    -------
    eff_eps : 1d-array
        Effective permittivity vector of calculated eigenmodes
    guided : 3d-array
        Field distributions of the guided eigenmodes
    """
    N = np.shape(prm)[0]

    # initializing the 2D mode operator
    L = sps.lil_matrix(sps.eye(N**2))
    for j in range(0,N): #
        for k in range(0,N):
            index = k*N + j # mein Punkt
            L[index,index] = -4/h**2 + k0**2 * prm_2D[k,j]
            if (index +1 < N**2):
                L[index,index+1] = 1
            if index -1 >= 0:
                L[index,index-1] = 1
            if index +N < N**2:
                L[index,index+N] = 1
            if index -N >= 0:
                L[index,index-N] = 1

    eigenvalues, eigenvectors = eigs(L, k = numb, which = "LR", maxiter = 1000)
    return eigenvalues, eigenvectors
    pass


# Test Values
grid_size     = 50 #100 # please choose appropriate value = N
number_points = 100 #200 # please choose appropriate value
# h             = ?? # please choose appropriate value
lam           = 0.78
k0            = 2*np.pi/lam
e_substrate   = 1.5**2
delta_e       = 1.5e-2
w             = 15.0

h = grid_size/ number_points
x = np.linspace(-grid_size/2, grid_size/2, number_points)
y = x

# Gaussian waveguide profile
prm_1D = e_substrate + delta_e * np.exp(- (x/w)**2)
prm_2D = np.array([[e_substrate + delta_e * np.exp(- (ix**2 + iy**2)/w**2) for ix in x] for iy in y])

# TASK 1
# checking the gaussian permittivity profile and therefor suitable x values
# plt.figure()
# plt.plot(x,prm_1D)
# plt.show()
# eff_eps, guided_modes = guided_modes_1DTE(prm_1D, k0, h)
# print(eff_eps)
# print(guided_modes)
# for i in range(len(eff_eps)):
#     f = plt.figure()
#     plt.plot(x, prm_1D - e_substrate, label = "permittivity")
#     plt.plot(x, guided_modes[:,0,i], label = "field")
#     plt.show()

# TASK 2
print(guided_modes_2D(prm_2D, k0, h, 10))