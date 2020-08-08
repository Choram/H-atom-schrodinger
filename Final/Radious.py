import numpy as np
from scipy import constants as const
from scipy import sparse as sparse
from scipy.sparse.linalg import eigs
from matplotlib import pyplot as plt
 
hbar = const.hbar
e = const.e
m_e = const.m_e
pi = const.pi
epsilon_0 = const.epsilon_0
joul_to_eV = e
 
def calculate_potential_term(r):
    potential = e**2 / (4.0 * pi * epsilon_0) / r
    potential_term = sparse.diags((potential))
    return potential_term
 
def calculate_angular_term(r, l):
    angular = l * (l + 1) / r**2
    angular_term = sparse.diags((angular))
    return angular_term
 
def calculate_laplace_three_point(r, N):
    h = r[1] - r[0]
     
    main_diag = -2.0 / h**2 * np.ones(N)     
    off_diag  =  1.0 / h**2 * np.ones(N - 1)
    laplace_term = sparse.diags([main_diag, off_diag, off_diag], (0, -1, 1))
    return laplace_term
     
def build_hamiltonian(r, l, N):
    laplace_term =   calculate_laplace_three_point(r, N)
    angular_term =   calculate_angular_term(r, l)
    potential_term = calculate_potential_term(r)
     
    hamiltonian = -hbar**2 / (2.0 * m_e) * (laplace_term - angular_term) - potential_term
 
    return hamiltonian

def get_eigs(n, l, N=2000):
    r = np.linspace(2e-9, 0.0, N, endpoint=False)
    hamiltonian = build_hamiltonian(r, l, N)

    number_of_eigenvalues = 30
    eigenvalues, eigenvectors = eigs(hamiltonian, k=number_of_eigenvalues, which='SM')

    eigenvectors = np.array([x for _, x in sorted(zip(eigenvalues, eigenvectors.T), key=lambda pair: pair[0])])
    eigenvalues = np.sort(eigenvalues)

    return eigenvalues[n - 1 - l], eigenvectors[n - 1 - l]

if __name__ == '__main__':
    w, v = get_eigs(2, 1, 2000)
    r = np.linspace(2e-9, 0.0, 2000, endpoint=False)
    print(w.real / e)
    plt.plot(np.flip(np.abs(v/r)**2))
    plt.show()