import numpy as np 
import matplotlib.pyplot as plt 
import scipy.special as sp 

def get_angular(l, m):
    PHI, THETA = np.mgrid[0:2*np.pi:200j, 0:np.pi:100j] #arrays of angular variables
    R = np.abs(sp.sph_harm(m, l, PHI, THETA))

    X = R * np.sin(THETA) * np.cos(PHI)
    Y = R * np.sin(THETA) * np.sin(PHI)
    Z = R * np.cos(THETA)

    print(R.shape, PHI.shape)
    r = []

    for i in range(200):
        for j in range(100):
            if PHI[i, j] == 0 or PHI[i, j] == PHI.max():
                print(R[i, j])
                r.append(R[i, j])

    r = np.array(r)
    return r



if __name__ == '__main__':
    m=0
    l=0
    r = get_angular(l, m)
    print(r)
    plt.plot(r)
    plt.show()

