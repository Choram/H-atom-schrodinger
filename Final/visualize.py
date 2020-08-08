import matplotlib.pyplot as plt 
import numpy as np 
from scipy import constants as const
import Angular
import Radious

n = 3
l = 2
m = 1

eig, r = Radious.get_eigs(n, l)
Y = Angular.get_angular(l, m)
domain = np.linspace(2e-9, 0.0, 2000, endpoint=False)
plt.plot(np.flip(np.abs(r/domain)**2), color='darkblue')
plt.show()
plt.plot(Y)
plt.show()


t = np.linspace(0, 2*np.pi, len(Y))
t1 = np.linspace(0, 3, len(r))
A = []

for i in range(len(t)):
    A.append(Y[i] * np.flip(np.abs(r/domain)**2))

print("Done!")
A = np.array(A)
A = A.T
plt.imshow(A)
plt.show()

R, P = np.meshgrid(t, t1)
X, Y = P*np.cos(R), P*np.sin(R)

plt.figure(figsize=(6, 6))
plt.pcolormesh(X, Y, A)
plt.show()
