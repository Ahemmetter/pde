#-*- coding: utf-8 -*-

"""Laplace solver using the inverse matrix method."""

import numpy as np
import scipy as sp
import scipy.linalg
import matplotlib.pyplot as plt

# --- Constants --- #

n = 10                  # size of the entire field (below 100)
size = n - 2            # size of matrix that needs to be calculated (m)
sol = np.zeros((n, n))  # solution grid
x, y = np.meshgrid(np.arange(0, n), np.arange(0, n))        # grid

# --- Boundary conditions --- #

phixb = np.reshape(np.sin(x[-1:, :])/np.sin(n), (n,))
phix0 = np.zeros((n,))
phiay = np.reshape(np.sinh(y[:, -1:])/np.sinh(n), (n,))
phi0y = np.zeros((n,))

sol[-1, :] = phixb
sol[0, :] = phix0
sol[:, -1] = phiay
sol[:, 0] = phi0y

# --- Functions --- #

def makeA(m):
    B_diag = -4 * np.eye(m)
    B_u = np.diag([1] * (m - 1), 1)
    B_l = np.diag([1] * (m - 1), -1)
    B = B_diag + B_u + B_l
    A = scipy.linalg.block_diag(*([B] * m))
    D_u = np.diag(np.ones(m * (m - 1)), m)
    D_l = np.diag(np.ones(m * (m - 1)), -m)
    A += D_u + D_l
    return A

def makeb(m):
    b = np.zeros(m ** 2)
    b[:m] = np.reshape(phix0, (n,) )[1:-1]                      # phix0
    b[(m-1)*m:m*m] = -np.reshape(phixb, (n,) )[1:-1]            # phixb
    b[(m - 1): m * m + 1 : m] = -np.reshape(phi0y, (n,) )[1:-1] # phi0y
    b[: m * m: m] = np.reshape(phiay, (n,) )[1:-1]              # phiay
    return b

def solve(makeA, makeb, m, s):
    A = makeA(m)
    b = makeb(m)
    U = sp.linalg.solve(A, b)
    s[1:-1, 1:-1] = U.reshape((m,m))
    return s

final = solve(makeA, makeb, size, sol)

# --- Plot --- #

colorinterpolation = 100
colourMap = plt.cm.jet
plt.contourf(x, y, final, colorinterpolation, cmap=colourMap)
plt.colorbar()

plt.show()