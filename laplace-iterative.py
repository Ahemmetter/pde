#-*- coding: utf-8 -*-

"""Laplace solver"""

import numpy as np
import matplotlib.pyplot as plt

# --- Constants --- #

steps = 300                 # number of iterations
a = 1.0                      # size in x
b = 1.5                      # size in y
delta = 0.01                   # fineness

x, y = np.meshgrid(np.arange(0, a, delta), np.arange(0, b, delta))        # grid
solution = np.empty((int(b/delta), int(a/delta)))                                 # holds solution
print solution.shape

# --- Boundary conditions --- #

# test boundary conditions
# phixb = x[-1:, :]
# phix0 = 0
# phiay = 10
# phi0y = 0

# sina = np.sin(a)
phixb = np.sin(x[-1:, :])/np.sin(a)
phix0 = 0
phiay = np.sinh(y[:, -1:])/np.sinh(b)
phi0y = 0

solution[-1:, :] = phixb
solution[:1, :] = phix0
solution[:, -1:] = phiay
solution[:, :1]  = phi0y

# --- Calculation --- #

for iteration in range(0, steps):
    for i in range(1, int(b/delta)-1):
        for j in range(1, int(a/delta)-1):
            solution[i, j] = 0.25 * (solution[i+1][j] + solution[i-1][j] + solution[i][j+1] + solution[i][j-1])

# --- Plot --- #

colorinterpolation = 50
colourMap = plt.cm.jet

plt.contourf(x, y, solution, colorinterpolation, cmap=colourMap)
plt.colorbar()

plt.show()
