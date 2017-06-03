#-*- coding: utf-8 -*-

"""Laplace solver"""

import numpy as np
import matplotlib.pyplot as plt

# --- Constants --- #

steps = 500                 # number of iterations
a = 10                      # size in x
b = 15                      # size in y
delta = 1                   # fineness

x, y = np.meshgrid(np.arange(0, a), np.arange(0, b))        # grid
solution = np.empty((b, a))                                 # holds solution

# --- Boundary conditions --- #

# test boundary conditions
# phixb = x[-1:, :]
# phix0 = 0
# phiay = 10
# phi0y = 0

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
    for i in range(1, b-1, delta):
        for j in range(1, a-1, delta):
            solution[i, j] = 0.25 * (solution[i+1][j] + solution[i-1][j] + solution[i][j+1] + solution[i][j-1])

# --- Analytical --- #





# --- Plot --- #

colorinterpolation = 50
colourMap = plt.cm.jet

plt.contourf(x, y, solution, colorinterpolation, cmap=colourMap)
plt.colorbar()

plt.show()
