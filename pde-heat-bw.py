#-*- coding: utf-8 -*-

"""Backward solver (implicit) for the 1D heat equation"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation as anim
import scipy.linalg

# --- Constants --- #

dt = 0.001                                              # time step
C = 0.4                                                 # stability factor <= 0.5 to be stable
dx = np.sqrt(dt/C)                                      # space step
length = 1                                              # domain size
duration = 6                                            # duration

print "time step: " + str(dt)
print "space step: " + str(dx)

# --- Arrays and Matrices --- #

x = np.linspace(0, length, length/dx)                   # position array
t = np.linspace(0, duration, duration/dt)               # time array

T = np.zeros(len(x))                                    # temperature array
T0 = np.sin(np.pi * x)                                  # initial condition
Ti = T0                                                 # copy to plot IC

A = np.zeros((len(x), len(x)))                          # empty matrix A, vector b
b = np.zeros(len(x))

for i in range(1, len(x)-1):                            # assemble matrix A
    A[i, i - 1] = -C
    A[i, i + 1] = -C
    A[i, i] = 1 + 2 * C

A[0, 0] = A[len(x)-1, len(x)-1] = 1

def init():                                             # creates empty line
    line.set_data([], [])
    return line

def animate(n):                                         # solver
    b[1:-1] = -Ti[1:-1]
    b[0] = b[len(x)-1] = 0
    T[:] = scipy.linalg.solve(A, b)
    Ti[:] = T
    line.set_data(x, T)
    return line,

# --- Animation --- #

fig = plt.figure()                                      # creates figure
ax = plt.axes(xlim = (0, length), ylim = (0, 1.0*1.2))  # creates axes
ax.plot(x, T0, "r--")                                   # initial condition
line, = ax.plot([], [])                                 # plots line objects

anim = anim.FuncAnimation(fig, animate, init_func = init,           # animation command
                               frames = len(t), interval = 40, blit = False, repeat = False)

#anim.save('pde-heat-bw.mp4',fps=30)                                # optional command for saving as mp4
plt.show()