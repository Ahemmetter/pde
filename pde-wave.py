#-*- coding: utf-8 -*-

"""Explicit solver for the 1D wave equation"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation as anim

# --- Constants --- #

dt = 0.01                                   # time step
dx = 0.01                                   # space step
c = 1                                       # speed of propagation
a = c * dt/dx                               # step factor
length = 1                                  # domain size
duration = 6                                # duration

# --- Arrays and Matrices --- #

x = np.linspace(0, length, length/dx)       # position array
t = np.linspace(0, duration, duration/dt)   # time array

u = np.zeros(len(x))                        # solution array
u1 = np.sin(np.pi * x)                      # IC
u2 = np.zeros(len(x))                       # IC
v = np.sin(np.pi * x)                       # IC

# initial setup

u[1:-1] = u1[1:-1] + dt * v[1:-1] + 0.5 * a ** 2 * (u1[0:-2] - 2 * u1[1:-1] + u1[2:])
u[0] = u[len(x)-1] = 0
u2[:], u1[:] = u1, u

def init():
    line.set_data([], [])
    return line

def animate(n):
    u[1:-1] = - u2[1:-1] + 2 * u1[1:-1] + a**2 * (u1[0:-2] - 2 * u1[1:-1] + u1[2:])

    # Insert boundary conditions
    u[0] = u[len(x)-1] = 0

    # Switch variables before next step
    u2[:], u1[:] = u1, u

    line.set_data(x, u)
    return line,

# --- Animation --- #

fig = plt.figure()                                          # creates figure
ax = plt.axes(xlim = (0, length), ylim = (-1.2, 1.2))       # creates axes
line, = ax.plot([], [])                                     # plots line objects

anim = anim.FuncAnimation(fig, animate, init_func = init,   # animation command
                               frames = len(t), interval = 40, blit = False, repeat = False)

#anim.save('pde-wave.mp4',fps=30)                           # optional command for saving as mp4
plt.show()


