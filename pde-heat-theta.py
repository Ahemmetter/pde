#-*- coding: utf-8 -*-

"""Laplace solver"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation as anim
import scipy.linalg

# --- Constants --- #

dt = 0.001                          # time step
C = 0.2                             # stability factor for explicit scheme (<= 0.5)
dx = np.sqrt(dt/C)                  # space step
length = 1                          # domain size
duration = 6                        # duration of simulation
theta = 2./3                        # 0: forward (explicit), 1: backward (implicit), 1/2: Crank-Nicolson
lx = length/dx                      # number of space steps
lt = duration/dt                    # number of time steps
t0 = 0                              # start time

print "time step: " + str(dt)
print "space step: " + str(dx)

# --- Arrays and Matrices --- #

x = np.linspace(0, length, lx)      # position array
t = np.linspace(0, duration, lt)    # time array

T = np.zeros(lx)                    # temperature array
T0 = np.sin(np.pi * x)              # initial condition
Ti = T0                             # copy to plot initial condition

A = np.zeros((len(x), len(x)))      # matrices
b = np.zeros(len(x))

for i in range(1, len(x)-1):        # assembles matrix A
    A[i, i - 1] = -C/(2*dt)         # bi
    A[i, i + 1] = -C/(2*dt)
    A[i, i] = (1+C)/dt              # ai

A[0, 0] = A[len(x)-1, len(x)-1] = 1

def init():
    """Initializes the objects line1, line2 and line3 as empty arrays, which will be filled with the data that should
    be plotted."""

    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    return line1, line2, line3,

def animate(n):
    """Main function. Calculates the values numerically, analytically and """

    # numerical solution
    b[1:-1] = (1 - theta) / dx ** 2 * Ti[0:-2] \
        + (1 / dt - 2 * (1 - theta) / dx ** 2) * Ti[1:-1] \
        + (1 - theta) / dx ** 2 * Ti[2:]
    b[0] = b[len(x)-1] = 0
    T[:] = scipy.linalg.solve(A, b)
    Ti[:] = T
    line1.set_data(x, T)

    # analytical solution
    Ta = np.sin(np.pi * x) * np.exp(-np.pi ** 2 * (t0 + n*dt))
    line2.set_data(x, Ta)

    # error
    Te = abs(Ta - T)
    line3.set_data(x, Te)

    return line1, line2, line3,

# --- Animation --- #

fig, (a0, a1) = plt.subplots(2, sharex = True, gridspec_kw = {"height_ratios": [3, 1]})
a0.set_xlim([0, length])
a0.set_ylim([0, 1.0 * 1.2])
a1.set_ylim(ymin = 0, ymax = 0.01)
a0.plot(x, T0, "r--")                                                           # initial condition
line1, = a0.plot([], [])                                                        # plots line objects
line2, = a0.plot([], [])                                                        # plots line objects
line3, = a1.plot([], [])

anim = anim.FuncAnimation(fig, animate, init_func = init,                       # animation command
                               frames = int(lt), interval = 40, blit = False, repeat = False)

#anim.save('pde-heat-theta.mp4',fps=30)                                         # optional command for saving as mp4
plt.show()