#-*- coding: utf-8 -*-

"""Forward solver (explicit) for the 1D heat equation"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation as anim

# --- Constants --- #

dt = 0.001                                                      # time step
C = 0.5                                                         # stability factor (<= 0.5)
dx = np.sqrt(dt/C)                                              # space step
length = 1                                                      # domain size
duration = 6                                                    # duration

print "time step: " + str(dt)
print "space step: " + str(dx)

# --- Arrays and Matrices --- #

x = np.linspace(0, length, length/dx)                           # position array
t = np.linspace(0, duration, duration/dt)                       # time array

T = np.zeros(len(x))                                            # temperature array
T0 = np.sin(np.pi * x)                                          # initial condition
Ti = T0                                                         # copy to plot IC

def init():
    line.set_data([], [])
    return line

def animate(n):
    T[1:-1] = Ti[1:-1] + C * (Ti[:-2] - 2 * Ti[1:-1] + Ti[2:])
    T[0] = 0
    T[-1] = 0
    Ti[:] = T
    line.set_data(x, T)
    return line,

# --- Animation --- #

fig = plt.figure()                                              # creates figure
ax = plt.axes(xlim = (0, length), ylim = (0, 1.0*1.2))          # creates axes
ax.plot(x, T0, "r--")                                           # initial condition
line, = ax.plot([], [])                                         # plots line objects

anim = anim.FuncAnimation(fig, animate, init_func = init,       # animation command
                               frames = len(t), interval = 40, blit = False, repeat = False)

#anim.save('pde-heat-fw.mp4',fps=30)                            # optional command for saving as mp4
plt.show()


