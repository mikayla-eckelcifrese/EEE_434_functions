# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 11:13:45 2023

@author: mgoryll
"""

# First Numerov algorithm script posted on  
# https://liu-group.github.io/interactive-Numerov-PIB/
#
# Adapted to use the slider in matplotlib
#
# Spyder needs to be configured from "inline" to "automatic" matplotlib plotting
#

# %matplotlib widget
# import ipywidgets as widgets
from matplotlib.widgets import Slider
import numpy as np
import matplotlib.pyplot as plt
import time

class NumerovSolverPIB:
    def __init__(self, xlower, xupper, npoints=1000):
        self.xlower = xlower
        self.xupper = xupper
        self.npoints = npoints
        self.x = np.linspace(self.xlower,self.xupper,self.npoints)
        self.delta = self.x[1]-self.x[0]
        self.x1 = 1e-2 # A small positive value as the initial guess the psi value at the second grid point
        self.psi_left = None
        self.psi_right = None
        self.k2 = np.pi**2
        self.prob_left = None
        self.prob_right = None

    def propagateNumerov(self,x0,x1,psi0,psi1):
        psi2 = 2*psi1 - psi0 - self.k2*psi1*self.delta**2
        return psi2

    def Numerov_left(self):
        self.psi_left = np.zeros(len(self.x))
        self.psi_left[1] = self.x1
        for i in range(1,len(self.x)-1):
            self.psi_left[i+1] = self.propagateNumerov(self.x[i-1],self.x[i],self.psi_left[i-1],self.psi_left[i])

    def Numerov_right(self):
        self.psi_right = np.zeros(len(self.x))
        self.psi_right[-2] = self.x1
        for i in range(len(self.x)-2,0,-1):
            self.psi_right[i-1] = self.propagateNumerov(self.x[i+1],self.x[i],self.psi_right[i+1],self.psi_right[i])

# Number of discrete points to use in the range
# n_options = widgets.IntSlider(
#     value=10,
#     min=3,
#     max=101,
#     step=5,
#     description=r'npoints:',
#     disabled=False,
# )
def numerov1(numPoints):
  # Create the figure and the line that we will manipulate
  fig, ax = plt.subplots()
  solver=NumerovSolverPIB(0,1,numPoints)
  start=time.time()
  solver.Numerov_left()
  solver.Numerov_right()
  end=time.time()
  timetaken=end-start

  fig.suptitle("Numerov solution to PIB\n" + "Time taken: {:.2f} ms".format(timetaken*1000))

  # plot both the lines and points to make it visually better
  line_l, = ax.plot(solver.x, solver.psi_left, c='b',label=r'$\psi_\mathrm{left}$')
  scatter_l = ax.scatter(solver.x, solver.psi_left, c='b')

  line_r, = ax.plot(solver.x, solver.psi_right, c='r',label=r'$\psi_\mathrm{right}$')
  scatter_r = ax.scatter(solver.x, solver.psi_right, c='r')
  ax.set_xlabel(r'$x$ (bohr)')
  ax.set_ylabel(r'$\psi(x)$')
  plt.legend()

  # adjust the main plot to make room for the sliders
  fig.subplots_adjust(bottom=0.25)

  # matplotlib slider
  axn= fig.add_axes([0.25, 0.1, 0.65, 0.03])
  n_slider = Slider(
      ax=axn,
      label='npoints:',
      valmin=3,
      valmax=101,
      valstep=1,
      valinit=10,
  )

  # The function to be called anytime a slider's value changes
  def update(val):
     # plt.clf()
      solver=NumerovSolverPIB(0,1,n_slider.val)
      start=time.time()
      solver.Numerov_left()
      solver.Numerov_right()
      end=time.time()
      timetaken=end-start

      # The wave function is not normalized.
      # prob = np.trapz(np.power(solver.psi_left,2),solver.x)
      fig.suptitle("Numerov solution to PIB\n" + "Time taken: {:.2f} ms".format(timetaken*1000))

      # plot both the lines and points to make it visually better
      line_l.set_xdata(solver.x)
      line_l.set_ydata(solver.psi_left)
      scatter_l.set_offsets(np.c_[solver.x, solver.psi_left])

      line_r.set_xdata(solver.x)
      line_r.set_ydata(solver.psi_right)
      scatter_r.set_offsets(np.c_[solver.x, solver.psi_right])

      fig.canvas.draw_idle()

  # register the update function with each slider
  n_slider.on_changed(update)

  plt.show()
