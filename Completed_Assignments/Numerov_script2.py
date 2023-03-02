# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 16:09:23 2023

@author: mgoryll
"""

# Third Numerov algorithm script posted on  
# https://liu-group.github.io/interactive-Numerov-PIB/
#
# Adapted to use the slider in matplotlib
#

from matplotlib.widgets import Slider
import numpy as np
import matplotlib.pyplot as plt
import time

# An improved version that handles normalization
class NumerovSolverPIB_v2:
    def __init__(self, xlower, xupper, npoints=1000, n = 1):
        self.n = n
        self.xlower = xlower
        self.xupper = xupper
        self.npoints = npoints
        self.x = np.linspace(self.xlower,self.xupper,self.npoints)
        self.delta = self.x[1]-self.x[0]
        self.x1 = 1e-2 # A small positive value as the initial guess the psi value at the second grid point
        self.psi_left = None
        self.psi_right = None
        self.k2 = n**2 * np.pi**2
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
        # Calculate the integral of probability distribution
        self.prob_left = np.trapz(np.power(self.psi_left,2),self.x)
        # Normalize the function
        self.psi_left = self.psi_left/np.sqrt(self.prob_left)

    def Numerov_right(self):
        if self.n % 2 == 0:
          Constant = -1
        else:
          Constant = 1
        self.psi_right = np.zeros(len(self.x))
        self.psi_right[-2] = self.x1
        for i in range(len(self.x)-2,0,-1):
            self.psi_right[i-1] = self.propagateNumerov(self.x[i+1],self.x[i],self.psi_right[i+1],self.psi_right[i])
        # Calculate the integral of probability distribution
        self.prob_right = np.trapz(np.power(self.psi_right,2),self.x)
        # Normalize the function
        self.psi_right = self.psi_right/np.sqrt(self.prob_right)*Constant

def numerov2(numPoints, n = 1, xlower = 0, xupper = 1):
  # Create the figure and the line that we will manipulate
  fig, ax = plt.subplots()
  solver=NumerovSolverPIB_v2(xlower,xupper,numPoints, n = n)
  start=time.time()
  solver.Numerov_left()
  solver.Numerov_right()
  end=time.time()
  timetaken=end-start

  prob_left = np.trapz(np.power(solver.psi_left,2),solver.x)
  prob_right = np.trapz(np.power(solver.psi_right,2),solver.x)

  fig.suptitle("Numerov solution to PIB, "
            + r"$\int\psi_\mathrm{left}$=" + "{:.2g}".format(prob_left)
            + r", $\int\psi_\mathrm{right}=$" + "{:.2g}".format(prob_right)
            + "\nTime taken: {:.2f} ms".format(timetaken*1000))

  # plot both the lines and points to make it visually better
  line_l, = ax.plot(solver.x, solver.psi_left, c='b',label=r'$\psi_\mathrm{left}$')
  scatter_l = ax.scatter(solver.x, solver.psi_left, c='b')

  line_r, = ax.plot(solver.x, solver.psi_right, c='r',label=r'$\psi_\mathrm{right}$')
  scatter_r = ax.scatter(solver.x, solver.psi_right, c='r')

  # Plot the reference analytical solution
  x = np.linspace(0,1,1000)
  L = 1
  psi = (np.sqrt(2/L))*np.sin(n*np.pi*x/L)
  ax.plot(x, psi, c='k',label=r'$\psi$')

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
      solver=NumerovSolverPIB_v2(xlower, xupper,n_slider.val, n = n)
      start=time.time()
      solver.Numerov_left()
      solver.Numerov_right()
      end=time.time()
      timetaken=end-start

      prob_left = np.trapz(np.power(solver.psi_left,2),solver.x)
      prob_right = np.trapz(np.power(solver.psi_right,2),solver.x)
      fig.suptitle("Numerov solution to PIB, "
                + r"$\int\psi_\mathrm{left}$=" + "{:.2g}".format(prob_left)
                + r", $\int\psi_\mathrm{right}=$" + "{:.2g}".format(prob_right)
                + "\nTime taken: {:.2f} ms".format(timetaken*1000))

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