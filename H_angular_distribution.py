# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 11:42:14 2023

@author: mgoryll
Based on the code from https://dpotoyan.github.io/Chem324/H-atom-wavef.html 
"""

import matplotlib.pyplot as plt

from matplotlib import cm

import numpy as np
from numpy import vectorize

# import ipyvolume as ipv

# Import special functions 
import scipy.special as spe

def psi_ang(phi, theta, l=0, m=0):
    sphHarm = spe.sph_harm(m, l, phi, theta)  # sph_harm switches theta and phi

    if m < 0 and (abs(m) % 2) == 0:
        return np.sqrt(2) * sphHarm.imag
    elif m < 0:
        return np.sqrt(2) * (-1) * sphHarm.imag
    elif m > 0 and (m % 2) == 0:
        return np.sqrt(2) * sphHarm.real
    elif m > 0:
        return np.sqrt(2) * (-1) * sphHarm.real

    return sphHarm.real


phi, theta = np.linspace(0, 2 * np.pi, 100), np.linspace(0, np.pi, 100)

phi, theta = np.meshgrid(phi, theta)

@vectorize
def plot_angular_dist(l, m):
    Ylm = psi_ang(phi, theta, l=l, m=m)

    x = np.sin(theta) * np.cos(phi) * abs(Ylm)
    y = np.sin(theta) * np.sin(phi) * abs(Ylm)
    z = np.cos(theta) * abs(Ylm)

    '''Set up the 3D Canvas'''

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    ''' Normalize color bar to [0,1] scale'''

    fcolors = (Ylm - Ylm.min()) / (Ylm.max() - Ylm.min())

    '''Make 3D plot of real part of spherical harmonic'''

    ax.plot_surface(x, y, z, facecolors=cm.seismic(fcolors), alpha=0.3)

    ''' Project 3D plot onto 2D planes'''

    cset = ax.contour(x, y, z, 20, zdir='z', offset=-1, cmap='summer')
    cset = ax.contour(x, y, z, 20, zdir='y', offset=1, cmap='winter')
    cset = ax.contour(x, y, z, 20, zdir='x', offset=-1, cmap='autumn')

    ''' Set axes limit to keep aspect ratio 1:1:1 '''

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()
