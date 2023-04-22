# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 10:23:15 2023

@author: mgoryll
Based on the code from https://dpotoyan.github.io/Chem324/H-atom-wavef.html
"""
import matplotlib.pyplot as plt

import numpy as np
from numpy import vectorize

# Import special functions 
import scipy.special as spe


def psi_R(r, n=1, l=0):
    coeff = np.sqrt((2.0 / n) ** 3 * spe.factorial(n - l - 1) / (2.0 * n * spe.factorial(n + l)))

    laguerre = spe.assoc_laguerre(2.0 * r / n, n - l - 1, 2 * l + 1)

    return coeff * np.exp(-r / n) * (2.0 * r / n) ** l * laguerre

def psi_ang(phi, theta, l=0, m=0):
    sphHarm = spe.sph_harm(m, l, phi, theta)
    if m < 0 and (abs(m) % 2) == 0:
        return np.sqrt(2) * sphHarm.imag
    elif m < 0:
        return np.sqrt(2) * (-1) * sphHarm.imag
    elif m > 0 and (m % 2) == 0:
        return np.sqrt(2) * sphHarm.real
    elif m > 0:
        return np.sqrt(2) * (-1) * sphHarm.real

    return sphHarm.real

def HFunc(r, theta, phi, n, l, m):
    '''
    Hydrogen wavefunction // a_0 = 1

    INPUT
        r: Radial coordinate
        theta: Polar coordinate
        phi: Azimuthal coordinate
        n: Principle quantum number
        l: Angular momentum quantum number
        m: Magnetic quantum number

    OUTPUT
        Value of wavefunction
    '''

    return psi_R(r, n, l) * psi_ang(phi, theta, l, m)


r, phi, theta = np.mgrid[0:5:36j, 0:2 * np.pi:36j, 0:np.pi:36j]

@vectorize
def plot_orbital(n, l, m):
    print(f'l = {l}, m = {m}')
    psi_nlm = HFunc(r, theta, phi, n=n, l=l, m=m)

    x = r * np.sin(theta) * np.cos(phi) * psi_nlm ** 2
    y = r * np.sin(theta) * np.sin(phi) * psi_nlm ** 2
    z = r * np.cos(theta) * psi_nlm ** 2

    '''Set up the 3D Canvas'''

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    '''Make 3D plot of the probability point cloud'''

    ax.scatter(x, y, z, c=r ** 2 * psi_nlm ** 2, cmap='viridis', marker='.')

    ''' Set axes limit to keep aspect ratio 1:1:1 '''
    ax.set_xlim3d(-0.015, 0.015)
    ax.set_ylim3d(-0.015, 0.015)
    ax.set_zlim3d(-0.015, 0.015)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()

@vectorize
def plot_hybrid_orbital(n = 1, l = 0, m = 0):
    print(f'n = {n}, l = {l}, m = {m}')

    psi_nlm_s = HFunc(r, theta, phi, n=n, l=l, m=m)
    psi_nlm_p = HFunc(r, theta, phi, n=n+1, l=l+1, m=m+1)
    psi_nlm = psi_nlm_s + psi_nlm_p

    x = r * np.sin(theta) * np.cos(phi) * psi_nlm ** 2
    y = r * np.sin(theta) * np.sin(phi) * psi_nlm ** 2
    z = r * np.cos(theta) * psi_nlm ** 2

    '''Set up the 3D Canvas'''

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    '''Make 3D plot of the probability point cloud'''

    ax.scatter(x, y, z, c=r ** 2 * psi_nlm ** 2, cmap='viridis', marker='.')

    ''' Set axes limit to keep aspect ratio 1:1:1 '''
    ax.set_xlim3d(-0.015, 0.015)
    ax.set_ylim3d(-0.015, 0.015)
    ax.set_zlim3d(-0.015, 0.015)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()