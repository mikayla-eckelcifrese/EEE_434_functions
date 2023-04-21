# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 12:45:37 2023

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


def plot_rad_dist(n, l, ax):
    r = np.linspace(0, 100, 1000)

    R = psi_R(r, n=n, l=l)

    label = f"n = {n}, l = {l}"
    ax.plot(r, r ** 2 * R ** 2, lw=3, label = label)
    ax.set_xlabel('$r [a_0]$', fontsize=12)
    ax.set_ylabel('$r^2 R^2_{nl}(r)$', fontsize=12)
    ax.grid(True)
    ax.legend()


