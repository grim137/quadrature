import matplotlib
import random
import numpy as np
import multiprocessing as mp
import time
import sys
from scipy.integrate import dblquad
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from math import pi,sqrt,cos,sin,atan,e,log,exp


import os
from scipy.interpolate import interp1d  #Beth 6/23/2020
import scipy.integrate as integrate


#----------------------------------------------------------------
# Define physical constants
#----------------------------------------------------------------
hbar = 6.582e-16          # Planck's constant (eV*s)
e0 = 8.85e-12             # epsilon 0 (C^2/Nm^2)
c = 2.998e8               # speed of light (m/s)
me = 9.109e-31            # electron mass (kg)
r = 2.818e-15             # electron radius (m)
ec = 511e3                # electron rest energy (eV)
q = 1.602e-19             # electron charge (C)
eVtoJ = 1.602176565e-19   # eV to Joules conversion
mc2 = 511e3                # electron mass [eV]
eV2J = 1.602176565e-19   # eV to Joules conversion
hz2eV = 241799050402293.0  # Convert Hz to eV
g = 9.81                  # acceleration due to gravity (m/s^2)


#------------------------------------------
#   Gaussian Function
#------------------------------------------
def func1(x, a, b, c):

    y = a * exp((-((x-b)**2 ))/( 2*c**2 )) 

    return y


#------------------------------------------
#   Gaussian Analytic Integral
#------------------------------------------
def Int1(a, b, c):

    I = a * c * sqrt(2*pi)
    return I

#------------------------------------------
#   Quadratic
#------------------------------------------
def func2(x, a, b, c):

    y = a * x**2 + b * x +c    

    return y



#=================
# MAIN PROGRAM
#=================

# Base Function

Nout = 100
x_min = -10
x_max = 10
step = (x_max-x_min)/Nout

y = []
x = []

A = 1
B = 1
C = 1

for i in range(0, Nout):
    x.append(x_min+ i * step)
    y.append(func2(x[i],A,B,C))


# Quadrature: Riemann Sum

integral = 0

for i in range(1, Nout):
    integral = integral + y[i] * step

#print('Analytic Integral = ', Int1(A,B,C))
#print('Quadrature Integral= ', integral)


# Generate Plot

fig, ax = plt.subplots(figsize =(7,7))

#ax.set_xlabel(r'Horizontal Position [m]')
#ax.set_ylabel(r'Vertical Position [m]')

ax.plot(x,y, 'r', label = 'Function')

ax.set_xlim(x_min, x_max)
#ax.set_ylim(x_min, UP])

ax.grid(linewidth = 0.5, color ='#4e5481')
#ax.legend(borderpad=0.5, fontsize = 12, loc =9)
plt.show()


plt.savefig('quadrature.pdf', format='pdf', dpi=2000)

plt.show()