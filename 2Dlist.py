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
def gauss(x, a, b, c):

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

x_out = 4
x_min = -10
x_max = 10
x_step = (x_max-x_min)/x_out

y_out = 4
y_min = -12
y_max = 12
y_step = (y_max-y_min)/y_out

x = []
z2D = []

sigma_x = 3
x_0 = 0
A_x = 1/(sigma_x * sqrt(2*pi))

sigma_y = 4
y_0 = 0
A_y = 1/(sigma_y * sqrt(2*pi))

z2D = []

for i in range (0, x_out): 
    x.append(x_min+ i * x_step)
    y = []
    z = []
    for j in range (0, y_out):
        y.append(y_min+ j * y_step)
        z.append(gauss(x[i],A_x, x_0,sigma_x)*gauss(y[j],A_y,y_0,sigma_y))
    z2D.append(z)
    

print(z2D)


