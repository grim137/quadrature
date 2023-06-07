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
import math
import pylab


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
#   Gaussian Integral
#------------------------------------------
def Int1(x_min, x_max, a, b, c):

    I = a * c * sqrt(2*pi)
    return I


#------------------------------------------
#   Quadratic
#------------------------------------------
def func2(x_0, a, b, c):

    y_0 = (a * x_0**2) + (b * x_0) + c 

    return y_0

#------------------------------------------
#   Quadractic Integral
#------------------------------------------
def Int2(x_min, x_max, a, b, c):

    I = (((a*x_max**3)/3) + ((b*x_max**2)/2) +c*x_max) - (((a*x_min**3)/3) + ((b*x_min**2)/2) +c*x_min)
    return I




#=================
# MAIN PROGRAM
#=================

# Base Function

Nout = 2000
x_min = -10
x_max = 10
step = (x_max-x_min)/Nout

y = []
x = []
y_0 = []
y_1 = []

a = 1
b = 1
c = 1
for i in range(0, Nout):
    x.append(x_min + i * step)
    y.append(func1(x[i], a, b, c))


# Quadrature: 

integral = 0

#for i in range(1, Nout):                           #Riemann Sums
    #integral = integral + y[i] * step


#for i in range (1, Nout):                           #Trapezoid Rule
    #integral = integral + ((y[i-1] + y[i])/2) * step

for i in range (0, Nout):                           #Simpson's Rule
    for s in range (2, Nout, 2):
        x_even = y[s] * 2
    for t in range (1, Nout, 2):
        x_odd = y[t] * 4
    
    integral = integral + step * (y[0] + y[i] + x_odd + x_even)/3
   #integral = integral + (step/3)((i[0] + i[Nout])(4(x_odds) + 2(x_even)))
    

print('Analytic Integral = ', Int1(x_min, x_max, a,b,c))
print("Simpson's Rule = ", integral)


# Generate Plot

fig, ax = plt.subplots(figsize =(7,7))

#ax.set_xlabel(r'Horizontal Position [m]')
#ax.set_ylabel(r'Vertical Position [m]')

ax.plot(x,y, 'r', label = 'Trajectory')

ax.set_xlim(x_min, x_max)
#ax.set_ylim(x_min, UP])

ax.grid(linewidth = 0.5, color ='#4e5481')
#ax.legend(borderpad=0.5, fontsize = 12, loc =9)
plt.show()


plt.savefig('quadrature.pdf', format='pdf', dpi=2000)

plt.show()