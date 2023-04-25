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
# Read in initial conditions from a file
#------------------------------------------
def read2C_ary(Location, Column1, Column2):
    x=[]
    y=[]
    crs=open(Location,"r")
    Np = 0
    for columns in (raw.strip().split() 
        for raw in crs):
            x.append(columns[Column1])
            y.append(float(columns[Column2]))
            Np = Np + 1
    return x, y, Np


#--------------------------------------------------
# Read in initial conditions from a file
#--------------------------------------------------
def read3C_ary(Location, Column1, Column2, Column3):
    x=[]
    y=[]
    z=[]
    crs=open(Location,"r")
    Np = 0
    for columns in (raw.strip().split() 
        for raw in crs):
            x.append(float(columns[Column1]))
            y.append(float(columns[Column2]))
            z.append(float(columns[Column3]))
            Np = Np + 1
    return x, y, z, Np


#--------------------------------------------------
# Read in initial conditions from a file
#--------------------------------------------------
def read4C_ary(Location, Column1, Column2, Column3, Column4):
    w=[]
    x=[]
    y=[]
    z=[]
    crs=open(Location,"r")
    Np = 0
    for columns in (raw.strip().split() 
        for raw in crs):
            w.append(float(columns[Column1]))
            x.append(float(columns[Column2]))
            y.append(float(columns[Column3]))
            z.append(float(columns[Column4]))
            Np = Np + 1
    return w, x, y, z, Np


#--------------------------------------------------
# Read in initial conditions from a file
#--------------------------------------------------
def read5C_ary(Location, Column1, Column2, Column3, Column4, Column5):
    v=[]
    w=[]
    x=[]
    y=[]
    z=[]
    crs=open(Location,"r")
    Np = 0
    for columns in (raw.strip().split() 
        for raw in crs):
            v.append(float(columns[Column1]))
            w.append(float(columns[Column2]))
            x.append(float(columns[Column3]))
            y.append(float(columns[Column4]))
            z.append(float(columns[Column5]))
            Np = Np + 1
    return v, w, x, y, z, Np


#------------------------------------------
#   Frequency Modulation (Chirping) 1D
#------------------------------------------
def FM1D(t,a,S,erft):
    f0 = 1/(1+((a**2)/2))
    FM = f0*(1+(sqrt(pi)*erft*((S*(a**2)))/(4*(c*t))))
    return FM


#-------------------------------------------
# Efficient evaluation of an error function
#-------------------------------------------
def erf(x):
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911
    sgn = 1
    if x < 0:
        sgn = -1
    x = abs(x)
    t = 1.0/(1.0+p*x)
    y = 1.0 - (((((a5*t + a4)*t)+a3)*t + a2)*t +a1) *t*exp(-x*x)
    return sgn*y;


#------------------------------------------
#   Frequency Modulation (Chirping) 1D
#------------------------------------------
def projectile(v_0,theta, N):

    x = []
    y = []
    vx = []
    vy = []

    t_f = (2*v_0*sin(theta))/g

    print('t_f = ', t_f)
    t = 0

    for i in range (0,N+1):
        t = t_f * ((i)/N)
        x.append(v_0 * cos(theta) * t)
        vx.append(v_0 * cos(theta))
        y.append((v_0 * sin(theta) * t )- (0.5 * g * t**2))
        vy.append(v_0 * sin(theta) * t - g * t)

    return x,y,vx,vy



#=================
# MAIN PROGRAM
#=================

Nout = 100
launch_speed = 15.3
launch_deg = 72
launch_angle = launch_deg*(2*pi/360)

v0_x = launch_speed*cos(launch_angle)
v0_y = launch_speed*sin(launch_angle)

x, y, vx, vy = projectile(launch_speed, launch_angle, Nout)

UP = max(max(x),max(y))*1.1
DOWN = max(y)*0.1
LEFT = 1
RIGHT = max(max(x),max(y))*1.1

soa = np.array([[0, 0, RIGHT, 0], [0, 0, 0, UP], [0, 0, -LEFT, 0],[0,0,0,-DOWN]])
X, Y, U, V = zip(*soa)
soa = np.array([[0, 0, (RIGHT*0.1*v0_x/launch_speed), (RIGHT*0.1*v0_y/launch_speed)], [0, 0, 0, 0], [0, 0, 0, 0],[0,0,0,0]])
A, B, C, D = zip(*soa)

fig, ax = plt.subplots(figsize =(7,7))

ax.tick_params(axis = 'x', which = 'both', bottom=False, top = False, labelbottom = False)
ax.tick_params(axis = 'y', which = 'both', left=False, right = False, labelleft = False)

v0color = '#f7022a'
trajcolor = '#56ae57'

#ax.set_xticks(np.arange(-LEFT,RIGHT+1, step = 1))
#ax.set_yticks(np.arange(-DOWN,UP+1, step = 1))
ax.set_xlabel(r'Horizontal Position [m]')
ax.set_ylabel(r'Vertical Position [m]')
#ax.plot(x_axis,x_axisy, linewidth=2, color = 'k')
#ax.plot(y_axisx,y_axis, linewidth=2, color = 'k')
ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1, color=['k','k','k','k'])
ax.quiver(A, B, C, D,  angles='xy', scale_units='xy', scale=1, color=[v0color,'g','b','c'])
ax.plot(0,0, color = v0color,label = '$v_0$=%.*f m/s'%(2, launch_speed))
ax.plot(x,y, ':', color = trajcolor, label = 'Trajectory')
ax.plot(0,0, 'k',label = 'Launch Angle %i $^0$ '%launch_deg)
#ax.plot(0,0, 'c',label = 'D')
ax.set_xlim([-LEFT, RIGHT])
ax.set_ylim([-DOWN, UP])
ax.grid(linewidth = 0.5, color ='#4e5481')
ax.legend(borderpad=0.5, fontsize = 12, loc =9)
plt.show()


plt.savefig('vector2D.pdf', format='pdf', dpi=2000)

plt.show()