# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 13:07:51 2021

@author: AYUSH
"""
from scipy.integrate import solve_ivp
from numpy import sin,cos
import matplotlib.pyplot as pt
import numpy as np
from mpl_toolkits import mplot3d

tspan=[0,8]
teval=np.linspace(tspan[0], tspan[1], 1000)
g=9.81

def e_l_eqns(t,u):
    v=-0.25
    r=6+(v*t)
    theta,omega,phi,psi=u
    du=[omega,(sin(theta)*cos(theta)*(psi**2))-((2*v*omega)/r)+((g*sin(theta))/r),psi,-((2*v*psi)/r)-(2*cos(theta)*omega*psi/sin(theta))]
    return du

res=solve_ivp(e_l_eqns, tspan, [2.96,0.1,0.1,0.1], t_eval=teval)
theta,omega,phi,psi=res.y
t=res.t
r1=6 - (0.25*t) 
x=r1*sin(theta)*cos(phi)
y=r1*sin(theta)*sin(phi)
z=r1*cos(theta)
pt.plot(t,theta,color='red')
pt.show()
pt.plot(t,omega,color='blue')
pt.show()
pt.plot(t,phi)
pt.show()
pt.plot(t,psi)
pt.show()
pt.style.use('seaborn-poster')
g = pt.figure(figsize = (12,12))
ax = pt.axes(projection='3d')
ax.grid()
ax.plot3D(x, y, z)
ax.set_title('3D Parametric Plot')
ax.set_xlabel('x', labelpad=20)
ax.set_ylabel('y', labelpad=20)
ax.set_zlabel('z', labelpad=20)
pt.show()