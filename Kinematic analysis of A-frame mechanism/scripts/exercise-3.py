# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 16:45:12 2021

@author: AYUSH
"""
import math
import matplotlib.pyplot as pt
import numpy as np


# fixed parameters
AB=5
BC=6
CD=5
DA=6
theta1=math.pi

#time limit
tspan=[0,45]
teval=np.linspace(tspan[0], tspan[1], 4501)

#independent coordinate
zetain=math.pi/6 # angle subtended by the A-frame [rad]
beta=90*math.pi/(180*45)
theta2=zetain+beta*teval

x=-AB*np.cos(theta2)-DA*np.cos(theta1)
y=-AB*np.sin(theta2)-DA*np.sin(theta1)

z=(np.power(x,2) + np.power(y,2) + CD**2 - BC**2)/(2 * CD)

w=z/(np.sqrt(np.power(x,2) + np.power(y,2)))
theta4=np.arccos(-w) + np.arctan(y/x) 

y_C=CD*np.sin(theta4)
y_B=AB*np.sin(theta2)

x_C=CD*np.cos(theta4)+DA
x_B=AB*np.cos(theta2)

theta3=np.arctan(y_C - y_B,x_C - x_B)

#plotting
pt.plot(teval,theta2,color='red')
pt.show()
pt.plot(teval,theta4,color='blue')
pt.show()