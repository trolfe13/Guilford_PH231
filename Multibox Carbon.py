# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 16:03:36 2023

Multibox Carbon Cycle Model

@author: trolfe
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd

"""
Numerical Constants 
"""

B2 = 9.4
B3 = 10.2
a_d = 0.230
a_r = 1.0
gamma = 62.0    # PgC
GAMMA = 198     # PgC
k12 = 0.0931    # 1/yr
k13 = 0.0311    # 1/yr
k15 = 147       # 1/yr
k21 = 58*(730**(-B2))   # PgC^(1-B2) / yr
k23 = 0.0781    # 1/yr
k24 = 0.0164    # 1/yr
k31 = 18*(140**(-B3))   # PgC^(1-B3) / yr
k34 = 0.714     # 1/yr
k42 = 0.00189   # 1/yr
k43 = 0.00114   # 1/yr
k51 = 0.0862    # 1/yr
k56 = 0.0862    # 1/yr
k61 = 0.0333    # 1/yr

"""
Governing Equations 
"""

def Carbon_Model(X,t): #defines the system of diff. eq.
    
    # reforestation term
    def Fr(t):
        return 0
    
    # deforestation term
    def Fd(t):
        if t < 171: 
            return 0.3 + 0.01*t
        if t >= 171:
            return 0.3 + 0.01*171
    
    # fossil fuel emissions term
    def Ff(t):
        if t < 100:
            return 0.014*t
        if t >= 100 and t < 171:
            return 1.4 + (4.6/40)*(t-100)
        if t >= 171:
            return 0.0 # 1.4 + (4.6/40)*(171-100)
    
    # initial inputs
    M1, M2, M3, M4, M5, M6, M7, M8 = X
    
    
    # coupled differential equations
    
    eq1 = -(k12+k13)*M1 - k15*M8*(M1-gamma)/(M1+GAMMA) + k21*M2**B2 + k31*M3**B3 + k51*M5 + k61*M6 + Ff(t) + Fd(t) - Fr(t)
    
    eq2 = k12*M1 - (k23+k24)*M2 - k21*M2**B2 + k42*M4
    
    eq3 = k13*M1 + k23*M2 - k34*M3 - k31*M3**B3 + k43*M4
    
    eq4 = k24*M2 + k34*M3 - (k42+k43)*M4
    
    eq5 = k15*M8*(M1-gamma)/(M1+GAMMA) - (k51+k56)*M5 - Fd(t) + Fr(t)
    
    eq6 = k56*M5 - k61*M6
    
    eq7 = -Ff(t)
    
    eq8 = -(a_d*Fd(t) - a_r*Fr(t))/M5_0
    
    return (eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8)



"""
Initial Masses
"""

M1_0 = 612      # PgC
M2_0 = 730      # PgC
M3_0 = 140      # PgC
M4_0 = 37000    # PgC
M5_0 = 580      # PgC
M6_0 = 1500     # PgC
M7_0 = 5300     # PgC
M8_0 = 1        


"""
Spin-up Stuff
"""

# time-span between 1850 and 2020
t = np.arange(0,171,1)


# initial values of Carbon Cycle
u0 = np.array([M1_0, M2_0, M3_0, M4_0, M5_0, M6_0, M7_0, M8_0])


# call odeint to solve our equation
M = odeint(Carbon_Model,u0,t)

# import CO2 data
CO2_data = pd.read_csv(r'G:\My Drive\trolfe@bates.edu 2023-07-18 16 32\Bates College Courses\EA-PH 220 (Dynamical Climate)\Lab Stuff\CO2_1850.csv')

# Create the figure
fig, ax = plt.subplots()
txt = 'black'
p1, = ax.plot(t + 1850, M[:,0]/2.13, label = "Model")
p2, = ax.plot(CO2_data['Year'], CO2_data['PPM'], color = 'red', label = r'Atm CO$_2$')
ax.set_ylabel(r'Atmospheric Carbon (ppm)', color=txt, fontsize=16)
ax.set_xlabel(r'Time (yr)', color=txt, fontsize=16)
ax.set_title(r'Atmospheric Carbon vs Time', color=txt, fontsize = 18)
ax.tick_params(axis='x', labelsize=16, colors=txt)
ax.tick_params(axis='y', labelsize=16, colors=txt)
ax.legend(handles = [p1, p2], fontsize = 14)
#ax.set_ylim([0, 0.05])
plt.rcParams["figure.figsize"] = (9,9)
plt.show()


# Scenario, leave emissions steady at 2020 level

# initial values of Carbon Cycle
M1_1 = M[-1,0]
M2_1 = M[-1,1]
M3_1 = M[-1,2]
M4_1 = M[-1,3]
M5_1 = M[-1,4]
M6_1 = M[-1,5]
M7_1 = M[-1,6]
M8_1 = M[-1,7]

u1 = np.array([M1_1, M2_1, M3_1, M4_1, M5_1, M6_1, M7_1, M8_1])

def Carbon_Model_2(X,t): #defines the system of diff. eq.
    
    # reforestation term
    def Fr(t):
        return 0
    
    # deforestation term
    def Fd(t):
        return 0 #0.3 + 0.01*170
    
    # fossil fuel emissions term
    def Ff(t):
        return 0 #1.4 + (4.6/40)*(170-100)
    
    # initial inputs
    M1, M2, M3, M4, M5, M6, M7, M8 = X
    
    
    # coupled differential equations
    
    eq1 = -(k12+k13)*M1 - k15*M8*(M1-gamma)/(M1+GAMMA) + k21*M2**B2 + k31*M3**B3 + k51*M5 + k61*M6 + Ff(t) + Fd(t) - Fr(t)
    
    eq2 = k12*M1 - (k23+k24)*M2 - k21*M2**B2 + k42*M4
    
    eq3 = k13*M1 + k23*M2 - k34*M3 - k31*M3**B3 + k43*M4
    
    eq4 = k24*M2 + k34*M3 - (k42+k43)*M4
    
    eq5 = k15*M8*(M1-gamma)/(M1+GAMMA) - (k51+k56)*M5 - Fd(t) + Fr(t)
    
    eq6 = k56*M5 - k61*M6
    
    eq7 = -Ff(t)
    
    eq8 = -(a_d*Fd(t) - a_r*Fr(t))/M5_0
    
    return (eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8)

# time-span between 2020 and 2500
t1 = np.arange(0,981,1)

# call odeint to solve our equation
M1 = odeint(Carbon_Model_2,u1,t1)

# Create the figure
fig, ax = plt.subplots()
txt = 'black'
p1, = ax.plot(t + 1850, M[:,0]/2.13, label = "Model 1")
p2, = ax.plot(t1 + 2020, M1[:,0]/2.13, label = "Model 2")
ax.set_ylabel(r'Atmospheric Carbon (ppm)', color=txt, fontsize=16)
ax.set_xlabel(r'Time (yr)', color=txt, fontsize=16)
ax.set_title(r'Atmospheric Carbon vs Time', color=txt, fontsize = 18)
ax.tick_params(axis='x', labelsize=16, colors=txt)
ax.tick_params(axis='y', labelsize=16, colors=txt)
ax.legend(handles = [p1], fontsize = 14)
#ax.set_ylim([0, 0.05])
plt.rcParams["figure.figsize"] = (9,9)
plt.show()