#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 18:33:03 2022

@author: trotabas
"""


import numpy             as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
from scipy import interpolate
# =============================================================================
# ---------------------: INPUT FILES
# =============================================================================
datContent = [i.strip().split() for i in open("input.dat").readlines()]
dict_input = {}
for i in open("input.dat").readlines():
    x = i.strip().split() 
    if x[0][0] != '/':
        dict_input[x[0]] = np.float64(x[1])
#===================================================================================
#---------------------:PHYSICAL-CONSTANT--------------------------------------------
#===================================================================================
r_earth     = dict_input.get("r_earth")        #-------Avogadro----Number---------[mole^-1]
G           = dict_input.get("G")              #-------Avogadro----Number---------[mole^-1]
m_earth     = dict_input.get("m_earth")        #-------Avogadro----Number---------[mole^-1]
mu          = G*m_earth;                       #-------Avogadro----Number---------[mole^-1]
g           = dict_input.get("g")              #-------Avogadro----Number---------[mole^-1]
#===================================================================================
#---------------------:SATELLITE-PARAMETERS-------------------------------
#===================================================================================
m           = dict_input.get("m")              #-------Avogadro----Number---------[mole^-1]
S           = dict_input.get("S")              #-------Avogadro----Number---------[mole^-1]    
#===================================================================================
#---------------------:THRUSTER-PARAMETERS-------------------------------
#===================================================================================
m_erg0      = dict_input.get("m_erg0")         #-------Avogadro----Number---------[mole^-1]
ISP         = dict_input.get("ISP")            #-------Avogadro----Number---------[mole^-1]
Ve          = ISP*g;                           #-------Avogadro----Number---------[mole^-1]
Tmax        = dict_input.get("Tmax")           #-------Avogadro----Number---------[mole^-1]
#===================================================================================
#---------------------:ORBIT-PARAMETERS----------------------------------
#===================================================================================    
H           = dict_input.get("H")              #-------Avogadro----Number---------[mole^-1]
r_min       = dict_input.get("r_min")             #-------Avogadro----Number---------[mole^-1]
r_max       = dict_input.get("r_max")             #-------Avogadro----Number---------[mole^-1]
#===================================================================================
#---------------------:DRAG-PARAMETERS----------------------------------
#===================================================================================    
Cd           = dict_input.get("Cd")              #-------Avogadro----Number---------[mole^-1]
#===================================================================================
#---------------------:TIME-PARAMETERS----------------------------------------------
#===================================================================================
t_max       = dict_input.get("t_max")          #-------Avogadro----Number---------[mole^-1]
dt          = dict_input.get("dt")             #-------Avogadro----Number---------[mole^-1]
dt_save    = dict_input.get("dt_save")       #-------Avogadro----Number---------[mole^-1]
#===================================================================================
#---------------------:TIME-PARAMETERS----------------------------------------------
#===================================================================================
Scenario_1     = int(dict_input.get("Scenario_1"))     #-------Keep-Orbit
Scenario_2     = int(dict_input.get("Scenario_2"))     #-------Reach new orbit

# %%
# =============================================================================
# Evolution of the air mass density 
# =============================================================================
path = "Output/"
df = pd.read_csv(path + 'Test_01.csv', delim_whitespace=True)#, header=None)#, delimiter=None)

data_h      = 1e-3*np.array(df['r(m)'].tolist())        # km
data_Vr     = np.array(df['Vr(m/s)'].tolist())          # m/s
data_time   = np.array(df['t(s)'].tolist())             # s
data_merg   = 1e3*np.array(df['m_erg(kg)'].tolist())    # g
data_theta  = np.array(df['theta(rad)'].tolist())
data_Omega  = np.array(df['Omega(s-1)'].tolist())
data_T      = np.array(df['T(mN)'].tolist())





# %%
fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize=(2*1.25*3.37,2*3), constrained_layout=False)
plt.gcf().subplots_adjust(right = 0.90, left = 0.095, wspace = 0.55, top = 0.99, hspace = 0.35, bottom = 0.075 )
ax00 = axes[0,0]
ax01 = axes[1,0]
ax10 = axes[0,1]
ax11 = axes[1,1]
# =============================================================================
# AXE 00
# =============================================================================
#-->
x=data_h*np.cos(data_theta)*1e3
y=data_h*np.sin(data_theta)*1e3
#-->
x_earth=r_earth*np.cos(np.linspace(0,2.0*np.pi,100))
y_earth=r_earth*np.sin(np.linspace(0,2.0*np.pi,100))
#-->
ax00.plot(x/r_earth, y/r_earth, color = "forestgreen", linestyle = "None", marker = ".", label = r'$\textrm{Trajectory polar}$')
ax00.plot(x_earth/r_earth, y_earth/r_earth, color = "blue", label = r'$\textrm{Earth surface}$')
ax00.plot(x[0]/r_earth, y[0]/r_earth,'d', color = "red", label = r'$\textrm{Start}$')
ax00.plot(x[-1]/r_earth, y[-1]/r_earth,'d', color = "darkorange", label = r'$\textrm{Stop}$')
ax00.legend(loc=0, frameon = False)
ax00.set_xlabel(r'$x/r_{earth}$')
ax00.set_ylabel(r'$y/r_{earth}$')
# =============================================================================
# AXE 10
# =============================================================================
ax10.plot(data_time/(24*3600), data_h-r_earth*1e-3, color = "forestgreen", label = r'$h$ $(\mathrm{km})$')
ax10.legend(loc=0, frameon = False)
ax10.set_xlabel(r'$\textrm{Day in orbit}$') 
ax10.set_ylabel(r'$\textrm{Altitude evolution}$ $(\mathrm{km})$') 
ax2 = ax10.twinx()
ax2.plot(data_time/(24*3600), data_merg, color = "crimson", label = r'$h$ $(\mathrm{km})$')
ax2.set_ylabel(r'$\textrm{Ergol mass evolution}$ $(\mathrm{g})$', color = "crimson") 
    
# =============================================================================
# AXE 01
# =============================================================================
V_norm = np.sqrt(data_Vr**2 + (data_h*1e3*data_Omega)**2)

ax01.plot(data_time/(24*3600), V_norm, color = "forestgreen", label = r'$\sqrt{ {V_r}^2 + \left ( r \Omega \right )^2}$ $(\mathrm{m/s})$')
ax01.legend(loc=0, frameon = False)
ax01.set_xlabel(r'$\textrm{Day in orbit}$') 
ax01.set_ylabel(r'$\sqrt{ {V_r}^2 + \left ( r \Omega \right )^2}$ $(\mathrm{m/s})$') 
ax2 = ax01.twinx()
ax2.plot(data_time/(24*3600), data_merg, color = "crimson", label = r'$h$ $(\mathrm{km})$')
ax2.set_ylabel(r'$\textrm{Ergol mass evolution}$ $(\mathrm{g})$', color = "crimson") 

# =============================================================================
# AXE 11
# =============================================================================
ax11.plot(data_time/(24*3600), data_T, color = "forestgreen", label = r'$T$ $(\mathrm{mN})$')
ax11.legend(loc=0, frameon = False)
ax11.set_xlabel(r'$\textrm{Day in orbit}$') 
ax11.set_ylabel(r'$\textrm{Trust evolution}$ $(\mathrm{mN})$') 
# ax11.set_yscale('log')
ax2 = ax11.twinx()
ax2.plot(data_time/(24*3600), data_merg, color = "crimson", label = r'$h$ $(\mathrm{km})$')
ax2.set_ylabel(r'$\textrm{Ergol mass evolution}$ $(\mathrm{g})$', color = "crimson") 