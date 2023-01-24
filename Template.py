# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 07:12:19 2021

@author: TROTABAS Baptiste
"""

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
#rc('xtick', labelsize=20) 
#rc('ytick', labelsize=20) 
#rc({'font.size': 22})

VERY_SMALL_SIZE = 7
SMALL_SIZE = 12
MEDIUM_SIZE = 16 #14
BIGGER_SIZE = 16
BIG_SIZE = 20
#-->
AIP_Text = 10 #8-10
AIP_LINEWIDTH = 1.5
rc('font', size=AIP_Text)          # controls default text sizes
rc('axes', titlesize=AIP_Text)     # fontsize of the axes title
rc('axes', labelsize=AIP_Text)     # fontsize of the x and y labels
rc('xtick', labelsize=AIP_Text)    # fontsize of the tick labels
rc('ytick', labelsize=AIP_Text)    # fontsize of the tick labels
rc('legend', fontsize=AIP_Text)    # legend fontsize
rc('figure', titlesize=AIP_Text, figsize = (1.25*3.37,3.0))   # fontsize of the figure title
rc('lines', linewidth = AIP_LINEWIDTH)