#!/usr/bin/env python2
# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import sys
import os
import shutil

# sys.path.insert(1, os.path.abspath('../'))
# sys.path.insert(1, os.path.abspath('../functions'))


from ExoLens import ExoLens
import functions_main as fm
from functions_plotting import plot_cmf


"""
For exercise #1
"""

# Calculates the CMF and its uncertainties (-sigma, +sigma) from Mass and Radius
#e.g.  mass = [1, 0.01] and radius = [1., 0.01] for Mpl = 1+/-0.1 and Rpl = 1+/-0.1
mass = [1, 0.01]; radius = [1., 0.01]

CMF_planet = fm.calc_cmfrho(mass,radius)
print (CMF_planet)


"""
For exercise #2
"""

# Calcute CMF star from Fe/Mg and Si/Mg abundance ratios by mole (for Earth: Si/Mg = 0.9, Fe/Mg = 0.9)
#
FeMg = 0.9; sigFeMg = 0.01; SiMg = 0.9; sigSiMg = 0.1
CMF_star = fm.calc_cmf_star(FeMg, sigFeMg, SiMg, sigSiMg)
print (CMF_star)


"""
For exercise #3
"""

# Based on the previous estimates of CMF_planet and CMF_star, calcute the probability that the CMF planet is consistent with the CMF star. 
# 

Prob = fm.calc_integral(CMF_planet, CMF_star)

print (Prob)

#Plot the distributions
planet_name = "my_planet"
plot_cmf(mass, radius, planet_name, CMF_star)
os.system('move ../Outputs/Plots/%s.png ./' % planet_name)