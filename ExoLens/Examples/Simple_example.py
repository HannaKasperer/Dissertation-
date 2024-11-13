#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 13:50:44 2020

@author: joeschulze

This example script takes inputs from a csv file and computes the core mass
fraction (CMF) inferred from a planet's mass and radius, including estimates for the
upper and lower uncertainties on CMF. It also calculates the probability that
the CMF inferred from mass and radius are consistent with the CMF expected
from the major refractory abundances, Fe/Mg and Si/Mg, of the host star.

This script also generates a plot of the estimated planetary CMF's vs Radius
for the input sample.
"""


import pandas as pd
import numpy as np
import sys
import os

sys.path.insert(1, os.path.abspath('../'))
sys.path.insert(1, os.path.abspath('../functions'))

from ExoLens import ExoLens

data = pd.read_csv('../DataSets/PlanetSample.csv')

#initialize arrays
cmfrho_arr = np.zeros(len(data)); cmfrho_uperr_arr = np.zeros(len(data)); cmfrho_lowerr_arr = np.zeros(len(data))
cmfstar_arr = np.zeros(len(data)); sigcmfstar_arr = np.zeros(len(data));
integral_arr = np.zeros(len(data))

for i in range(0, len(data)):
    
    #get input parameters from input file
    mass = [data['Mass'][i], data['sigM'][i]]
    radius = [data['Radius'][i], data['sigR'][i]]
    sratios = [data['FeMg'][i], data['sigFeMg'][i], data['SiMg'][i], data['sigSiMg'][i]] 
    

    #If you do not wish to generate plots, set the optional plot_flag to zero.
    #Default is to do the plotting.
    cmf_rho, cmf_star, P_ho = ExoLens(mass, radius, data['Planet'][i], sratios, plot_flag = 1)
    
    
    #append to output arrays. If no abundances are provided to ExoLens, the integral
    #array is given a "nan".
    cmfrho_arr[i] = cmf_rho[0]; cmfrho_uperr_arr[i] = cmf_rho[1]; cmfrho_lowerr_arr[i] = cmf_rho[2]
    cmfstar_arr[i] = cmf_star[0]; sigcmfstar_arr[i] = cmf_star[1];
    integral_arr[i] =  P_ho
    

out = pd.DataFrame({'Planet': data['Planet'],'cmfrho': cmfrho_arr, 'sigcmfp': cmfrho_uperr_arr, 'sigcmfm': cmfrho_lowerr_arr, 'cmfstar': cmfstar_arr, 'sigcmfstar': sigcmfstar_arr, 'P(H0)':integral_arr})

out.to_csv('../Outputs/Results/PlanetSample_ExoLens_simple_example_results.csv')   
    




