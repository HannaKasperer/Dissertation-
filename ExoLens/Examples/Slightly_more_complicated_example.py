#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 13:53:52 2020

@author: joeschulze

This example runs exolens on a sample with of likely rocky exos around stars 
where the data file only has metallicity measurements. This script uses some
of the miscellaneous functions to first approximate [Si/H] and [Mg/H] from 
[Fe/H] and linear fits to Hypatia [Mg/H] vs [Fe/H] and [Si/H] vs [Fe/H] data.
It then uses these approximated values to estimate the major refractory
abundances, Fe/Mg and Si/Mg, when then are used to calc CMF star.
"""

import pandas as pd
import numpy as np

import os
import sys


sys.path.insert(1, os.path.abspath('../'))
sys.path.insert(1, os.path.abspath('../functions'))

from ExoLens import ExoLens
import functions_misc as fm
import functions_plotting as fplt


data = pd.read_csv('../DataSets/Dai_USP_sample.csv')

#initialize arrays
cmfrho_arr = np.zeros(len(data)); cmfrho_uperr_arr = np.zeros(len(data)); cmfrho_lowerr_arr = np.zeros(len(data))
cmfstar_arr = np.zeros(len(data)); sigcmfstar_arr = np.zeros(len(data));
integral_arr = np.zeros(len(data))

for i in range(0, len(data)):
    
    mass = [data['Mass'][i], data['sigM'][i]]
    radius = [data['Radius'][i], data['sigR'][i]]
    
    FeH = data['FeH'][i]; sigFeH = data['sigFeH'][i]
    
    #In this example the input file only contains stellar metallicity values. 
    #In this case, we will use [Fe/H] and trends observed in the Hypatia Catalogue 
    #to approximate CMF star.
    SiH, sigSiH, MgH, sigMgH = fm.approx_SiH_MgH(FeH, sigFeH)
    FeMg, sigFeMg, SiMg, sigSiMg = fm.calc_FeMg_SiMg(FeH, sigFeH, SiH,  sigSiH, MgH, sigMgH)
    sratios = [FeMg, sigFeMg, SiMg, sigSiMg]
    
    
    #If you do not wish to generate plots, set the optional plot_flag to zero.
    #Default is to do the plotting.
    cmf_rho, cmf_star, P_ho = ExoLens(mass, radius, data['Planet'][i], sratios, plot_flag = 1)
    
    
    #append to output arrays. If no CMF star is provided to ExoLens, the integral
    #array is given a "nan".
    cmfrho_arr[i] = cmf_rho[0]; cmfrho_uperr_arr[i] = cmf_rho[1]; cmfrho_lowerr_arr[i] = cmf_rho[2]
    cmfstar_arr[i] = cmf_star[0]; sigcmfstar_arr[i] = cmf_star[1];
    integral_arr[i] =  P_ho
    


out = pd.DataFrame({'Planet': data['Planet'],'cmfrho': cmfrho_arr, 'sigcmfp': cmfrho_uperr_arr, 'sigcmfm': cmfrho_lowerr_arr, 'cmfstar': cmfstar_arr, 'sigcmfstar': sigcmfstar_arr, 'P(H0)':integral_arr})
out.to_csv('../Outputs/Results/PlanetSample_ExoLens_somewhat_more_complicated_example_results.csv')   
    