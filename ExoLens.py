#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:27:08 2020

@author: joeschulze

This is the main ExoLens function. Most of the work is done in the subfunctions
found in ./functions/functions_main.py. Please see the read me file for more details.
"""

import functions_main as fm
from functions_plotting import plot_cmf



def ExoLens(mass, radius, planet_fname = ' ', sratios = [], plot_flag = 1):
    
    if len(sratios)>0:
        FeMg = sratios[0]; sigFeMg = sratios[1]
        SiMg = sratios[2]; sigSiMg = sratios[3] 
    
        #calc CMFstar
        cmfstar = fm.calc_cmf_star(FeMg, sigFeMg, SiMg, sigSiMg)
        
    else:
        cmfstar = []
    
    #Generate plots for each planet in sample if specified by user. Default
    #is to generate the plots.
    if plot_flag !=0:
        plot_cmf(mass, radius, planet_fname, cmfstar)
        
    
    #Calculate the core mass fraction inferred from planetary M and R
    cmfrho = fm.calc_cmfrho(mass, radius)
    
    #If a value for the expected CMF is given then calculate the probability 
    #that the CMF calculated from planetary M and R is consistent with the
    #expected value. If no expected value is given, return 'NaN'.
    if len(sratios)>0:
        integral = fm.calc_integral(cmfrho, cmfstar)
    else:
        integral = float("NaN")
        
 
    return cmfrho, cmfstar, integral
    
   
        
    
        
    
    