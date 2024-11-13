#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 16:13:37 2020

@author: joeschulze
"""

import numpy as np



#Function to calc Fe/Mg and Si/Mg from [Fe/H], [Mg/H], and [Si/H]. [X/H] should
#be given to this function in dex!
def calc_FeMg_SiMg(FeH, sigFeH, SiH,  sigSiH, MgH, sigMgH):
    #Ref solar values from Lodders 2009
    FeH_sun = 7.45; SiH_sun = 7.52; MgH_sun = 7.54
    #sigFeH_sun = 0.03; sigSiH_sun = 0.02; sigMgH_sun = 0.02
    sigFeH_sun = 0.08; sigSiH_sun = 0.06; sigMgH_sun = 0.06
    #sigFeH_sun = 0.0; sigSiH_sun = 0.0; sigMgH_sun = 0.0
        
    FeMg = 10**(FeH+FeH_sun-MgH-MgH_sun)
    SiMg = 10**(SiH+SiH_sun-MgH-MgH_sun)
    
    sigFeMg = FeMg*np.log(10.0)*np.sqrt(sigFeH**2 + sigMgH**2 + sigFeH_sun**2 + sigMgH_sun**2)
    sigSiMg = SiMg*np.log(10.0)*np.sqrt(sigSiH**2 + sigMgH**2 + sigSiH_sun**2 + sigMgH_sun**2)
    
    return FeMg, sigFeMg, SiMg, sigSiMg


#approximate [Si/H] and [Mg/H] from Hypatia Catalogue. Only have stellar metalicity?
#Try approximating [Si/H] and [Mg/H] from the Hypatia Catalogue!
def approx_SiH_MgH(FeH, sigFeH):
    SiH = 0.03492595 + 0.85403399*FeH
    MgH = 0.06313159 + 0.77379491*FeH
    
    sigSiH = ((0.85403399*sigFeH)**2 + 0.1344355**2)**0.5
    sigMgH = ((0.77379491*sigFeH)**2 + 0.121829**2)**0.5
    
    return SiH, sigSiH, MgH, sigMgH



