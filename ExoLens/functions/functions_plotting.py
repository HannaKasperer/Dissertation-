#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:05:06 2020

@author: joeschulze
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sp
from functions_main import cmf_rho
import matplotlib as mpl


#initialize some arrays.
cmfout = []
sigcmfp = []
sigcmfm = []
integ = []


angle = np.linspace(0.0,2*np.pi,100)
# outputpath = '../Outputs/Plots' 

def plot_cmf(mass, radius, planet_fname, cmfstar = []):
    
    #This function generates the planetary output plots.
 
    fig, axs = plt.subplots(2, 2, figsize = (20,20))
    axs = axs.flatten()
    
    cmfstar_flag = 0
        
    Mass = mass[0]; sigM = mass[1]
    Radius = radius[0]; sigR = radius[1]
    

    
    if len(cmfstar)>0:
        CMFstar = cmfstar[0]; sigCMFstar = cmfstar[1]
        cmfstar_flag = 1
    
    #This chunk of code generates CMF vs Radius plots for a given planet. The
    #plot consists of 10-90% confidence ellipses plotted as solid black lines,
    #and the 1-sigma (68%) and 2-sigma (95%) confidence ellipses plotted as dashed
    #cyan lines.
    for k in np.linspace(0,0.9,10):
        
        csquared = sp.chi2.ppf(k,2)
        rrrrx =(csquared**0.5)*sigR*np.cos(angle)+Radius
        rmrmx = (csquared**0.5)*sigM*np.sin(angle)+Mass
        cmf = cmf_rho(rmrmx, rrrrx)
        axs[0].plot(cmf,rrrrx, 'k-')
        
        csquared68 = sp.chi2.ppf(0.68,2)
        rrrr68 =(csquared68**0.5)*sigR*np.cos(angle)+Radius
        rmrm68 = (csquared68**0.5)*sigM*np.sin(angle)+Mass
        cmf68 = cmf_rho(rmrm68, rrrr68)
        axs[0].plot(cmf68,rrrr68, 'c--')
        
        csquared95 = sp.chi2.ppf(0.95,2)
        rrrr95 =(csquared95**0.5)*sigR*np.cos(angle)+Radius
        rmrm95 = (csquared95**0.5)*sigM*np.sin(angle)+Mass
        cmf95 = cmf_rho(rmrm95, rrrr95)
        axs[0].plot(cmf95,rrrr95, 'c--')
        
        
        if cmfstar_flag == 1:
            axs[0].axvline(CMFstar, color = 'm')
            axs[0].axvline(CMFstar+sigCMFstar, linestyle = '--', color = 'm')
            axs[0].axvline(CMFstar-sigCMFstar, linestyle = '--', color = 'm')

        axs[0].set_ylabel(r'Radius (R$_\oplus)$', fontsize = 32)
        axs[0].set_xlabel('CMF', fontsize = 32)
        axs[0].tick_params('both', direction = 'in', top = True, right = True, labelsize = 24, length = 10)
        axs[0].set_xlim([0.0,100])

        
    #------------------------------------------------------------------------------------------------#        
    #This chunk of code does the same thing as above only with CMF vs Mass    

        axs[1].plot(cmf,rmrmx, 'k-')
        axs[1].plot(cmf68,rmrm68, 'c--')
        axs[1].plot(cmf95,rmrm95, 'c--')

        if cmfstar_flag == 1:
            axs[1].axvline(CMFstar, color = 'm')
            axs[1].axvline(CMFstar+sigCMFstar, linestyle = '--', color = 'm')
            axs[1].axvline(CMFstar-sigCMFstar, linestyle = '--', color = 'm')
        axs[1].set_ylabel(r'Mass (M$_\oplus)$', fontsize = 32)
        axs[1].set_xlabel('CMF', fontsize = 32)
        axs[1].tick_params('both', direction = 'in', top = True, right = True, labelsize = 24, length = 10)
        axs[1].set_xlim([0.0,100.0])
   
       
    #------------------------------------------------------------------------------------------------#        
    #This chunk of code generates a contour map of CMF as a function of both mass and radius.
    #The M-R confidence intervals described above are overlain on this plot as well.

        massmat = np.linspace(Mass-3*sigM, Mass+3*sigM, 100)
        radmat = np.linspace(Radius-3*sigR, Radius+3*sigR, 100)
        RX, MY = np.meshgrid(radmat, massmat)
        CMFZ = cmf_rho(MY,RX)
        cmap = plt.get_cmap('coolwarm')
        norm = mpl.colors.Normalize(vmin=0.0,vmax=100.0)
    
    
        axs[2].contourf(RX,MY,CMFZ,  20, norm = norm, cmap = cmap)
        CS = axs[2].contour(RX,MY,CMFZ, 20, colors = 'w')

        axs[2].plot(rrrrx,rmrmx, 'k-', linewidth = 3)
        axs[2].clabel(CS, fontsize=12, inline=1)
        axs[2].plot(rrrr68,rmrm68, 'c--', linewidth = 3)
        axs[2].plot(rrrr95,rmrm95, 'c--', linewidth = 3)
    
        axs[2].set_ylabel(r'Mass (M$_\oplus$)', fontsize = 32)
        axs[2].set_xlabel(r'Radius (R$_\oplus$)', fontsize = 32)
        axs[2].tick_params('both', labelsize = 24, length = 10)
    
    
       
    
    
    #------------------------------------------------------------------------------------------------#
     
    #This chunck of code approximates the upper and lower 1-sigm uncertainties in
    #CMF rho from the 68% confidence ellipse and plots the resultant 1D PDF for CMF rho.
    #If CMF star is input, the expected CMF PDF if plotted as well.
    angle10 = np.linspace(0.0,2*np.pi,1000)
    rrrr =(csquared68**0.5)*sigR*np.cos(angle10)+Radius
    rmrm = (csquared68**0.5)*sigM*np.sin(angle10)+Mass
    cmfdif = cmf_rho(rmrm, rrrr) - cmf_rho(Mass, Radius)*np.ones(len(rrrr))


    sigmp = 0
    sigmm = 0
    sigrp = 0
    sigrm = 0

    countp = 0
    countm = 0

    for j in range(0, len(rrrr)):
        if cmfdif[j]>0:
            sigmp = sigmp + rmrm[j]
            sigrp = sigrp + rrrr[j]
            countp = countp+1
        if cmfdif[j]<0:
            sigmm = sigmm + rmrm[j]
            sigrm = sigrm + rrrr[j]
            countm = countm+1
        

        
    sigmm = sigmm/countm
    sigrm = sigrm/countm

    sigrp = sigrp/countp
    sigmp = sigmp/countp
    
    xarr = np.linspace(0,100, 1000)
    
    cmf = cmf_rho(Mass, Radius)
    
    
    if cmfstar_flag == 1:
        if cmf-CMFstar>0:
            sigcmf = abs(cmf - cmf_rho(sigmm,sigrm))
        
        else: 
            sigcmf = abs(cmf-cmf_rho(sigmp,sigrp))
        
    #if no CMF star value is given just use a nominal value of 35 (median FGK value)
    #to determine whether to use the upper or lower uncertainty in CMF rho
    #for plotting purposes
    else:
        if cmf-35.0>0:
            sigcmf = abs(cmf - cmf_rho(sigmm,sigrm))
        
        else: 
            sigcmf = abs(cmf-cmf_rho(sigmp,sigrp))

    
    planet = sp.norm.pdf(xarr, loc = cmf, scale = sigcmf)
    
    
    axs[3].plot(xarr, planet, 'k-', label = r'CMF$_\rho$')
    
    if cmfstar_flag == 1:
        star = sp.norm.pdf(xarr, loc = CMFstar, scale = sigCMFstar)
        axs[3].plot(xarr, star, 'm-', label = r'CMF$_\star$')

    
    axs[3].set_ylabel(r'$\phi (H_0)$', fontsize = 32)
    axs[3].set_xlabel('CMF', fontsize = 32)
    axs[3].tick_params('both', direction = 'in', top = True, right = True, labelsize = 24, length = 10)
    
    axs[3].legend(fontsize = 24)
    

    plt.suptitle(planet_fname, fontsize = 40)
    
    # plt.tight_layout()
    # plt.savefig(outputpath + '/' + planet_fname, facecolor = 'w' )
    # plt.close('all')
    plt.show()
    
  