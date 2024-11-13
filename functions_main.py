#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 13:54:11 2020

@author: joeschulze
"""

import numpy as np
import scipy.stats as sp
from scipy.integrate import quad



#Define CMF function and density functions

#General CMF function in terms of planetary observables -- Mass and Radius. 
def cmf_rho(M, R):
    cmfrho = (4./3.)*np.pi*(R**3)*rhoc(M) - M*rhoc(M)/rhom(M)
    cmfrho = cmfrho/(M*(1-rhoc(M)/rhom(M)))
    return 100*cmfrho


#general density function 
def rho(M,a0,a1,a2):
    rho = a0+(a1*(M**a2))
    return rho


#Build functions for the bulk core density and bulk mantle density using the 
#general density function. The parameters being fed into the rho function were
#found from fits of rhoc vs Mp and rhom vs Mp using exoplex over the mass
#range of 0.01 Earth masses to 11 Earth masses.
#   rhoc: bulk core density
#   rhom: bulk mantle density
def rhoc(Mp):
    rhoc = rho(Mp, 0.36125593204992873, 0.153364551060784, 0.5351858792798839)
    return rhoc


def rhom(Mp):

    rhom = rho(Mp, 0.12372222237561939, 0.06565578408094341, 0.4096898631741826)
    return rhom



#calc_cmfrho mostly calculates the 1-sigma uncertainty in CMF rho. It does so
#by first generating the 68% confidence ellipse of the joint mass-radius distribution.
#It then feeds this ellipse into cmf_rho to generate the 68% confidence interval
#on CMF rho. Finally it averages all of the values along this ellipse that yield
#a CMF that is greater than the mean CMF value to give a scalar value for the 
#upper uncertainty in CMF rho. Similarly, it averages all of the CMF values along
#this ellipse that yield a value that is less than the mean CMF to calculate
#a scalar approximation for the lower uncertainty on CMF rho.
    
def calc_cmfrho(mass, radius):
    
    #unpack mass and radius
    Mass = mass[0]; sigM = mass[1]
    Radius = radius[0]; sigR = radius[1]
    
    #ellipse angle array
    angle = np.linspace(0.0,2*np.pi,1000)
    
    #68% confidence ellipse
    csquared68 = sp.chi2.ppf(0.68,2)
    
    
    #Calc semi-major and semi-minor axes.
    rrrr =(csquared68**0.5)*sigR*np.cos(angle)+Radius
    rmrm = (csquared68**0.5)*sigM*np.sin(angle)+Mass
    
    #calculate the difference of between all points on the CMF ellipse 
    #and the mean CMF value.
    cmfdif = cmf_rho(rmrm, rrrr) - cmf_rho(Mass, Radius)*np.ones(len(rrrr))

    #The next block of code takes the averages above and below the mean
    #CMF value to approximate the upper and lower uncertanties in CMF rho.
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
    
    #print 'M upper: ', sigmp
    #print 'R upper: ', sigrp
    
    #print 'M lower: ', sigmm
    #print 'R lower: ', sigrm
    
    
    cmfrho = cmf_rho(Mass, Radius)
    
    sigcmfp_val = abs(cmfrho-cmf_rho(sigmp,sigrp))
    sigcmfm_val = abs(cmfrho-cmf_rho(sigmm,sigrm))
    
    return cmfrho, sigcmfp_val, sigcmfm_val


#This function calculates the probability that CMF rho is consistent with CMF star
#under the assumption that the planet can be described with a purely rocky composition
#consisting of a pure Fe core and pure MgSiO3 mantle.
def calc_integral(cmfrho, cmfstar, xarr = np.linspace(0,100,10**6)):
    
    CMFrho = cmfrho[0]; sigCMFrho_p = cmfrho[1]; sigCMFrho_m = cmfrho[2]
    CMFstar = cmfstar[0]; sigCMFstar = cmfstar[1]

    
    if CMFrho-CMFstar>0:
        sigcmfrho = sigCMFrho_m
        
    else: 
        sigcmfrho = sigCMFrho_p

    
    def jointref(x):
        jr = sp.norm.pdf(x, loc = CMFstar, scale = sigcmfrho)*sp.norm.pdf(x, loc = CMFstar, scale = sigCMFstar)
        return jr
    
    def joint(x):
        j = sp.norm.pdf(x, loc = CMFrho, scale = sigcmfrho)*sp.norm.pdf(x, loc = CMFstar, scale = sigCMFstar)
        return j
    

    #Find Normalization constant. This constant comes from the assumption that the two distributions are the same (they aren't).
    #This will allow us to find how similar the distributions actually are relative to the ideal case that they are the same.

    normconst = quad(jointref, -200, 200)[0]
    integral = quad(joint, -200, 200)[0]

    integral = 100*integral/normconst

    
    return  integral


#calc radius as a function of CMF and planetary mass. Useful for making theoretical
#MR diagrams.
def calc_Rp(CMF, Mp):
    Rp = ((CMF*Mp*(1-rhoc(Mp)/rhom(Mp)) + Mp*rhoc(Mp)/rhom(Mp))*3.0/(4.0*np.pi*rhoc(Mp)))**(1.0/3.0)
    return Rp



#Calcute CMF star from Fe/Mg and Si/Mg.
def calc_cmf_star(FeMg, sigFeMg, SiMg, sigSiMg):
    mfe = 55.845; msio2 = 60.08; mmgo = 40.3044
    cmfstar = FeMg*mfe/(FeMg*mfe + SiMg*msio2 + mmgo)
    
    sigcmfstar = (FeMg*msio2*sigSiMg)**2 + (sigFeMg*(mmgo + msio2*SiMg))**2
    sigcmfstar = mfe*np.sqrt(sigcmfstar)/((FeMg*mfe + mmgo + msio2*SiMg)**2)
    
    return 100*cmfstar, 100*sigcmfstar

    