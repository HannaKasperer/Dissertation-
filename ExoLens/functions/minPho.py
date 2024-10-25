# -*- coding: utf-8 -*-
"""
Created on Thu May 28 14:51:51 2020

@author: schulze.61
"""

import numpy as np
import matplotlib.pyplot as plt

from functions_main import calc_integral

from functions_main import calc_cmf_star


CMFrho = [26.0, 4.9, 5.1]
FeMg = 0.71; SiMg = 0.90


sigabund = np.linspace(0.06, 0.2, 200)

pho = []

for i in range(0, len(sigabund)):
    
    cmfstar = calc_cmf_star(FeMg, sigabund[i]*FeMg, SiMg, sigabund[i]*SiMg)
    integ = calc_integral(CMFrho, [cmfstar[0], cmfstar[1]])
    pho = np.append(pho, integ)
    
    
    
idx = np.argmin(abs(pho-5.0))

plt.axvline(sigabund[idx])
plt.axhline(pho[idx])
    
plt.plot(sigabund, pho)
plt.show()

cmfstar2 = calc_cmf_star(FeMg, sigabund[idx]*FeMg, SiMg, sigabund[idx]*SiMg)

cmfstarbc = calc_cmf_star(FeMg, 0.06*FeMg, SiMg, 0.06*SiMg)

print ("PHo: ", pho[idx])
print ("sigabund: ", sigabund[idx])
print ("CMFstar: ", cmfstar2)
print ("PHo dbl check: ", calc_integral(CMFrho, [cmfstar2[0], cmfstar2[1]]))

print("Best case Pho: ", calc_integral(CMFrho, [cmfstarbc[0], cmfstarbc[1]]))