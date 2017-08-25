# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 22:39:46 2016

@author: alex
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

data=np.loadtxt("input.dat")
sigma1=1
sigma2=2
phases = data[:,0]
vels1 = data[:,1]
vels2 = data[:,2]


def vrad(gamma, kamp, eks, omega, tper, t):
    v = vanom(eks, tper, t)
    return gamma+kamp*(eks*np.cos(omega)+np.cos(v+omega))

def vanom(eks, tper, t):
    eanom = eksanom(eks, tper, t)
    tg = np.sqrt((1+eks)/(1-eks))*np.tan(0.5*eanom)
    return 2*np.arctan(tg)
    
def eksanom(eks, tper, t):
    m=2*np.pi*(t-tper)
    ee =0
    for i in range(0, 150):
        ee = m+eks*np.sin(ee)
    return ee
    
def chi20(pars):
    res1=0
    for i in range(0, len(phases)):
        res1 = res1 + (vrad(pars[0], pars[1], pars[3], pars[4], pars[5], phases[i])-vels1[i])**2
    res2=0
    for i in range(0, len(phases)):
        res2 = res2 + (vrad(pars[0], -pars[2], pars[3], pars[4]+0.5*np.pi, pars[5]+0.5, phases[i])-vels2[i])**2
    res=res1/sigma1**2+res2/sigma2**2
    return res    


    
xx=[30, 20, 20, 0.5, 0.5, 0.5]

pars = minimize(chi20, xx, method="Nelder-Mead", tol=1e-15)
print(pars.x)

pp=np.linspace(0,1,50)
vv1=np.zeros(len(pp))
vv2=np.zeros(len(pp))

for i in range(0,len(pp)):
    vv1[i]=vrad(pars.x[0], pars.x[1], pars.x[3], pars.x[4], pars.x[5], pp[i]+0.0)
    vv2[i]=vrad(pars.x[0], pars.x[2], pars.x[3], pars.x[4], pars.x[5], pp[i]+0.0)

dev1=np.zeros(len(phases))
dev2=np.zeros(len(phases))
for i in range(0, len(phases)):
    dev1[i]=vels1[i]-vrad(pars.x[0], pars.x[1], pars.x[3], pars.x[4], pars.x[5], phases[i]+0.0)
    dev2[i]=vels2[i]-vrad(pars.x[0], pars.x[2], pars.x[3], pars.x[4], pars.x[5], phases[i]+0.0)


f, aa = plt.subplots(2,1, gridspec_kw = {'height_ratios':[2, 1]})
f.set_size_inches(4,7)
aa[1].axes.set_xlim([0,1])
aa[1].axes.set_xlim([0,1])
aa[0].set_xlabel("$\\varphi$")
aa[1].set_xlabel("$\\varphi$")
aa[0].set_ylabel('RV, km/s')
aa[1].set_ylabel('RV - model, km/s')
aa[0].plot(phases, vels1, '.r', markersize=10)
aa[0].plot(phases, vels2, '.b', markersize=10)
aa[1].plot(phases, dev1, '.r', markersize=10)
aa[1].plot(phases, dev2, '.b', markersize=10)
aa[0].plot(pp, vv1, '-k')
aa[0].plot(pp, vv2, '-k')
aa[0].grid(True)
aa[1].grid(True)

# Monte-Carlo
#print('Starting Monte-Carlo simulations')
#nn=200
#std1=np.std(dev1)
#std2=np.std(dev2)
#print('std1 = ', std1)
#print('std2 = ', std2)
#xx=[pars.x[0], pars.x[1], pars.x[2], pars.x[3]]
#yy01=np.zeros(len(vels1))
#yy02=np.zeros(len(vels2))
#for i in range(0,len(vels1)):
#    yy01[i]=vrad(xx[0], xx[1], xx[3], phases[i]+0.0)
#for i in range(0,len(vels2)):
#    yy02[i]=vrad(xx[0], xx[2], xx[3], phases[i]+0.0)
#spars=np.zeros(len(xx))
#for i in range(0,nn):    
#    for j in range(0,len(vels1)):
#        vels1[j]=yy01[j]+np.random.normal(0, std1)
#    for j in range(0,len(vels2)):
#        vels2[j]=yy02[j]+np.random.normal(0, std2)
#    pars = minimize(chi20, xx, method="Nelder-Mead", tol=1e-15)
#    for j in range(0, len(spars)):
#        spars[j]=spars[j]+(pars.x[j]-xx[j])**2

#for i in range(0, len(spars)):
#    spars[i]=np.sqrt(spars[i]/nn)
    
#print('dV0 = ', spars[0])
#print('dK1 = ', spars[1])
#print('dK2 = ', spars[2])
#print('dt0 = ', spars[3])