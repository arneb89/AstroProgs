# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 20:49:37 2017

@author: alex
"""
import numpy as np

class simplex:
    
    def NMopt(func, x0, ll, alpha, gamma, ro, sigma, niter):    
        hh=[]
        hf=[]
        
        nn = len(x0)+1
        dd = len(x0)
        xx = []; yy = []
        for i in range(nn):
            xx.append(np.zeros(dd))
        for i in range(len(xx)):
            xx[i] = x0 + np.array([np.random.normal(0, ll) for i in range(nn-1)])
            yy.append(func(xx[i]))
            
        
        for k in range(niter):
            # Sorting        
            ss=sorted(list(zip(xx, yy)), key=lambda x: x[1])
            xx = np.array(ss)[:,0]
            yy = np.array(ss)[:,1]
            hh.append(xx[0])
            hf.append(yy[0])
            
            print('fmin = ', yy[0])
            # grav center
            x0 = np.mean(xx[0:len(xx)-1])
            # reflection
            xr = x0+alpha*(x0-xx[-1])
            yr = func(xr)
            if yr<yy[-2] and yr>=yy[0]:
                xx[-1] = xr        
                yy[-1] = yr
                print(k, 'reflection')
                continue
            if yr<yy[0]:
                xe = x0 + gamma*(xr-x0)
                ye = func(xe)
                if ye<yr:
                    xx[-1] = xe
                    yy[-1] = ye
                    print(k, 'expantion')
                    continue
                else:
                    xx[-1] = xr        
                    yy[-1] = yr
                    print(k, 'reflection')
                    continue
            else:
                xc = x0+ro*(xx[-1]-x0)
                yc = func(xc)
                if yc < yy[-1]:
                    xx[-1] = xc
                    yy[-1] = yc
                    print(k, 'contraction')
                    continue
                else:
                    for i in range(1, len(xx)):
                        xx[i] = xx[0] + sigma*(xx[i]-xx[0])
                    print(k, 'shrinkage')
                    continue
        return np.array(hh), np.array(hf)