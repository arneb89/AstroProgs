import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

def psi1(ro, delta, r1):
    res=np.arccos((ro**2+delta**2-r1**2)/(2*ro*delta))
    return res

def phi1(ksi, delta, r2):
    return np.arccos((ksi**2+delta**2-r2**2)/(2*ksi*delta))

def Delta(inc, theta):
    return np.sqrt(np.cos(inc)**2+np.sin(inc)**2*np.sin(theta)**2)
    
def F1(inc, r1, r2, X01, X11, X02, X12, delta):
    res=0
    if (delta>=r1) and (delta>=(r1-r2)):
        res = X02*integrate.quad(lambda ro: 2*ro*psi1(ro, delta, r1), delta-r1, r2)[0]+\
        X12*integrate.quad(lambda ro: 2*ro*np.sqrt(1-ro**2/r2**2)*psi1(ro, delta, r1), delta-r1, r2)[0]
    if (delta<r1) and (delta>=(r1-r2)):
        res = X02*np.pi*(r1-delta)**2+ \
        X02*integrate.quad(lambda ro: 2*psi1(ro, delta, r1)*ro, r1-delta, r2)[0]+ \
        X12*(2/3.0)*np.pi*(r1-delta)**2+ \
        X12*integrate.quad(lambda ro: 2*ro*np.sqrt(1-ro**2/r2**2)*psi1(ro, delta, r1), r1-delta, r2)[0]
    if delta<(r1-r2):
        res = X02*np.pi*r2**2+X12*(2/3.0)*np.pi*r2**2
    return res
    
def F2(inc, r1, r2, X01, X11, X02, X12, delta):
    res=0
    if delta>=r2 and delta>=(r1-r2):
        res=X01*integrate.quad(lambda ksi: 2*ksi*phi1(ksi, delta, r2), delta-r2, r1)[0]+\
        X11*integrate.quad(lambda ksi: 2*ksi*np.sqrt(1-ksi**2/r1**2)*phi1(ksi, delta, r2), delta-r2, r1)[0]
    if delta>=r2 and delta<(r1-r2):
        print(delta)
        res=X01*integrate.quad(lambda ksi: 2*ksi*phi1(ksi, delta, r2), delta-r2, delta+r2)[0]+\
        X11*integrate.quad(lambda ksi: 2*ksi*np.sqrt(1-ksi**2/r1**2)*phi1(ksi, delta, r2), delta-r2, delta+r2)[0]
    if delta<r2 and delta>=(r1-r2):
        res=X01*np.pi*(r2-delta)**2+\
        X01*integrate.quad(lambda ksi: 2*phi1(ksi, delta, r2)*ksi, r2-delta, r1)[0]+\
        X11*(2/3.0)*np.pi*(r2-delta)**2+X11*integrate.quad(lambda ksi: 2*ksi*np.sqrt(1-ksi**2/r1**2)*phi1(ksi, delta, r2), r2-delta, r1)[0]
    if delta<r2 and delta<(r1-r2):
        res=X01*np.pi*(r2-delta)**2+\
        X01*integrate.quad(lambda ksi: 2*phi1(ksi, delta, r2)*ksi, r2-delta, r2+delta, epsabs=1e-40)[0]+\
        X11*(2/3.0)*np.pi*(r2-delta)**2+X11*integrate.quad(lambda ksi: 2*ksi*np.sqrt(1-ksi**2/r1**2)*phi1(ksi, delta, r2), r2-delta, r2+delta, epsabs=1e-40)[0]
    return res

def GetFlux(inc, r1, r2, intes1, intes_rat, x1, x2, phase):
    theta=phase*2*np.pi    
    delta=Delta(inc, theta)
    X01=intes1*(1-x1)
    X11=intes1*x1
    X02=(intes1/intes_rat)*(1-x2)
    X12=(intes1/intes_rat)*x2    
    res=np.pi*intes1*(1-x1/3)+np.pi*(intes1/intes_rat)*(1-x2/3)
    if delta<r1+r2 and np.cos(theta)>0:
        res=res*(1-F1(inc, r1, r2, X01, X11, X02, X12, delta))
    if delta<r1+r2 and np.cos(theta)<0:
        res=res*(1-F2(inc, r1, r2, X01, X11, X02, X12, delta))
    return res       
    

xx=np.linspace(0, 1, 100)
yy=np.zeros(len(xx))
for i in range(len(yy)):
    yy[i]=GetFlux(0.9991*np.pi/2, 0.7, 0.1, 1, 0.7, 0.6, 0.6, xx[i])
    
plt.plot(xx, yy, '.')