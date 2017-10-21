import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from simplex import simplex as sim
from scipy.optimize import minimize

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
        X01*integrate.quad(lambda ksi: 2*phi1(ksi, delta, r2)*ksi, r2-delta, r2+delta)[0]+\
        X11*(2/3.0)*np.pi*(r2-delta)**2+X11*integrate.quad(lambda ksi: 2*ksi*np.sqrt(1-ksi**2/r1**2)*phi1(ksi, delta, r2), r2-delta, r2+delta)[0]
        
        s1=X01*np.pi*(r2-delta)**2
        s2=X01*integrate.quad(lambda ksi: 2*phi1(ksi, delta, r2)*ksi, r2-delta, r2+delta)[0]
        s3=X11*2*np.pi*integrate.quad(lambda ksi: np.sqrt(1-ksi**2/r1**2)*ksi, 0, r2-delta)[0]
        s4=X11*integrate.quad(lambda ksi: 2*ksi*np.sqrt(1-ksi**2/r1**2)*phi1(ksi, delta, r2), r2-delta, r2+delta)[0]
        res=s1+s2+s3+s4        
        print('0000', delta, res, s1, s2, s3, s4)        
        #if np.abs(delta)<1e-9:
        #    print('sssss')
        #    res=X01*np.pi*(r2-delta)**2+X11*(2/3.0)*np.pi*(r2-delta)**2
    return res

def GetFlux(inc, r1, r2, unecl, flux_rat, x1, x2, phase):
    theta=phase*2*np.pi    
    delta=Delta(inc, theta)
    
    intes1=np.pi*r1**2*(1-x1/3)*(1+1/flux_rat)
    intes1=1/intes1
    intes2=np.pi*r2**2*(1-x2/3)*(1+flux_rat)
    intes2=1/intes2
    
    X01=intes1*(1-x1)
    X11=intes1*x1
    X02=intes2*(1-x2)
    X12=intes2*x2    
    res=unecl
    if delta<r1+r2 and np.cos(theta)<0:
        res=res*(1-F1(inc, r1, r2, X01, X11, X02, X12, delta))
    if delta<r1+r2 and np.cos(theta)>=0:
        res=res*(1-F2(inc, r1, r2, X01, X11, X02, X12, delta))
    return res

def Chi2(x):
    global xx0, yy0
    global linx1, linx2
    s=0
    for i in range(len(xx0)):
        s=s+(yy0[i]-GetFlux(x[0], x[1], x[2], x[3], x[4], linx1, linx2, xx0[i]))**2
    return s

delta=0.125333233564
r2=0.2
r1=0.4
def f(ksi, delta, r2):
    if 1-ksi**2/r1**2 <0:
        print('hhhh')
    return 2*ksi*np.sqrt(1-ksi**2/r1**2)*phi1(ksi, delta, r2)
xxx=np.linspace(r2-delta, r2+delta, 100)
yyy=np.zeros(len(xxx))
for i in range(len(xxx)):
    yyy[i]=f(xxx[i], delta, r2)

#plt.plot(xxx, yyy)
#plt.show()

    
#data=np.loadtxt('LC_VZHya.dat', skiprows=2)
#xx0=data[:,0]; yy0=data[:,1]
#yy0 = 10**(-0.4*yy0)
linx1=0.5
linx2=0.5

#hx, hy = sim.NMopt(Chi2, [85*np.pi/180, 0.1, 0.1, 0.224, 1], 0.01, 1.0, 2.0, 0.5, 0.5, 500)
#hx=hx[-1,:]

#res = minimize(Chi2, [85*np.pi/2, 0.1, 0.1, 0.224, 1], method='nelder-mead')
#hx=res.x

#print('inc   = ', hx[0]*180/np.pi - np.floor(hx[0]*180/np.pi/360)*360)
#print('r1/a  = ', hx[1])
#print('r2/a  = ', hx[2])
#print('L0    = ', hx[3])
#print('L1/L2 = ', hx[4])

#xx1=np.linspace(0,1,1000)
#yy1=np.zeros(len(xx1))
#for i in range(len(xx1)):
#    yy1[i]=GetFlux(hx[0], hx[1], hx[2], hx[3], hx[4], linx1, linx2, xx1[i])
    
 
    
    
xx1=np.arange(-0.27,0.77,0.01)
yy1=np.zeros(len(xx1))
for i in range(len(xx1)):
    yy1[i]=GetFlux(1*np.pi/2, 0.2, 0.01, 8, 0.1, linx1, linx2, xx1[i])
    
#plt.plot(xx0, yy0, '.b')
plt.plot(xx1, yy1, '-r')
plt.show()