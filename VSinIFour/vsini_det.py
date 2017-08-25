import numpy as np
from scipy.interpolate import LSQUnivariateSpline
import matplotlib.pyplot as plt
prof = np.loadtxt("aver_373.lsd")
xx = prof[:,0]
yy = prof[:,1]
plt.plot(xx, yy)
dx=xx[1]-xx[0]

t = np.linspace(-200, 200, 8)
spl = LSQUnivariateSpline(xx, yy, t)
xx_mod=np.linspace(-250, 250, 100)
yy_mod = spl(xx_mod)
plt.plot(xx_mod, yy_mod)
std=0
for i in range(0, len(xx)):
    std=std+(yy[i]-spl(xx[i]))**2
std=np.sqrt(std/len(xx))



print('std = ', std)
print('S/N = ', np.mean(yy)/std)

n_samples=50
yys = np.zeros([n_samples, len(yy)])
for i in range(0, n_samples):
    for j in range(0, len(xx)):
        yys[i,j]=spl(xx[j])+np.random.normal(0, std)
    #plt.plot(xx, yys[i,:])
    
plt.show()

yys1=np.zeros([n_samples, len(xx)+6000])
for k in range(0, n_samples):
    for i in range(0, len(yys1[k,:])):
        yys1[k,i]=1.0
    for i in range(0, len(xx)):
        yys1[k,i+3000]=yys[k, i]


nu = np.fft.fftfreq(len(yys1[0,:]), dx)
nu = np.fft.fftshift(nu)
min_ff=np.zeros(n_samples)
for k in range(0, n_samples):
    ff = np.fft.fft(yys1[k,:])/len(xx)
    ff = np.fft.fftshift(ff)
    ff = np.absolute(ff)**2
    plt.plot(nu, np.log(ff), '-')
    
    ff_lim_1=0.003
    ff_lim_2=0.005
    mini = 1e30
    mini_x = 0
    for i in range(0, len(nu)):
        if((nu[i]>=ff_lim_1)&(nu[i]<=ff_lim_2)):
            if(ff[i]<mini):
                mini=ff[i]
                mini_x=nu[i]
    
    min_ff[k]=mini_x
plt.xlim(0,0.01)
plt.show()
    
print(min_ff)
print('mean x_min = ', np.mean(min_ff))
print('std x_min = ', np.std(min_ff))

x_lin = 0.567
vsini=np.zeros(n_samples)
for i in range(0, n_samples):
    vsini[i] = (1/np.mean(min_ff[i]))*(0.610+0.062*x_lin+0.027*x_lin**2+0.012*x_lin**3+0.004*x_lin**4)

print('mean vsini = ', np.mean(vsini))
print('dvsini = ', np.std(vsini))