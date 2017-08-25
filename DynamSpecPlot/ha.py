import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as ss

folds = np.loadtxt("files.dat", dtype=bytes, delimiter='\t').astype(str)

order_str = "_008.norm.dat"

specs=[]
corrs=[]
phases=[]
for i in range(0, len(folds[:,0])):
    specs.append(np.loadtxt("../EXTRACTED/"+folds[i,0]+"/"+folds[i,1]+"/"+folds[i,1]+order_str))
    corrs.append(float(folds[i,3]))
    phases.append(float(folds[i,2]))
    
ews = np.zeros(len(folds))

v0=-11.8
for k in range(0, len(folds)):
    print(phases[k])
    wing=1.2
    xx=[]
    yy=[]
    x0=6562.8*(1.0+(v0-corrs[k])/300000.0)
    xlim1=x0-wing; xlim2=x0+wing
    print(xlim1, xlim2)
    for i in range(0, len(specs[k][:,0])):
        if((specs[k][i,0]>=xlim1) & (specs[k][i,0]<=xlim2)):
            xx.append(specs[k][i,0])
            yy.append(specs[k][i,1])
    yyf = ss.savgol_filter(yy, 5, 2)

    #plt.plot(xx, yyf)
    #plt.show()
    
    s=0
    for i in range(0, len(xx)-1):
        s=s+0.5*(yyf[i]+yyf[i+1])*(xx[i+1]-xx[i])
    ews[k]=s
   
    
print('Ave = ', np.mean(ews))
plt.plot(phases, ews, 'r+', markersize=16)
plt.xlabel("$\\varphi$")
plt.ylabel("EW H$\\alpha$")
plt.show()