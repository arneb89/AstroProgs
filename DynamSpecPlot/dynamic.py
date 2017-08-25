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
    
xx=[]
yy=[]
v0=-11.8
lam0=6562.8
wing=4 # wing in angstrems
for k in range(0, len(folds)):
    xx1=[]
    yy1=[]
    x0=lam0*(1.0+(v0-corrs[k])/300000.0)
    xlim1=x0-wing; xlim2=x0+wing
    for i in range(0, len(specs[k][:,0])):
        if((specs[k][i,0]>=xlim1) & (specs[k][i,0]<=xlim2)):
            xx1.append(specs[k][i,0])
            yy1.append(specs[k][i,1])
    yyf = ss.savgol_filter(yy1, 5, 2)
    xx.append(xx1)
    yy.append(yyf)

for i in range(0, len(xx)):
    x0=lam0*(1.0+(v0-corrs[i])/300000.0)
    for j in range(0, len(xx[i])):
        xx[i][j] = (xx[i][j] - x0)*300000/x0

period = 1.76159 #days
expos = (30/60)/24/period
pp =  np.linspace(0,1,100)
vv = np.linspace(-100, 100, 200)
xxm, yym = np.meshgrid(vv, pp)
zz=np.zeros((len(pp),len(vv)))

for i in range(0,len(pp)):
    filled = False
    for p in range(0, len(phases)):
        if (phases[p]<=(pp[i]+expos/2)) & (phases[p]>(pp[i]-expos/2)):
            for k in range(0, len(vv)):
                zz[i,k]=np.interp(vv[k], xx[p], yy[p])
            filled = True
    if filled==False:
        for k in range(0, len(vv)):
                zz[i,k]=1
                
                

plt.imshow(zz, interpolation="nearest", cmap=plt.get_cmap("seismic"), extent=[vv[0],vv[len(vv)-1],1,0], aspect=300)
plt.xlabel("Velosity, km/s")
plt.ylabel("Phase")
plt.colorbar()
plt.show()