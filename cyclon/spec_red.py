import numpy as np
import matplotlib.pyplot as plt


dir = "071126/"
file = "s9120706.txt"

bad_areas = [[2000, 4000], [6498,6599], [6650,6702], [5839, 5895], [4816, 4883], [6868, 6919], \
            [4450, 4486], [4650, 4713], [4312, 4371], [7000, 9000], [4076, 4115], [5383, 5434]]

spectra = np.loadtxt(dir+file)

x=[]
y=[]
for i in range(0, len(spectra[:,0])):
    b  = True
    for j in range(0, len(bad_areas)):
        if((spectra[i, 0]> bad_areas[j][0]) & (spectra[i, 0]< bad_areas[j][1])):
            b = False
    if(b == True):
        x.append(spectra[i,0])
        y.append(spectra[i,1])

z = np.polyfit(x, y, 25)
p = np.poly1d(z)
xbin = []
ybin = []
lam_start = 4000
lam_stop = 7000
index_start =0
index_stop = len(spectra[:,0])-1
for i in range(0, len(spectra[:,0])):
    if(lam_start<=spectra[i,0]):
        index_start = i
        break
for i in range(len(spectra[:,0])-1, 0, -1):
    if(lam_stop>=spectra[i,0]):
        index_stop = i
        break
i=index_start
win = 50
while(i<index_stop):
    mean_inten = 0
    mean_lambd = 0
    print(i)
    for j in range(i, i++win):
        
        b  = True
        for k in range(0, len(bad_areas)):
            if((spectra[j, 0]> bad_areas[k][0]) & (spectra[j, 0]< bad_areas[k][1])):
                b = False
        if(b == False):        
            mean_inten = mean_inten+p(spectra[j,0])
        else:
            mean_inten = mean_inten + spectra[j,1]
        mean_lambd = mean_lambd + spectra[j,0]
    mean_inten = mean_inten/ win
    mean_lambd = mean_lambd/ win
    xbin.append(mean_lambd)
    ybin.append(mean_inten)
    i=i+win

    
plt.plot(spectra[:,0], spectra[:,1])
plt.plot(xbin, ybin, 'r')
plt.plot(x, p(x), 'r.')
plt.ylim(0.2e-15, 0.2e-14)
plt.show()

out = open(dir+file+'.red', 'w')
for i in range(0, len(xbin)):
    out.write(str(xbin[i])+" "+ str(ybin[i])+ "\n")
out.close()