import numpy as np
data = np.loadtxt('LC.txt', skiprows=3)
data1=data[0:len(data[:,0]):3, :]
np.savetxt('out.dat', data1)