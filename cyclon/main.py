import matplotlib.pyplot as plt
import numpy as np
from phys_const import phys_const as pc 
from grid import abs_coeffs
from scipy.optimize import minimize

def intensity(coeff_x, coeff_o, te, mag_str, size_par, lamb):
    mag_str = mag_str*1e6
    omega_c = pc.e * mag_str / (pc.me * pc.c)
    nu_rat = (2*np.pi* pc.c/(lamb*1e-8))/omega_c
    factRJ = pc.k * te / (8 * np.pi**3 * pc.c**2)
    io = factRJ * ((nu_rat * omega_c)**2) * (1.0 - np.exp(-coeff_o * size_par))
    ix = factRJ * ((nu_rat * omega_c)**2) * (1.0 - np.exp(-coeff_x * size_par))
    return (io + ix)/lamb**2
     
def scale(mas1, mas2):
    sum1=0
    sum2=0
    for i in range(0, len(mas1)):
        sum1 = sum1+mas2[i]
        sum2= sum2+mas1[i]
    scale = sum2/sum1
    return scale
##########################################################################

spectrum = np.loadtxt("071126/s9120704.txt.red")
lambs = spectrum[:,0]
intes_obs =  spectrum[:,1]
intes = np.zeros(len(lambs))
mean = np.mean(intes_obs)
for i in range(0, len(lambs)):
    intes_obs[i] = intes_obs[i]/mean
te_num = 2
obj = abs_coeffs()
path = 'output2.txt'
obj.openfile(path)

def chi2(x):
    mag_str = x[0]
    theta = x[1]
    log_size_par = x[2]
    size_par = 10**(log_size_par)
    omega_c = mag_str*1e6*pc.e/ (pc.me*pc.c)
    for i in range(0, len(lambs)):
        omega = 2*np.pi*pc.c/(lambs[i]*1e-8)
        nu_rat = omega/omega_c
        coeff_x = obj.interp('x', te_num, theta, nu_rat)
        coeff_o = obj.interp('o', te_num, theta, nu_rat)
        intes[i] = intensity(np.exp(coeff_x), np.exp(coeff_o), obj.tes[te_num]*11.6e6, mag_str, size_par, lambs[i])
    mean = np.mean(intes)
    for i in range(0, len(lambs)):
        intes[i]=intes[i]/mean
    #s = scale(intes, intes_obs)
    chi = 0
    for i in range(0, len(lambs)):
        chi = chi + (intes[i] - intes_obs[i])**2
    #chi = chi / len(lambs)
    #print "[mag, theta, log_size] = ", x, ", CHI = ", chi
    return chi    

#mags = np.linspace(30, 40, 10)
#thetas = np.linspace(50*np.pi/180, 90*np.pi/180, 30)
#log_size = np.linspace(3,5, 10)
#chi_min=1e30
#mags_0=0
#thetas_0 =0
#log_size_0 =0
#for mm in mags:
#    for tt in thetas:
#        for ss in log_size:
#            chi = chi2([mm, tt, ss])
#            print mm, tt, ss, chi
#            if(chi<chi_min):
#                chi_min = chi
#                mags_0=mm
#                thetas_0 =tt
#                log_size_0 =ss

bnds = ((30,55),(40*np.pi/180, 89*np.pi/180),(2.5, 6.0))
x1mas = np.linspace(bnds[0][0], bnds[0][1], 3)
x2mas = np.linspace(bnds[1][0], bnds[1][1], 3)
x3mas = np.linspace(bnds[2][0], bnds[2][1], 3)
xinits = []
for i in x1mas:
    for j in x2mas:
        for k in x3mas:
           xinits.append([i,j,k])
x_mins = []
y_mins = []
i=0
for x in xinits:
    print "***************************************************************"
    print "TRY NUMBER ", i+1, "/", len(xinits) 
    print "INIT = ", x[0], x[1]*180/np.pi, x[2]    
    res = minimize(chi2, x, method='L-BFGS-B', bounds=bnds, tol = 1e-10)
    x_mins.append(res.x)
    y_mins.append(res.fun)
    print "X_MIN = ", np.round(x_mins[i][0],3), np.round(x_mins[i][1]*180/np.pi,3), \
    np.round(x_mins[i][2],3), "CHI = ", y_mins[i]
    i=i+1
x_min_glob = np.zeros(3)
y_min_glob = y_mins[0]
for i in range(0, len(xinits)):
    if(y_mins[i]<=y_min_glob):
        x_min_glob[0] = x_mins[i][0]
        x_min_glob[1] = x_mins[i][1]
        x_min_glob[2] = x_mins[i][2]
        y_min_glob = y_mins[i]
        
    

print "-----------------------------------------------------------------------------"
mags_0 = x_min_glob[0]
thetas_0 = x_min_glob[1]
log_size_0 = x_min_glob[2]               
print 'MAG = ', mags_0, ', THETA = ', thetas_0*180/np.pi, ', LOGLAM = ', log_size_0

size_par = 10**(log_size_0)
omega_c = mags_0*1e6*pc.e/ (pc.me*pc.c)
for i in range(0, len(lambs)):
    omega = 2*np.pi*pc.c/(lambs[i]*1e-8)
    nu_rat = omega/omega_c
    coeff_x = obj.interp('x', te_num, thetas_0, nu_rat)
    coeff_o = obj.interp('o', te_num, thetas_0, nu_rat)
    intes[i] = intensity(np.exp(coeff_x), np.exp(coeff_o), obj.tes[te_num], mags_0, size_par, lambs[i])
mean = np.mean(intes)
for i in range(0, len(lambs)):
    intes[i]=intes[i]/mean
plt.plot(lambs, intes)
plt.plot(lambs, intes_obs)
plt.show()