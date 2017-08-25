import matplotlib.pyplot as plt
import numpy as np
from phys_const import phys_const as pc 
from grid import abs_coeffs

def intensity(coeff_x, coeff_o, te, mag_str, size_par, lamb):
    mag_str = mag_str*1e6
    omega_c = pc.e * mag_str / (pc.me * pc.c)
    nu_rat = (2*np.pi* pc.c/(lamb*1e-8))/omega_c
    factRJ = pc.k * te / (8 * np.pi**3 * pc.c**2)
    io = factRJ * ((nu_rat * omega_c)**2) * (1.0 - np.exp(-coeff_o * size_par))
    ix = factRJ * ((nu_rat * omega_c)**2) * (1.0 - np.exp(-coeff_x * size_par))
    return (io + ix)/lamb**2
    
lambs = np.linspace(4000, 7000, 50)
intes = np.zeros(len(lambs))

obj = abs_coeffs()
path = 'output1.txt'
obj.openfile(path)

te_num = 6
mag_str = 40
theta = 55.68775*np.pi/180
log_size_par = 3.013535

size_par = 10**log_size_par
omega_c = mag_str*1e6*pc.e/ (pc.me*pc.c)
for i in range(0, len(lambs)):
    omega = 2*np.pi*pc.c/(lambs[i]*1e-8)
    nu_rat = omega/omega_c
    coeff_x = obj.interp('x', te_num, theta, nu_rat)
    coeff_o = obj.interp('o', te_num, theta, nu_rat)
    intes[i] = intensity(np.exp(coeff_x), np.exp(coeff_o), obj.tes[te_num], mag_str, size_par, lambs[i])  
plt.plot(lambs, intes, '-')
plt.show()