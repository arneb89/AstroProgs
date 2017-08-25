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
    return (io + ix)
     
##########################################################################

lambs=np.linspace(4000, 7000, 1000)
intes=np.zeros(len(lambs))
theta=52*np.pi/180
mag_str=51
te_num = 0
size_par=10**(2.5)
obj = abs_coeffs()
path = 'GRID_15_19_1.dat'
obj.openfile(path)

omega_c = mag_str*1e6*pc.e/ (pc.me*pc.c)
for i in range(0, len(lambs)):
        omega = 2*np.pi*pc.c/(lambs[i]*1e-8)
        nu_rat = omega/omega_c
        coeff_x = obj.interp('x', te_num, theta, nu_rat)
        coeff_o = obj.interp('o', te_num, theta, nu_rat)
        intes[i] = intensity(np.exp(coeff_x), np.exp(coeff_o), obj.tes[te_num]*11.6e6, mag_str, size_par, lambs[i])
		
plt.plot(lambs, intes, '-')
plt.show()