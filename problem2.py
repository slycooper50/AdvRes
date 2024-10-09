#!/bin/python
import numpy as np
import scipy.special as sc
import matplotlib.pyplot as plt

# Field Unit Conversion factors
mD = 9.869233E-13 * 1.E-3
ft = 0.3048
cp = 1.E-3
psi = 6894.76
hr = 3600.0
day = 86400.0
barrel = 0.158987
# Reservoir Properties (SI units)
k = 40 * mD
phi = 0.15
h = 10. 
r_w = 0.375 * ft
mu_o = 0.55E-3
c_o = 1.90E-9
c_w = 0.50E-9
So = 0.78
Sw = 0.22
c_f = 1.00E-9
c_t = c_o*So + c_w*Sw + c_f
Pi = 4.E7
q = - 185. * (barrel/day)

# Well flowing pressure
def P_wf(r, t):
    return Pi - (q*mu_o*1.2/(4*np.pi*k*h)) * sc.expi(-(phi*mu_o*c_t*r**2)/(4*k*t))
    #return Pi + (q*mu_o/(4*np.pi*k*h)) * np.log(4*k*t/(1.78*phi*mu_o*c_t*r**2))

# Plotting of pressure vs time prior to shut-in 
t_1 = 1. * day
t_10 = 10. * day
t_100 = 100. * day
r = np.arange(r_w, 1000.*ft+r_w, 1)

p_1 = P_wf(r, t_1)
p_2 = P_wf(r, t_10)
p_3 = P_wf(r, t_100)
fig, ax = plt.subplots()
ax.set_xlabel('radius (ft)')
ax.set_ylabel('Pressure (psi)')
l1, = ax.plot(r/ft, p_1/psi, color='r')
l2, = ax.plot(r/ft, p_2/psi, color='b')
l3, = ax.plot(r/ft, p_3/psi, color='g')
ax.set_ylim(top=Pi/psi)
ax.legend((l1,l2,l3), ('t = 1 day','t = 10 days','t = 100 days'), loc='lower right')
plt.show()
