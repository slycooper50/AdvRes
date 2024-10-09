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
k = 200 * mD
phi = 0.20 
h = 30 * ft
r_w = 0.5 * ft
mu_o = 0.8 * cp
c_t = 3.0E-5 * (1./psi)
Pi = 5000. * psi
q = - 500. * (barrel/day)
tp = 3100.0 * hr
dt = 1. * hr

# Well flowing pressure
def P_wf(r, t):
    return Pi - (q*mu_o/(4*np.pi*k*h)) * sc.expi(-(phi*mu_o*c_t*r**2)/(4*k*t))
    #return Pi + (q*mu_o/(4*np.pi*k*h)) * np.log(4*k*t/(1.78*phi*mu_o*c_t*r**2))
# Well shut-in pressure
def P_ws(r, dt):
    return Pi - (q*mu_o/(4*np.pi*k*h)) * np.log((tp+dt)/dt)

# (a)
p1 = P_wf(r_w, tp)
# (b)
p2 = P_ws(r_w, dt)
print(p1/psi)
print(p2/psi)

# Plotting of pressure vs time prior to shut-in 
t = np.arange(0.001, tp+dt, dt)
p = P_wf(r_w, t)
fig, ax = plt.subplots()
ax.set_xlabel('time (hr)')
ax.set_ylabel('Pressure (psi)')
#ax.set_xscale('log', base=np.exp(1))
ax.plot(t/hr, p/psi, color='r')
ax.set_ylim(top=Pi/psi)
plt.show()
