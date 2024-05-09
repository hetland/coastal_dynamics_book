
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt


x = np.linspace(0, 1, 50)
u_ini = np.ones_like(x)
dx = np.diff(x).mean()

def dudx(u, dx):
    dudx = np.empty_like(u)
    dudx[1:-1] = (u[2:]-u[:-2])/(2*dx)
    dudx[0] = (u[1]-u[0])/dx
    dudx[-1] = (u[-1]-u[-2])/dx
    return dudx

def cost(u):
    I1 = np.trapz( dudx(u, dx)**2 * dx, x )
    I2 = np.trapz( u**2 * dx, x )
    Iortho = 0.0
    for us in u_soln:
        Iortho += np.abs(np.trapz( u * us * dx, x ))
    return I1/I2 + 100.0 * (u[0] + (1.0 - np.trapz( u**2 * dx, x )) + Iortho)

u_soln = []
u1 = optimize.fmin_bfgs(cost, u_ini)
u_soln.append(u1)
u2 = optimize.fmin_bfgs(cost, u_ini)


plt.figure()
plt.plot(x, u1)
plt.plot(x, np.sin(x*np.pi/2.0), '--r')
plt.plot(x, u2, '-k')
plt.plot(x, -np.sin(x*3.0*np.pi/2.0), '--r')
plt.show()
