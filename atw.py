import numpy as np
import matplotlib.pyplot as plt

###############################################################################

alpha = 3e-3
Ho = 1.0

f  = 1e-4
g  = 9.8
r  = 1e-4
Lx = 1000e3
Ly = 100e3

Nx = 20000
Ny = 200

###############################################################################

x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
dx = x[1] - x[0]
dy = y[1] - y[0]

# Specify bathymetry
h = Ho + alpha*y

# Calculate bathymetric gradients
h_y = np.zeros_like(h)
h_y[1:-1] = (h[2:] - h[:-2]) / (2*dy)

# Specify wind forcing
X = 1e-4 * np.exp( -((x-100e3)/50e3)**2 )
# X = 1e-4 * np.sin(x/50e3)

# Initialize eta array
eta = np.zeros((Ny, Nx), 'd')

# set 'initial' upcoast condition at x=0
eta[:, 0] = np.zeros((Ny,), 'd')

fac1 = (dx / dy**2) * (r / (f * h_y[1:-1]))
fac2 = dy * f / (r * g)
for n in range(1,Nx):
    # The interior solution governed by heat equation
    eta[1:-1, n] = eta[1:-1, n-1] +  fac1 * np.diff(eta[:,n-1], n=2)
    # The coastal boundary condition
    eta[0, n] = eta[1, n] + fac2 * X[n]
    # The offshore-boundary condition
    # eta[-1, n] = eta[-2, n]
    eta[-1, n] = 0.0

# calculate the transport streamfunction
hu = 0.5 * (h[1:] + h[:-1])
psi = (g / f) * np.cumsum(hu[:, np.newaxis] * np.diff(eta, axis=0), axis=0)
psi = np.vstack((np.zeros((1, Nx)), psi))

###############################################################################

fig = plt.figure()
ax = fig.add_subplot(111)
pc = ax.pcolormesh(x/1e3, y/1e3, eta, cmap=plt.cm.RdBu_r)
eta0 = np.maximum(eta.max(), -eta.min())
pc.set_clim(-eta0, eta0)
cb = plt.colorbar(pc)
psi0 = np.maximum(psi.max(), -psi.min())
plt.contour(x/1e3, y/1e3, psi, np.linspace(-psi0, psi0, 9))
cb.set_label(r'$\eta$ [m]')
ax.set_xlabel('Along-shore direction [km]')
ax.set_ylabel('Cross-shore direction [km]')
plt.show()