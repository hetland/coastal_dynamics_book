
import numpy as np
import matplotlib.pyplot as plt

alpha = 1e-3
Ho    = 1e-6  # maintain a small, positive depth at the coast 
            # to eliminate a singularity in the definition of b_i

Ndays = 10

f  = 1e-4
g  = 9.8
Lx = 2000e3
Ly = 100e3
T  = Ndays * 86400.0

Nx = 1000
Ny = 500
Nt = Ndays * 24



def calculate_cross_shore_modes(h, dy):
    '''
    Return cont. shelf wave cross-shore modal structures, phase speed, and weights
    
    Parameters
    ----------
    h   :   Array defining the depths (must be positive-definite everywhere)
    dy  :   Float defining the (uniform) cross-shore coordinate grid resolution
    
    Returns
    -------
    F   :   Array of cross-shore modes
    ci  :   Phase-speeds associated with each mode
    bi  :   Normalizing weights associated with each mode
    
    '''
    
    assert np.all(h > 0.0), 'Depth must be positive definite everywhere'
    
    # calculate the depth gradients.  Note, the two boundary values are not used.
    h_y= np.ones_like(h)
    h_y[1:-1] = (h[2:]-h[:-2])/(2*dy)
    
    # Geometric factors based on the bathymetry and Coriolis parameter
    fac1 = h / (f * h_y * dy**2)
    fac2 = h_y / (2 * h * dy)
    
    u = -fac1 + fac2    # upper (i+1) diagonal
    d = 2*fac1          # center (i) diagonal
    l = -fac1 - fac2    # lower (i-1) diagonal
    
    # construct tri-diagonal matrix representing a discrete representation
    # of the boundary value problem given by equation 2.12.14 only for the
    # interior of the domain.
    A = np.diag(u[1:-2], k=1) + np.diag(d[1:-1], k=0) + np.diag(l[2:-1], k=-1)
    
    # Apply the boundary condition that dF/dy=0 at y=L.
    # (Note, nothing else needs to be done for the F=0 at y=0 boundary condition.)
    A[-1, -1] += -fac1[-2] + fac2[-2]
    
    # Calculate the eigenvalues and eigenvectors, and sort by eigenvalue.
    lam, v = np.linalg.eig(A)
    idx = np.argsort(lam)
    ci = 1.0 / lam[idx]   # The phase speed is the inverse eigenvalue.
    
    # Add in the two boundary values to the eigenvectors.
    print v.shape
    F = np.vstack((np.zeros((1, v.shape[1])), v[:, idx], v[-1, idx]))
    
    # Add boundary values to the cross-shore coordinate, to match.
    # y = np.hstack((0, y, y[-1]+dy))
    bi = np.trapz(F[:].T * h_y / h**2, dx=dy)
    
    return F*np.sign(bi), ci, abs(bi)


# Create the cross-shore coordinate
x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
t = np.linspace(0, T, Nt)

# Define the bathymetry as a function of y array.
h1 = Ho + alpha * y             # Linear gradient for deep ocean
# h2 = Ho + 30.0*np.exp(y/30e3)    # Exponential profile for shelf/slope  
# h = np.minimum(h1, h2)
# Nsmooth = 5
# h[Nsmooth:-Nsmooth] = np.convolve(np.ones((2*Nsmooth+1,), 'd')/(2.0*Nsmooth+1.0), h, mode='valid')
h = h1

# Calculate the 
F, ci, bi = calculate_cross_shore_modes(h, dy=y[1]-y[0])

def display_modes(N, ax=None):
    if ax is None:
        ax=plt.gca()
    for mode in range(N):
        ax.plot(y, bi[mode] * F[:, mode])
        print 'mode ', mode, ', phase speed = ', ci[mode]

fig = plt.figure()

ax1 = fig.add_subplot(211)
display_modes(5, ax=ax1)

ax2 = fig.add_subplot(212)
ax2.plot(np.abs(bi), 'ko')

plt.show()


def X(x, t):
    X = 5e-3 * np.exp(-((x-200e3)/100e3)**2 -((t-2.0*86400.0)/(1.0*86400.0))**2)
    X = X * (t-2.0*86400.0)/(1.0*86400.0)
    return X


def forced_modes(Ntime, Nnodes=100):
    '''
    calcualte stream-function at t=to based on first Nnodes nodes
    '''
    # time over which to integrate the characteristic, i.e., all times up to t[Ntime]
    ts = t[:Ntime+1]              
    
    psi = np.zeros((Ny, Nx), 'd')
    for n in range(Nnodes):
        # calculate the (const) value of the characteristic for each gridpoint in x
        s = ts[-1] - x/ci[n]
        # calculate the x-positions along the characteristic for each value of ts
        xs = ci[n] * (ts[:, np.newaxis] - s[np.newaxis, :])
        fac = np.sqrt(1.0 + ci[n]**2)  # dxi = fac * dt
        psi += np.trapz(-ci[n] * bi[n]*X(xs, ts[:, np.newaxis])/f, fac*ts, axis=0)*F[:, n, np.newaxis]
    
    return psi

# Get the streamfunction
# 
# fig = plt.figure()
# 
# hu = 0.5 * (h[1:] + h[:-1])
# 
# Nsubsample = 3
# for n in range(Nt/Nsubsample):
#     psi = forced_modes(n*Nsubsample, Nnodes=100)
#     eta = (f / g) * np.cumsum(np.diff(psi, axis=0) / hu[:, np.newaxis], axis=0)
#     eta = np.vstack((eta, eta[-1, :]))
#     eta -= eta[-1, :]
#     fig.clf()
#     ax = fig.add_subplot(111)
#     pc = ax.pcolormesh(x/1e3, y/1e3, eta, vmin=-0.02, vmax=0.02, cmap=plt.cm.RdBu_r)
#     plt.colorbar(pc, ax=ax)
#     ax.contour(x/1e3, y/1e3, psi, range(-100000, 100000, 1000), colors='k')
#     
#     plt.savefig('frame_%04d.png' % n)
#     print ' ### Wrote frame %04d' % n
