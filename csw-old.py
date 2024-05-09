
import numpy as np
import matplotlib.pyplot as plt

alpha = 3e-3
Ho = 1e-6   # maintain a small, positive depth at the coast 
            # to eliminate a singularity in the definition of b_i
f = 1e-4
L = 100e3
N = 999

def h(y):
    return Ho + alpha * y

def h_y(y):
    return alpha * np.ones_like(y)

# Create the cross-shore coordinate
y = np.linspace(0, L, N)
y = y[1:-1]           # Only consider the interior of the domain
dy = y[2] - y[1]    # Define dy based on L and N

# Geometric factors based on the bathymetry and Coriolis parameter
fac1 = h(y) / (f * h_y(y) * dy**2)
fac2 = h_y(y) / (2 * h(y) * dy)

u = -fac1 + fac2    # upper (i+1) diagonal
d = 2*fac1          # center (i) diagonal
l = -fac1 - fac2    # lower (i-1) diagonal

# construct tri-diagonal matrix representing a discrete representation
# of the boundary value problem given by equation 2.12.14
A = np.diag(u[:-1], k=1) + np.diag(d, k=0) + np.diag(l[1:], k=-1)
# Apply the boundary condition that dF/dy=0 at y=L.
A[-1, -1] += -fac1[-1] + fac2[-1]
# (Note, nothing else needs to be done for the F=0 at y=0 boundary condition.)

# Calculate the eigenvalues and eigenvectors, and sort by eigenvalue.
lam, v = np.linalg.eig(A)
idx = np.argsort(lam)
ci = 1.0 / lam[idx]   # The phase speed is the inverse eigenvalue.
# Add in the two boundary values to the eigenvectors.
v = np.vstack((np.zeros((1, v.shape[1])), v[:, idx], v[-1, idx]))
# Add boundary values to the cross-shore coordinate, to match.
y = np.hstack((0, y, y[-1]+dy))


bi = np.trapz(v[:].T * h_y(y) / h(y)**2, x=y)


def display_modes(N):
    for mode in range(N):
        plt.plot(y, bi[mode] * v[:, mode])
        print 'mode ', mode, ', phase speed = ', ci[mode]
    plt.show()

display_modes(5)


