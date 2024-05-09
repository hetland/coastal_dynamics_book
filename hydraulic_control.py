
import numpy as np
import matplotlib.pyplot as plt

N = 1001    # Number of bathymetric points
x = np.linspace(-3, 3, N)

# The bathymetry is a Gaussian bump, zero at the peak, negative elsewhere.
m = np.exp( -x**2 ) - 1.0

# The Bernoulli constant.
B = 1.5     # Exactly critical (Fr=1) at the peak (m=0)

# Solve the Bernoulli equation numerically for each bathymetric point
#   h**3 + (m - B) * h**2 + 0.5 = 0
r = []
for mi in m:
    r.append( np.roots((1, (mi-B), 0, 0.5)) )

r = np.asarray(r)   # This Nx3 array contains the three roots for each bathymetric point

fig = plt.figure()

h1 = r[:, 0]
h2 = r[:, 1]

Fr1 = h1**-1.5
u1 = Fr1 / h1

Fr2 = h2**-1.5
u2 = Fr2 / h2

# Plot the fluid thickness, h
ax_h = fig.add_subplot(311)
ax_h.plot(x, h1 + m, '-k')
ax_h.plot(x, h2 + m, '-k')
ax_h.fill(x, m, '0.5')
ax_h.plot(x, m, '-k', lw=1)

# Plot the flow speed, u / sqrt(g D)
ax_u = fig.add_subplot(312)
ax_u.plot(x, u1, '-k')
ax_u.plot(x, u2, '-k')

# Plot the Froude number  Fr = h**-3/2
ax_Fr = fig.add_subplot(313)
ax_Fr.plot(x, Fr1, '-k')
ax_Fr.plot(x, Fr2, '-k')



# The Bernoulli constant.
B = 2.0     # Super- or subcritical at the peak

# Solve the Bernoulli equation numerically for each bathymetric point
#   h**3 + (m - B) * h**2 + 0.5 = 0
r = []
for mi in m:
    r.append( np.roots((1, (mi-B), 0, 0.5)) )

r = np.asarray(r)   # This Nx3 array contains the three roots for each bathymetric point

# The first two roots are the real roots, and correspond to the fluid thickness.
h1 = r[:, 0]
h2 = r[:, 1]

Fr1 = h1**-1.5
u1 = Fr1 / h1

Fr2 = h2**-1.5
u2 = Fr2 / h2

# Plot the fluid thickness, h
ax_h.plot(x, h1 + m, '-k', lw=0.5)
ax_h.plot(x, h2 + m, '-k', lw=0.5)

# Plot the flow speed, u / sqrt(g D)
ax_u.plot(x, u1, '-k', lw=0.5)
ax_u.plot(x, u2, '-k', lw=0.5)

# Plot the Froude number  Fr = h**-3/2
ax_Fr.plot(x, Fr1, '-k', lw=0.5)
ax_Fr.plot(x, Fr2, '-k', lw=0.5)


ax_h.set_xticks([])
ax_u.set_xticks([])
ax_Fr.set_xticks([])

ax_h.set_yticks([0, 1, 2])
ax_u.set_yticks([0, 2, 4, 6, 8])
ax_Fr.set_yticks([0, 1, 2, 3])


plt.show()

