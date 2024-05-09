
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

################################################################################
# Define velocity profile

N  = 20
xi = np.linspace(-0.5, 0.5, N)
vhat = -np.cos(np.pi * xi)**2

# Uncomment these three lines to see what happens when the jet is not uniform.
# vhat[14] = -3
# vhat[15] = -5
# vhat[16] = -3

################################################################################
# Function definitions

def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr


def split(alpha):
    'return xi position of splitting streamfunction, percent transport and momentum'
    ip = np.trapz(vhat**2, xi)
    ipi = np.empty_like(xi)
    for n in range(1, len(xi)):
        ipi[n-1] = np.trapz(vhat[:n]**2, xi[:n])
    
    ipi[n] = np.trapz(vhat**2, xi)
    ipi -= 0.5 * (1 - np.sin(alpha)) * ip
    
    idx = np.argwhere(ipi >= 0)
    
    if idx is None:
        return xi[-1], 0.0, 0.0 # all transport goes right
    
    i = idx.min()
    if i == 0: 
        return xi[0], 1.0, 1.0 # all transport goes left
    
    yu = ipi[i]
    yl = ipi[i-1]
    xu = xi[i]
    xl = xi[i-1]
    xo = xl - (xu-xl)*yl/(yu-yl)
    
    vu = vhat[i]
    vl = vhat[i-1]
    vo = vl + (xo-xl)*(vu-vl)/(xu-xl)
    ivl = np.trapz(np.hstack((vhat[:i], vo)), np.hstack((xi[:i], xo))) / np.trapz(vhat, xi)       
    ipl = np.trapz(np.hstack((vhat[:i], vo))**2, np.hstack((xi[:i], xo))) / ip
    
    return xo, 100.0*ivl, 100.0*ipl

################################################################################
# Initial plotting

fig = plt.figure()
ax = fig.add_axes([0, 0.1, 1, 0.9])

# plot domain parameters
a = 2
b = 2
vscale = 10*np.abs(vhat).max()

ax.plot([-a, -a, a, a], [0, b, b, 0], '-k', alpha=0.25, lw=3.0)
ax.plot([-a, a], [0, 0], '-k', lw=3.0)


alpha = 0.0
x, y = rot2d(xi, 0.0, alpha)
u, v = rot2d(0, vhat, alpha)
pl_off = ax.plot(x, y+b, '-k')
q_off = ax.quiver(x, y+b, u, v, scale=vscale, width=0.003)

xo, lower_trans_percent, lower_momentum_percent = split(alpha)
v_right = np.ma.masked_where(xi < xo, vhat)
v_left = np.ma.masked_where(xi > xo, vhat)
q_right = ax.quiver(a*np.ones_like(xi), xi-xo, -v_right, 0, scale=vscale, width=0.003)
q_left = ax.quiver(-a*np.ones_like(xi), xo-xi, v_left, 0, scale=vscale, width=0.003)

txt_left = ax.text(-a+a/10.0, a/10.0, 
    '%4.1f%% Transport\n%4.1f%% Momentum' % (lower_trans_percent, lower_momentum_percent), 
    horizontalalignment='left', verticalalignment='bottom', fontsize=10)

txt_right = ax.text(a-a/10.0, a/10.0, 
    '%4.1f%% Transport\n%4.1f%% Momentum' % (100-lower_trans_percent, 100-lower_momentum_percent), 
    horizontalalignment='right', verticalalignment='bottom', fontsize=10)

txt_angle = ax.text(0, b/2.0, 
    r'$\alpha$ = %4.1f$^{\circ}$' % (alpha*180/np.pi,), 
    horizontalalignment='center', verticalalignment='bottom', fontsize=12)


ax.set_aspect(1.0)
ax.set_xlim(-a-1, a+1)
ax.set_ylim(-1, b+1)
ax.set_axis_off()

################################################################################
# Update figure based on user interaction

ax_alpha = fig.add_axes([0.13, 0.07, 0.77, 0.03])
salpha = Slider(ax_alpha, r'$\alpha$', -np.pi/2.0, np.pi/2.0, valinit=alpha)

def update(alpha):
    x, y = rot2d(xi, 0.0, alpha)
    u, v = rot2d(0, vhat, alpha)
    
    pl_off[0].set_data(x, y+b)
    q_off.set_offsets(zip(x, y+b))
    q_off.set_UVC(u, v)
    
    xo, lower_trans_percent, lower_momentum_percent = split(alpha)
    v_right = np.ma.masked_where(xi < xo, vhat)
    v_left = np.ma.masked_where(xi > xo, vhat)
    
    q_right.set_offsets(zip(a*np.ones_like(xi), xi-xo))
    q_left.set_offsets(zip(-a*np.ones_like(xi), xo-xi))
    q_right.set_UVC(-v_right, 0)
    q_left.set_UVC(v_left, 0)
    
    txt_left.set_text('%4.1f%% Transport\n%4.1f%% Momentum' % 
                      (lower_trans_percent, lower_momentum_percent))
    txt_right.set_text('%4.1f%% Transport\n%4.1f%% Momentum' % 
                       (100-lower_trans_percent, 100-lower_momentum_percent))
    txt_angle.set_text(r'$\alpha$ = %4.1f$^{\circ}$' % (alpha*180/np.pi,))

    
    plt.draw()

salpha.on_changed(update)

plt.show()
