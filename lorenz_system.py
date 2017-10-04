import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

def rhs_lorenz(old_solution, t, sigma, rho, beta):
    new_solution = np.zeros(3)
    x, y, z = old_solution[0], old_solution[1], old_solution[2]
    new_solution[0] = sigma*(y - x)
    new_solution[1] = x*(rho - z) - y
    new_solution[2] = x*y - beta*z 

    return new_solution
    
# lorenz values
sigma_l, beta_l, rho_l = 10, 8/3, 28

# initial condition
sol_init = np.array([1, 1, 1])

# time vector
dt = 0.01
t = np.arange(0, 150, dt)

# use an explicit 4-5 embedded RK pair
sol = odeint(rhs_lorenz, sol_init, t, args=(sigma_l, rho_l, beta_l))
sol2 = odeint(rhs_lorenz, sol_init+.001, t, args=(sigma_l, rho_l, beta_l))

plt.close('all')
fig = plt.figure()
ax = fig.gca(projection="3d")
ax.set_xlim((np.min(sol[:,0]), np.max(sol[:,0])))
ax.set_ylim((np.min(sol[:,1]), np.max(sol[:,1])))
ax.set_zlim((np.min(sol[:,2]), np.max(sol[:,2])))


for i, ti in enumerate(t):
    # ptlot the initial condition
    ax.cla()
    ax.set_xlim((np.min(sol[:,0]), np.max(sol[:,0])))
    ax.set_ylim((np.min(sol[:,1]), np.max(sol[:,1])))
    ax.set_zlim((np.min(sol[:,2]), np.max(sol[:,2])))
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_zlabel(r'$z$')
    ax.scatter(sol[0,0], sol[0,1], sol[0,2], color='C0')
    ax.scatter(sol2[0,0], sol2[0,1], sol2[0,2], color='C1')
    # plot the trajectory
    ax.plot(sol[:i,0], sol[:i,1], sol[:i,2], color='C0', linewidth=0.5)
    ax.plot(sol2[:i,0], sol2[:i,1], sol2[:i,2], color='C1', linewidth=0.5)
    # plot the latest point
    ax.scatter(sol[i,0], sol[i,1], sol[i,2], color='C0')
    ax.scatter(sol2[i,0], sol2[i,1], sol2[i,2], color='C1')
    plt.pause(1e-9)
