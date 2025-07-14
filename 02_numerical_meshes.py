# ===========================================================================
# Numerical_meshes
# ===========================================================================
import numpy as np
import matplotlib.pyplot as plt
#FIRST:

# Domain caracteristics
Lx, Ly = 0.6, 0.3
delta_x = 0.6/70                # spatial step size (x - direction)
delta_y = 0.3/35                # spatial step size (y - direction)
Nx = int(round(Lx/delta_x)) + 1 # Number of points on the x axis. (Nx=71)
Ny = int(round(Ly/delta_y)) + 1 # Number of points on the y axis. (Ny=36)
Ni = (Nx-2)*(Ny-2)              # Number of internal points (Number of Unknowns: 2006 )
Nt = 1405                       # Number of time steps.

# Funcion to generate the first mesh
def plot_mesh(Nx, Ny, ax, title):
    x = np.linspace(0, Lx, Nx)
    y = np.linspace(0, Ly, Ny)
    X, Y = np.meshgrid(x, y)
    
    ax.plot(X, Y, 'k-', lw=0.5)   
    ax.plot(X.T, Y.T, 'k-', lw=0.5)  
    ax.set_title(title)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_xticks([0, Lx])
    ax.set_yticks([0, Ly])
    ax.set_aspect('equal')
    ax.scatter(X.flatten(), Y.flatten(), color='black', s=2)

# figure
fig, ax = plt.subplots( figsize=(10, 5))
plot_mesh(70, 35, ax, "Numerical Mesh 70 × 35")

plt.tight_layout()
plt.show()



#SECOND
# Domain caracteristics
Lx, Ly = 0.6, 0.3
delta_x = 0.6/100                 # spatial step size (x - direction)
delta_y = 0.3/50                  # spatial step size (y - direction)
Nx = int(round(Lx/delta_x)) + 1   # Number of points on the x axis. (Nx=101)
Ny = int(round(Ly/delta_y)) + 1   # Number of points on the y axis. (Ny=51)
Ni = (Nx-2)*(Ny-2)                # Number of internal points (Number of Unknowns: 4851 )
Nt = 1210                         # Number of time steps.

# Funcion to generate the second mesh
def plot_mesh(Nx, Ny, ax, title):
    x = np.linspace(0, Lx, Nx)
    y = np.linspace(0, Ly, Ny)
    X, Y = np.meshgrid(x, y)
    
    ax.plot(X, Y, 'k-', lw=0.5)  
    ax.plot(X.T, Y.T, 'k-', lw=0.5) 
    ax.set_title(title)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_xticks([0, Lx])
    ax.set_yticks([0, Ly])
    ax.set_aspect('equal')
    ax.scatter(X.flatten(), Y.flatten(), color='black', s=1)

# figure
fig, ax = plt.subplots( figsize=(10, 5))
plot_mesh(100, 50, ax, "Numerical Mesh 100 × 50")

plt.tight_layout()
plt.show()
