import numpy as np
import matplotlib.pyplot as plt

# Combustion chamber dimensions
Lx, Ly = 0.6, 0.3
nx, ny = 70, 35
x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
X, Y = np.meshgrid(x, y)

#Parameters:
vx_max = 100         # Maximum average velocity --ajusted (m/s)
vy_max = vx_max/5
A = 1.7
B = (5*A*Ly/(np.pi*Lx))

# Velocity distribution considering injection and mass flow
vx = vx_max * (1 - np.exp(-A * X/(1*Lx))) * np.sin(np.pi * Y / (Ly))
vy = B*vy_max* np.exp(-A * X / Lx) *( np.cos(np.pi * Y / Ly))

# velocity magnitude
magnitude = (1/4)* np.sqrt(vx**2 + vy**2)

# Plotting the veloxity field
fig, ax = plt.subplots(figsize=(10, 5))
plt.quiver(X, Y, vx, vy, magnitude, scale=1000, cmap='jet')
plt.colorbar(label='Velocity Magnitude (m/s)')
plt.title('Prescribed Velocity field (v_max = 20 m/s)')
plt.xlabel('x')
plt.ylabel('y')

ax.set_xlim([0,Lx])
ax.set_ylim([0,Ly])
ax.set_aspect(1)

plt.show()

