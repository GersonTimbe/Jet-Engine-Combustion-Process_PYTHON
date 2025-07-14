import matplotlib.pyplot as plt
import numpy as np

# Chamber dimensions
Lx = 6  # Lenght
Ly = 4.8  # Height
fig, ax = plt.subplots(figsize=(6, 4))

# vertex of the chamber
x = np.array([0, 0.3, 0.3, 4.7, 4.7, 5])  
y1 = np.array([(Ly-2)/3, (Ly-2)/3, 0, 0, (Ly-2)/3, (Ly-2)/3])  # Parte inferior
y2 = np.array([2*(Ly-2)/3, 2*(Ly-2)/3, Ly-2, Ly-2, 2*(Ly-2)/3, 2*(Ly-2)/3])  # Parte superior

ax.plot(x, y1, 'k', linewidth=3)
ax.plot(x, y2, 'k', linewidth=3)

# random blue dots
x_min, x_max = 0.15,1.7  
y_min, y_max = 1.10, 1.9  

num_pontos = 30  # amount of dots
x_pontos = np.random.uniform(x_min, x_max, num_pontos)
y_pontos = np.random.uniform(y_min, y_max, num_pontos)

# dots plotting
ax.scatter(x_pontos, y_pontos, s=30, marker=".", color="blue", label="Air + Fuel")

# mixture entrance
ax.arrow(-1.8, (Ly-2)/2, 1.5, 0, head_width=0.2, head_length=0.2, fc='black', ec='black')
ax.text(-1.5, (Ly-2)/2 + 0.4, "Entry", fontsize=8, color='black')
ax.text(-1.8, (Ly-2)/2 - 0.4, "Air + Fuel", fontsize=8, color='black')

# Ponto de ignição 
ax.scatter(1.3, (Ly-2)/2, s=30, facecolors='red', edgecolors='red', linewidths=2, label="Ignition Point" )

#Lx
ax.annotate("", xy=(4.7+0.075, -0.4), xytext=(0.3-0.075, -0.4),
            arrowprops=dict(arrowstyle="<->", linewidth=1.5))
ax.text((Lx-1)/2, -0.75, r"$L_x$", fontsize=12, ha="center")

#Ly
ax.annotate("", xy=(5.4, (Ly-1)-0.93), xytext=(5.4 , -0.07),
            arrowprops=dict(arrowstyle="<->", linewidth=1.5))
ax.text(Lx-0.5, (Ly-2)/2, r"$L_y$", fontsize=12, va="center")


# Ajustando limites do gráfico
ax.set_xlim(-0.2, 1.2)
ax.set_ylim(-0.2, 1.2)

# Ajustes finais
ax.set_xlim(-2, 6)
ax.set_ylim(-2, 4)
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
# Exibir o gráfico
plt.legend(loc="upper left")
plt.show()