import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("dados_simulacaoVLCTF.csv",delimiter=",",skiprows=1)
time_steps = dados[:,0]
T_mean_list = dados[:,1]
Yf_mean_list = dados[:,2]
YO_mean_list = dados[:,3]
Ymist_mean_ignition_list = data[:,4]

# specific timing to plot
time_steps_to_plotT = [t for t in range(0,1399) if t%10==0] 
time_steps_to_plotT = np.array(time_steps_to_plotT)

time_steps_to_plotF = [t for t in range(0,701) if t%5==0]
time_steps_to_plotF = np.array(time_steps_to_plotF )

time_steps_to_plotO = [t for t in range(0,701) if t%5==0 ]
time_steps_to_plotO = np.array(time_steps_to_plotO)

# Filter only the values of the coresponding time
T_mean_selected = [T_mean_list[t] for t in time_steps_to_plotT]
Yf_selected = [Yf_mean_list[t] for t in time_steps_to_plotF]
YO_selected = [YO_mean_list[t] for t in time_steps_to_plotO]
Ymist_ignition_selected = [Ymist_mean_ignition_list[t] for t in time_steps_to_plotF]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Temperature avarage
markers_pointsM = [200,324,330,365,370]
indicesM = [np.argmin(np.abs(time_steps_to_plotF -t )) for t in markers_pointsM]

markers_pointsF = [200,325]
indicesF = [np.argmin(np.abs(time_steps_to_plotF -t )) for t in markers_pointsF]

markers_pointsO = [200,325]
indicesO = [np.argmin(np.abs(time_steps_to_plotF -t )) for t in markers_pointsO]

axes[0].plot(time_steps_to_plotT, T_mean_selected, 'r-', label="Average Temperature", linewidth=2)
axes[0].set_xlabel("Teme (steps)")
axes[0].set_ylabel("Temperature (K)")
axes[0].set_title("Average Temperature evolution")
axes[0].legend()
axes[0].grid(True)
axes[0].set_xlim(left=0)

# Gráfico da variação do combustível e oxidante
axes[1].plot(time_steps_to_plotF, Yf_selected, 'm-o',markevery = indicesF, label="F. média do combustivel", linewidth=2,markersize=4)
axes[1].plot(time_steps_to_plotO, YO_selected, 'b-s', markevery = indicesO,label="F. média do oxidante", linewidth=2,markersize=4)
axes[1].plot(time_steps_to_plotF, Ymist_ignition_selected, 'g-^', markevery = indicesM,label="F. mistura na região de ignição", linewidth=2)
axes[1].axvline(x=326, color='red', linestyle='--',label="ignição", linewidth=1.5)
axes[1].set_xlabel("Tempo (passos)")
axes[1].set_ylabel("Fração mássica")
axes[1].set_title("Variação da média do combustivel e oxidante")
axes[1].legend(loc='upper right')
axes[1].grid(True)
axes[1].set_xlim(left=0)

plt.tight_layout()
plt.show()