import numpy as np
import matplotlib.pyplot as plt

dados = np.loadtxt("dados_simulacao2nd2.csv",delimiter=",",skiprows=1)
time_steps = dados[:,0]
T_mean_list = dados[:,1]
Yf_mean_list = dados[:,2]
YO_mean_list = dados[:,3]
Ymist_mean_ignition_list = dados[:,4]



time_steps_to_plotT = [t for t in range(0,1201) if t%40==0] 
time_steps_to_plotT = np.array(time_steps_to_plotT)

time_steps_to_plotF = [t for t in range(0,601) if t%20==0]
time_steps_to_plotF = np.array(time_steps_to_plotF )

time_steps_to_plotO = [t for t in range(0,601) if t%20==0 ]
time_steps_to_plotO = np.array(time_steps_to_plotO)

# Filtrar apenas os valores correspondentes a esses tempos
T_mean_selected = [T_mean_list[t] for t in time_steps_to_plotT]
Yf_selected = [Yf_mean_list[t] for t in time_steps_to_plotF]
YO_selected = [YO_mean_list[t] for t in time_steps_to_plotO]
Ymist_ignition_selected = [Ymist_mean_ignition_list[t] for t in time_steps_to_plotF]

# Criando os gráficos no final da simulação
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Gráfico da temperatura média
markers_pointsM = [100, 180,200,220,500]
indicesM = [np.argmin(np.abs(time_steps_to_plotF -t )) for t in markers_pointsM]

markers_pointsF = [100,200,220,240,500]
indicesF = [np.argmin(np.abs(time_steps_to_plotF -t )) for t in markers_pointsF]

markers_pointsO = [ 100,200,220,240,500]
indicesO = [np.argmin(np.abs(time_steps_to_plotF -t )) for t in markers_pointsO]
axes[0].plot(time_steps_to_plotT, T_mean_selected, 'r-', label="Temperatura Média", linewidth=2)
axes[0].set_xlabel("Tempo (passos)")
axes[0].set_ylabel("Temperatura (K)")
axes[0].set_title("Evolução da Temperatura Média")
axes[0].legend()
axes[0].grid(True)
axes[0].set_xlim(left=0)

# Gráfico da variação do combustível e oxidante
axes[1].plot(time_steps_to_plotF, Yf_selected, 'm-o',markevery = indicesF, label="F. média do combustivel", linewidth=2)
axes[1].plot(time_steps_to_plotO, YO_selected, 'b-s', markevery = indicesO,label="F. média do oxidante", linewidth=2)
axes[1].plot(time_steps_to_plotF, Ymist_ignition_selected, 'g-^' ,markevery = indicesM ,label="F. mistura na região de ignição", linewidth=2) #'g-^'
axes[1].axvline(x=178, color='red', linestyle='--',label="ignição", linewidth=1.5)
axes[1].set_xlabel("Tempo (passos)")
axes[1].set_ylabel("Fração mássica")
axes[1].set_title("Variação da média do combustivel e oxidante")
axes[1].legend(loc='upper right')
axes[1].grid(True)
axes[1].set_xlim(left=0)

plt.tight_layout()
plt.show()