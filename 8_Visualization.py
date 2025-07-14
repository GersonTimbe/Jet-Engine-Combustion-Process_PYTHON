# ===========================================================================
#Visualization
# ===========================================================================
#Building a personalised colormap the mixture faction
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors 
if not combustion_started:             # before combustion   
    colorsf = [
        (0.00, '#440154'),  
        (0.10, '#472a7a'), 
        (0.20, '#3b518b'), 
        (0.30, '#2c6e8f'),  
        (0.40, '#25858e'),  
        (0.50, '#1f9e89'),
        (0.60, '#35b779'), 
        (0.70, '#6ece58'), 
        (0.80, '#b5de2b'), 
        (0.90, '#fde725'), 
        (1.00, '#ffffe5')  
    ]
    custom_cmapf = mcolors.LinearSegmentedColormap.from_list("custom_fuel", colorsf)

    colorst = [
        (0.0, "#000000"), 
        (0.1, "#200000"), 
        (0.2, "#400000"),  
        (0.3, "#600000"),  
        (0.4, "#800000"),
        (0.5, "#990000"),
        (0.6, "#a00000"),
        (0.7, "#b00000"),
        (0.8, "#c00000"),
        (0.9, "#d00000"),
        (1.0, "#e30808")
    ]

    custom_cmapt = mcolors.LinearSegmentedColormap.from_list("custom_temperature", colorst)
    T_min, T_max = 300, 800 
    Ymist_min, Ymist_max = 0.0, 1 

    #Image resolution definitions
    dpi = 100
    width_inches =1366/dpi
    height_inches =673/dpi
    plt.figure(figsize=(width_inches,height_inches),dpi=dpi)
    plt.clf()

    plt.subplot(1, 2, 1)                 #TEMPERATURE
    img1 = plt.imshow(T, cmap=custom_cmapt, extent=[0, Lx, 0, Ly], vmin=T_min, vmax=T_max)
    cbarT = plt.colorbar(img1)
    cbarT.label=''
    ticks = cbarT.get_ticks()                          
    labels = [f"{tick:.0f}" for tick in ticks]         
    labels[-1] += " K"                   
    cbarT.set_ticklabels(labels)        
    plt.xticks([0.0, Lx])                 
    plt.yticks([0.0, Ly])
    plt.title('Distribuição de Temperatura')

    plt.subplot(1, 2, 2)                #MIXTURE FRACTION
    plt.imshow(Ymist, cmap=custom_cmapf, extent=[0, Lx, 0, Ly], vmin=Ymist_min, vmax=Ymist_max)
    plt.colorbar(label='')
    plt.xticks([0.0, Lx])              
    plt.yticks([0.0, Ly])
    plt.title('Distribuição da Mist. Ar-Combustível')

    # Exibir passo de tempo na tela
    plt.suptitle(f"Passo de tempo: {t}")

else:                               #During and after ignition
    colorsf = [
        (0.00, '#440154'), 
        (0.02, '#472a7a'),  
        (0.04, '#3b518b'), 
        (0.06, '#2c6e8f'),  
        (0.08, '#25858e'), 
        (0.10, '#1f9e89'),  
        (0.20, '#35b779'),  
        (0.40, '#6ece58'),  
        (0.60, '#b5de2b'),  
        (0.80, '#fde725'),  
        (1.00, '#ffffe5')  
    ]
    custom_cmapf = mcolors.LinearSegmentedColormap.from_list("custom_fuel",
    colorsf)

    colorst = [
        (0.0, "#000000"),  
        (0.1, "#540000"),  
        (0.2, "#9f0000"),  
        (0.3, "#e90000"),  
        (0.4, "#ff3500"), 
        (0.5, "#ff7f00"),  
        (0.6, "#ffc900"), 
        (0.7, "#ffff1f"),  
        (0.8, "#ffff8f"), 
        (0.9, "#ffffdd"),  
        (1.0, "#ffffff")   
    ]
    custom_cmapt = mcolors.LinearSegmentedColormap.from_list("custom_fuel", colorst)

    T_min, T_max = 300, 3000         
    Ymist_min, Ymist_max = 0.0, 0.8 

    #Image resolution definitions
    dpi = 100
    width_inches =1366/dpi
    height_inches =673/dpi
    plt.figure(figsize=(width_inches,height_inches),dpi=dpi)
    plt.clf()

    plt.subplot(1, 2, 1)                 #TEMPERATURE
    img1 = plt.imshow(T, cmap=custom_cmapt, extent=[0, Lx, 0, Ly], vmin=T_min, vmax=T_max)
    cbarT = plt.colorbar(img1)
    cbarT.label=''
    ticks = cbarT.get_ticks()                         
    labels = [f"{tick:.0f}" for tick in ticks]        
    labels[-1] += " K"                  
    cbarT.set_ticklabels(labels)        
    plt.xticks([0.0, Lx])               
    plt.yticks([0.0, Ly])
    plt.title('Distribuição de Temperatura')

    plt.subplot(1, 2, 2)                #MIXTURE FRACTION
    plt.imshow(Ymist, cmap=custom_cmapf, extent=[0, Lx, 0, Ly], vmin=Ymist_min, vmax=Ymist_max)
    plt.colorbar(label='')
    plt.xticks([0.0, Lx])              
    plt.yticks([0.0, Ly])
    plt.title('Distribuição da Mist. Ar-Combustível')

    # Show time step
    plt.suptitle(f"Passo de tempo: {step}")
import os
if t%5==0 or t in [324,325, 326, 327, 328]:
    os.makedirs(f"imagens_da_simulacaoVLCTF",exist_ok=True)
    plt.savefig(f"imagens_da_simulacaoVLCTF/simulation_step_{step:04d}.png",dpi=150,bbox_inches='tight')
plt.close()