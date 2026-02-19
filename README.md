# Jet-Engine-Combustion-Process-Simulation (PYTHON) 
<p align ="justify">
This project simulates a combustion process in a combustion chamber of a Jet-engine using conservation laws, chemical kinetics and the Crank-Nicolson numerical method.
</p>

**Goal:** <p align ="justify"> Solve the 2D coupled system of partial diferential equations (mass, energy and species conservation equations) to describe temperature and chemical species variations during the injection and combustion process. </p>

**Keywords:** <p align ="justify"> Computational Simulation, Python Programming, Numerical Methods, Finite Differences,  Mathematical Modeling, Chemical Kinetics, Heat Transfer, Combustion Chamber, Jet Engine.</p>
</p>
   
**Language:** Python  
**Main Libraries:** NumPy, Scipy, Matplotlib, os  
**Applications:** <p align ="justify"> Thermal behavior and Air-Fuel Mixture consumption in jet engines (TCC - Undergraduate Thesis-Applied Mathematics)</p>

**Codes:**\
**Code 1:** *Ideal Chamber* \
 Generates the chamber used to do the simulations, that is an approximation of the real one.
This retangular chamber was needed to make it possible to apply the numerical method used. </p>

**Code 2:** *Numerical Meshes* \
Generates the two meshes used to do the partition of the domain and discretize the equations.

**Code 3:** *Modeled Velocity Field* \
Generates a visualization of the created variable velocity field $**v**=(v(x,y),u(x,y))$ used in simulation 3.

**Code 4:** <p align ="justify"> *Simulation 2* \
Generates the **Second Final Result of the Research**: How temperature and air-fuel mixture varies in the combustion chamber during  injection and during the combustion process. The code produces and store images of the simulation to a folder named "simulation2" and displays it on the screen to each 5th time step. These images were then used to produce the video of the simulation. In this code was used a constant velocity fiel. </p>

**Code 5:** *Simulation 3* \
<p align ="justify"> Generates the **Third Final Result of the Research**: How temperature and air-fuel mixture varies in the combustion chamber during  injection and during the combustion process, now using the *modeled velocity field*, that is more realistic to allow comparisons.</p>
 
## Summary in Portuguese
Este projeto faz parte do meu Trabalho de Conclusão de Curso (Matemática Aplicada) e simula o processo de combustão em uma câmara de combustão ideal de um Motor a jacto usando leis de conservação, cinética química e o método numérico de Crank-Nicolson. 


Objectivo: Resolver o sistema acoplado de equações diferencias parciais (equação de conservação da massa, energia e de espécies químicas) para descrever a variação térmica e de espécies químicas na câmara durante a injeccao e durante o processo de combustão.

