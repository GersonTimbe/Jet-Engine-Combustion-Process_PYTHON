# Mathematical Modeling: Combustion Process Simulation in a Jet Engine Chamber (PYTHON) 
This project simulates a combustion process in a combustion chamber of a Jet-engine using conservation laws, chemical kinetics and the Crank-Nicolson numerical method.

**Goal:** Solve the 2D coupled system of partial diferential equations (mass, energy and species conservation equations) to describe temperature and chemical species variations during the injection and combustion process
.<p align ="justify"> **Keywords:** Computational Simulation, Python Programming, Numerical Methods, Finite Differences,  Mathematical Modeling, Chemical Kinetics, Heat Transfer, Combustion Chamber, Jet Engine.
</p>
   
**Language:** Python\
**Main Libraries:** NumPy, Scipy, Matplotlib, os\
**Applications:** Thermal behavior and Air-Fuel Mixture consumption in jet engines (TCC - Undergraduate Thesis-Applied Mathematics)

**Codes:**\
**Code 1:** *Ideal Chamber* \
 Generates the chamber used to do the simulations, that is an approximation of the real one.
This rectangular chamber was needed to make it possible to apply the numerical method used. </p>

**Code 2:** *Numerical Meshes* \
Generates the two meshes used to do the partition of the domain and discretize the equations.
The results of these codes can be seen without running the codes inside the folder "results" above.

**Code 3:** *Modeled Velocity Field* \
Generates a visualization of the created variable velocity field $\mathbf{v}=(v(x,y),u(x,y))$ used in simulation 3
.<p align ="justify"> **Code 4:** *Simulation 2* \
Generates the **Second Final Result of the Research**: How temperature and air-fuel mixture varies in the combustion chamber during  injection and during the combustion process. When the code is run it produces and store in the computer images of the simulation in a folder named "simulation2" and displays it on the screen to each 5th time step. These images were then used to produce the video of the simulation, named "simulation2(download to play).mp4" , (see the "result" folder above) , and need to be downloaded to be played as it is around 1.5mb, (heavy to github). In this code was used a constant velocity field
.<p align ="justify">**Code 5:** *Simulation 3* \
Generates the **Third Final Result of the Research**: How temperature and air-fuel mixture varies in the combustion chamber during  injection and during the combustion process, now using the *modeled velocity field*, that is more realistic to allow comparisons between the simulations and check the stability of the solutions. The video of the simulation also can be downloaded in the "result" folder above, named "simulation3(download to play).mp4" </p>
 

