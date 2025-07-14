# ===========================================================================
# Physical_parameters (Ajusted)
# ===========================================================================
##Fluid mechanics parameters (S.I) 
rho = 1                            # density of the mixture
Cp  = 1000                         # specific heat capacity 
kT   = 350                         # thermal conductivity (ajusted)
alpha = kT/(rho*Cp)                # thermal difusity

#constant velocity field
v_x=25                             #Maximum Avarage velocity (ajusted)                  
v_y = 0

#variable velocity field  v=(v_x,v_y) (m/s)
vx_max = 25                        #Maximum Avarage velocity (ajusted)
vy_max = vx_max/5
A = 1.7
B = (5*A*Ly/(np.pi*Lx))
v_x = vx_max * (1 - np.exp(-A * x/(1*Lx))) * np.sin(np.pi * y / (Ly))
v_y = B*vy_max* np.exp(-A * y / Lx) *( np.cos(np.pi * y / Ly))
magnitude = np.sqrt(v_x**2 + v_y**2)
max_V_magnitude =20 


# parameters for species transport:  Fuel C10H22(n-decane) and O_2
MW_fuel = 0.142                   # molar weight of the fuel (kg/mol)
MW_O2   = 0.032                   # molar weight of the oxygen (kg/mol)
S_ratio = 15.5*(MW_O2/MW_fuel)    #stoichiometric ratio
Phi     = 1                       # equivalent ratio   
Ae      = 10**8                   # pre-exponential factor for reaction rate
E_a     = 147000                  # activation energy (J/mol)
R       = 8.314                   # universal gas constant (J/mol·K)
Delta_H =-6345000                 # heat of combustion (J/mol)
D_fuelx = 0.35                    # diffusion coefficient for the fuel
D_fuely = 0.17
D_Ox    = 0.35                    # diffusion coefficient for the oxygen
D_Oy =0.17
K       = np.zeros((Ny, Nx))      #reaction constant (updated later) 

mPex = max_V_magnitude*Lx/D_fuel  #peclet number to mass diffusion
tPex = max_V_magnitude*Lx/alpha   #peclet number to thermal diffusion

# Parameters for Initial and Boundries conditions (Heat Transfer) (simplified notations) 
# to follow the the oriatation of python i vertical index we will call i=0 battom (Upper of the combustor) and i=Ny upper (bottom of the combustor) 
he = 10**4                       # convective heat transfer coefficient (Left)
hi = 10**6                       # (Bottom)
hd = 200                         # (Right)
hs = 10**6                       # (Upper)
Te = 800                         # temperature of the left side of the chamber
Tenv = 320                       # temperature of the external environment

#auxiliaries
Ce = 1+(he/kT)*delta_x
Ci = 1+(hi/kT)*delta_y
Cd = 1+(hd/kT)*delta_x
Cs = 1+(hs/kT)*delta_y
Ke = ((he/kT)*delta_x*Tenv)/Ce
Ks = ((hs/kT)*delta_y*Tenv)/Cs
Kd = ((hd/kT)*delta_x*Tenv)/Cd
Ki = ((hi/kT)*delta_y*Tenv)/Ci
mfr_fuel = 1                      # mass flow rate to the fuel (kg/m2*s - variable VF)
mfr_O = (1/Phi)*S_ratio*mfr_fuel  # mass flow rate to the oxidizer (kg/m2*s) 