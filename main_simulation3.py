"""
@author: Gerson Joao Timbe

Script to simulate the combustion Process in a combustion chamber of a jet engine (assumptions a made) 
The main goal is to ilustrate the temperature and air-fuel injection and consumption variations

Simulation3
"""

import numpy as np
import scipy
import matplotlib
import os
import gc


#I. COEFFICIENTS OF THE MESH AND TIME STEP
Lx, Ly = 0.6, 0.3       
delta_x = 0.6/70      
delta_y = 0.3/35      
Nx = int(round(Lx/delta_x)) + 1 
Ny = int(round(Ly/delta_y)) + 1 
Ni = (Nx-2)*(Ny-2)             

Nt = 1405           
Dt = delta_x**2/(4)                 
x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
x, y = np.meshgrid(x, y)


#II. PARAMETERS
#Fluid Mechanics PARAMETERS ajusted to avoid instabilities
#more accurate parameters demand a longer simulation time (more than 20h in a normal computer)
rho = 1            
Cp  = 1000         
kT   = 350         
alpha = kT/(rho*Cp)                #diffusity

#Variable Velocity field  v=(v_x,v_y)
vx_max = 25       
vy_max = vx_max/5
A = 1.7
B = (5*A*Ly/(np.pi*Lx))
v_x = vx_max * (1 - np.exp(-A * x/(1*Lx))) * np.sin(np.pi * y / (Ly))
v_y = B*vy_max* np.exp(-A * y / Lx) *( np.cos(np.pi * y / Ly))
magnitude = np.sqrt(v_x**2 + v_y**2)
max_V_magnitude =20 


# Parameters for species transport:  Fuel C10H22(n-decane) and O_2
MW_fuel = 0.142                  # fuel molar weight
MW_O2   = 0.032                  # fuel molar weight       
S_ratio = 15.5*(MW_O2/MW_fuel)   # equivalence ratio
Phi     = 1              
Ae      = 10**8                  #|parameters of th erenhius model     
E_a     = 147000                 #|  
R       = 8.314                  #|  
Delta_H =-6345000                # combustion entalphy 
D_fuelx = 0.35                   #|Difusion constants 
D_Ox    = 0.35                   #|      
D_fuely = 0.17                   #|
D_Oy =0.17                       #|
K       = np.zeros((Ny, Nx))    


# Parameters for Initial and Boundries conditions (Heat Transfer) (simplified notations)  
he = 10**4  
hi = 10**6 
hd = 200   
hs = 10**6  
Te = 800  
Tenv = 320                       #Temperature of the environment

#auxiliaries
Ce = 1+(he/kT)*delta_x
Ci = 1+(hi/kT)*delta_y
Cd = 1+(hd/kT)*delta_x
Cs = 1+(hs/kT)*delta_y
Ke = ((he/kT)*delta_x*Tenv)/Ce
Ks = ((hs/kT)*delta_y*Tenv)/Cs
Kd = ((hd/kT)*delta_x*Tenv)/Cd
Ki = ((hi/kT)*delta_y*Tenv)/Ci
mfr_fuel = 1                   
mfr_O = (1/Phi)*S_ratio*mfr_fuel 


#1. CONSTRUCT THE MATRICES OF THE SYSTEM OF EQUATIONS.
A = np.zeros((Ni, Ni))
M = np.zeros((Ni, Ni))

# Function to fill matrices A and M for the temperature equation
def filling_matrixes_temperature(v_x,v_y):
    for k in range(Ni):
        i = 1 + k // (Nx-2)
        j = 1 + k % (Nx-2)
        # Diffusion and convection coefficients(energy equation)
        Dx = (Dt * alpha) / (2 * (delta_x)**2)         
        Dy = (Dt * 0.4*alpha) / (2 * (delta_y)**2)     
        Cxij = -(v_x[i,j] * Dt) / (2 * delta_x)      
        Cyij = -(v_y[i,j] * Dt) / (4 * delta_y)       

        # Coefficients for the discretized equation(energy equation)
        Sbx = Dx - Cxij   
        Smx = Dx + 0*Cxij  
        Sby = Dy - Cyij    
        Smy = Dy + Cyij   

        a = 1 + 2 * Dx + 2 * Dy - Cxij   
        b = 1 - 2 * Dx - 2 * Dy + Cxij   
       
        if i == 1:  
            if j == 1:
                A[k, k] =  a - Sby / Ci - Sbx/Ce
                M[k, k] =  b + Sby / Ci + Sbx/Ce
            elif j == Nx - 2:  
                A[k, k] =  a - Sby / Ci - Smx / Cd
                M[k, k] =  b + Sby / Ci + Smx / Cd
            else:  
                A[k, k] =  a - Sby / Ci
                M[k, k] =  b + Sby / Ci

        elif i == Ny - 2:  
            if j == 1:  
                A[k, k] = a - Smy / Cs - Sbx/Ce
                M[k, k] = b + Smy / Cs + Sbx/Ce
            elif j == Nx - 2:
                A[k, k] = a - Smy / Cs - Smx / Cd
                M[k, k] = b + Smy / Cs + Smx / Cd
            else:  
                A[k, k] = a - Smy / Cs
                M[k, k] = b + Smy / Cs

        elif j == 1 and i!=1 and i!=Ny-2:  
            i_min = max(1, int(round(Ny / 3))-1)          
            i_max = min(Ny - 2, int(round(2 * Ny / 3))+1)  
            if i_min < i < i_max:
                A[k, k] = a
                M[k, k] = b
            else:
                A[k, k] = a - Sbx/Ce
                M[k, k] = b + Sbx/Ce
        elif j == Nx - 2 and i!=1 and i!=Ny-2:  
             A[k, k] = a - Smx / Cd
             M[k, k] = b + Smx / Cd

        else:  
             A[k, k] = a
             M[k, k] = b

        #Non-Diagonal Terms    
        if i != 1:     
            A[k, (i-2)*(Nx-2) + (j-1)] = -Sby
            M[k, (i-2)*(Nx-2) + (j-1)] =  Sby

        if i != Ny-2: 
            A[k, i*(Nx-2) + (j-1)] = -Smy
            M[k, i*(Nx-2) + (j-1)] =  Smy

        if j != 1:    
            A[k, k-1] = -Sbx
            M[k, k-1] =  Sbx

        if j != Nx-2:  
            A[k, k+1] = -Smx
            M[k, k+1] =  Smx

# Initializing matrices for fuel and oxydizer transport (A_f, A_O, M_f and M_O)

# Function to fill matrices for f transport (A_f and M_f)
def filling_matrixes_f(A_f,M_f,T, YO,v_x,v_y):
    for k in range(Ni):
        i = 1 + k // (Nx-2)
        j = 1 + k % (Nx-2)
        #Arrenhius Model
        K[i,j]=Ae*np.exp(-E_a/(R*T[i,j]))

        #Coefficients for the discretized equation(fuel equation)
        Dfx = D_fuelx* Dt / (2 * (delta_x)**2)
        Dfy = D_fuely * Dt / (2 * (delta_y)**2)
        Cfxij = -(v_x[i,j] * Dt) / (2 * delta_x)  
        Cfyij = -(v_y[i,j] * Dt) / (4 * delta_y)      

        Sbxf = Dfx - Cfxij 
        Smxf = Dfx + 0*Cfxij 
        Sbyf = Dfy - Cfyij  
        Smyf = Dfy + Cfyij
 
        #Coefficients for fuel equation(2)
        r_ij_fuel = ((Dt )/(2*rho))*K[i,j]*(rho**2/(MW_O2))*YO[i,j]  
        a_fij = 1 + 2 * Dfx + 2 * Dfy - Cfxij +r_ij_fuel
        b_fij = 1 - 2 * Dfx - 2 * Dfy + Cfyij -r_ij_fuel

        # Fuel injection (as part of the mixture) 
        i_min = max(1, int(round(Ny / 3)) -1)          
        i_max = min(Ny - 2, int(round(2 * Ny / 3))+1 )  

        # Only to the left boundary:
        if j == 1 and i_min < i < i_max:
            A_f[k, k] = a_fij - Sbxf 
            M_f[k, k] = b_fij + Sbxf
        else:
            A_f[k, k] = a_fij
            M_f[k, k] = b_fij

        # Non-Diagonal Terms    
        if i != 1:     
            A_f[k, (i-2)*(Nx-2) + (j-1)] = -Sbyf
            M_f[k, (i-2)*(Nx-2) + (j-1)] =  Sbyf

        if i != Ny-2:  
            A_f[k, i*(Nx-2) + (j-1)] = -Smyf
            M_f[k, i*(Nx-2) + (j-1)] = Smyf

        if j != 1:    
            A_f[k, k-1] = -Sbxf
            M_f[k, k-1] =  Sbxf

        if j != Nx-2: 
            A_f[k, k+1] = -Smxf
            M_f[k, k+1] =  Smxf

# Function to fill matrices for O_2 transport (A_O and M_O)
def filling_matrixes_O(A_O,M_O,T, Yf,v_x,v_y):
    for k in range(Ni):
        i = 1 + k // (Nx-2)
        j = 1 + k % (Nx-2) 

        #Arrenhius Model
        K[i,j]=Ae*np.exp(-E_a/(R*T[i,j])) 

        #Coefficients for the discretized equation(oxidizer equation)
        DOx = D_Ox* Dt / (2 * (delta_x)**2)
        DOy = D_Oy * Dt / (2 * (delta_y)**2)
        COxij = -(v_x[i,j] * Dt) / (2 * delta_x)        
        COyij = -(v_y[i,j] * Dt) / (4 * delta_y)        

        SbxO = DOx - COxij    
        SmxO = DOx + 0*COxij  
        SbyO = DOy - COyij    
        SmyO = DOy + COyij    

        r_ij_O = ((Dt )/(2*rho))*K[i,j]*S_ratio*(rho**2/(MW_fuel))*Yf[i,j] 
        a_Oij = 1 + 2 * DOx + 2 * DOy - COxij +r_ij_O
        b_Oij = 1 - 2 * DOx - 2 * DOy + COyij -r_ij_O

        # Oxidizer injection limits
        i_min = max(1, int(round(Ny / 3)) -1)           
        i_max = min(Ny - 2, int(round(2 * Ny / 3))+1 )  

        if j == 1 and i_min < i < i_max:
            A_O[k, k] = a_Oij -SbxO 
            M_O[k, k] = b_Oij +SbxO
        else:
            A_O[k, k] = a_Oij
            M_O[k, k] = b_Oij

        #Non-Diagonal Terms    
        if i != 1:    
            A_O[k, (i-2)*(Nx-2) + (j-1)] = -SbyO
            M_O[k, (i-2)*(Nx-2) + (j-1)] =  SbyO

        if i != Ny-2:  
            A_O[k, i*(Nx-2) + (j-1)] = -SmyO
            M_O[k, i*(Nx-2) + (j-1)] = SmyO

        if j != 1:     
            A_O[k, k-1] = -SbxO
            M_O[k, k-1] =  SbxO

        if j != Nx-2:  
            A_O[k, k+1] = -SmxO
            M_O[k, k+1] =  SmxO


# 2. FUNCTIONS FOR INITIAL  AND BOUNDARY CONDITIONS (Species and Energy)
def define_initial_conditions_Yf():
    Yf_initial= np.zeros((Ny, Nx))
    return Yf_initial

def define_initial_conditions_YO():
    YO_initial= np.zeros((Ny, Nx))
    return YO_initial

def define_initial_conditions_T():
    xign = (Lx/9+2*Lx/9)/2    
    yign = 0.5*Ly           
    sigma = 0.01
    delta_T = 1200
    T_initial = np.ones((Ny, Nx))*320
    T_ignition = delta_T*np.exp(-((x-xign)**2+(y-yign)**2)/(sigma**2))  

    return T_initial, T_ignition

def apply_boundries_conditions_Yf(Yf):
    for k in range(Ni):
        i = 1 + k // (Nx-2)
        j = 1 + k % (Nx-2)
        i_min = max(1, int(np.round(Ny / 3)) - 1)           
        i_max = min(Ny - 2, int(np.round(2 * Ny / 3)) + 1) 
        if j == 1:
            if i_min < i < i_max:
                Yf[i, 0] = Yf[i, 1] + delta_x*(mfr_fuel/(v_x[i,j]*rho*Lx))  # Fuel injected from the left side
            else:
                Yf[0:i_min, 0] = Yf[i_max:-1, 0]= 0.0
        else:
            Yf[0,:]  = 0.0  
            Yf[:,-1] = 0.0 
            Yf[-1,:] = 0.0  
    return Yf

def stop_injection_Yf(Yf,step):
    for k in range(Ni):
        i = 1 + k // (Nx-2)
        j = 1 + k % (Nx-2)
        Yf[i, 0] = 0.0  
        Yf[0,:]  = 0.0  
        Yf[:,-1] = 0.0  
        Yf[-1,:] = 0.0  
    return Yf

def apply_boundries_conditions_YO(YO):
    for k in range(Ni):
        i = 1 + k // (Nx-2)
        j = 1 + k % (Nx-2)
        i_min = max(1, int(np.round(Ny / 3)) - 1)           
        i_max = min(Ny - 2, int(np.round(2 * Ny / 3)) + 1) 
        if j == 1:
            if i_min < i < i_max:
                YO[i, 0] = YO[i, 1] + delta_x*(mfr_O/(v_x[i,j]*rho*Lx)) # Fuel injected from the left side
            else:
                YO[0:i_min, 0] = YO[i_max:-1, 0]= 0.0
        else:
            YO[0,:]  = 0.0 
            YO[:,-1] = 0.0  
            YO[-1,:] = 0.0  
    return YO

def stop_injection_YO(YO, step):
    for k in range(Ni):
        i = 1 + k // (Nx-2)
        j = 1 + k % (Nx-2)
        YO[i, 0] = 0.0  
        YO[0,:]  = 0.0 
        YO[:,-1] = 0.0  
        YO[-1,:] = 0.0  
    return Yf

def apply_boundries_conditions_T(T):
    for k in range(Ni):
        i = 1 + k // (Nx-2)
        j = 1 + k % (Nx-2)
        i_min = max(1, int(np.round(Ny / 3)) - 1)
        i_max = min(Ny - 2, int(np.round(2 * Ny / 3)) + 1) 
        if j == 1:
            if i_min < i < i_max:
                T[i,0] = Te  
            else:
                T[i,0] = Ke + (1/Ce)*T[i,1]
        else:
            T[0,:]  = Ki + (1/Ci)*T[1, :]  
            T[:,-1] = Kd + (1/Cd)*T[i,-2]  
            T[-1,:] = Ks + (1/Cs)*T[-2,:]  
    return T


#3.MODELING THE SOURCE TERM Q (Heat Generation)
def modelate_source_term_Q(Yf, YO, Temp):
    #Modeling terms to be added in Q (BCs_Term)
    BCs_Term = np.zeros((Ny, Nx))
    for k in range(Ni):
        i = 1 + k // (Nx - 2)
        j = 1 + k % (Nx - 2)
        # Diffusion and convection coefficients of the enery discretized equation
        Dx = (Dt * alpha) / (2 * (delta_x)**2)        
        Dy = (Dt * 0.4*alpha) / (2 * (delta_y)**2)     
        Cxij = -(v_x[i,j] * Dt) / (2 * delta_x)       
        Cyij = -(v_y[i,j] * Dt) / (4 * delta_y)       

        # Coefficients for the discretized equation(energy equation)
        Sbx = Dx - Cxij    
        Smx = Dx + 0*Cxij  
        Sby = Dy - Cyij    
        Smy = Dy + Cyij   

        if i == 1:        
           if j == 1:         
                BCs_Term[i, j] = 2*(Sbx*Ke + Sby*Ki)
           elif j == Nx - 2:  
                BCs_Term[i, j] = 2*(Sby*Ki + Smx*Kd)
           else:              
                BCs_Term[i, j] = 2*Sby*Ki

        elif i == Ny - 2:  
            if j == 1:         
               BCs_Term[i, j] = 2*(Sbx*Ke + Smy*Ks)
            elif j == Nx - 2:  
               BCs_Term[i, j] = 2*(Smx*Kd + Smy*Ks)
            else:              
               BCs_Term[i, j] = 2*Smy*Ks 

        elif j == 1 and i != 1 and i != Ny - 2:      
            i_min = max(1, int(Ny / 3) - 1)
            i_max = min(Ny - 2, int(2 * Ny / 3) + 1)
            if i_min < i < i_max:
               BCs_Term[i, j] = 2*Sbx*Te 
            else:
               BCs_Term[i, j] = 2*Sbx*Ke

        elif j == Nx - 2 and i != 1 and i != Ny - 2:  
            BCs_Term[i, j] = 2*Smx*Kd

        else: 
            BCs_Term[i, j] = 0

    Q = np.zeros((Ny, Nx))  
    for k in range(Ni):
        i = 1 + k // (Nx - 2)
        j = 1 + k % (Nx - 2)
        q = Dt*Delta_H / (rho*MW_fuel*Cp)
        r_ij_energy = - (q)*(rho**2/(MW_fuel*MW_O2))*Ae*Yf[i,j]*YO[i,j]*( np.exp(-E_a / (R * Temp[i, j]))) #+ np.exp(-E_a / (R * T_n1[i, j])))
        Q[i, j] = r_ij_energy + BCs_Term[i,j]
        Q = np.maximum(Q, 0)
        Q = np.minimum(Q, 900)
    return BCs_Term, Q

#Source created by the discretization of species equation (Fuel)
BCs_Term_fuel = np.zeros((Ny, Nx))
for k in range(Ni):
    i = 1 + k // (Nx - 2)
    j = 1 + k % (Nx - 2)

    Dfx = D_fuelx* Dt / (2 * (delta_x)**2)
    Cfxij = -(v_x[i,j] * Dt) / (2 * delta_x)       
    Sbxf = Dfx - Cfxij    

    # Fuel injection Limits
    i_min = max(1, int(np.round(Ny / 3)) - 1)           
    i_max = min(Ny - 2, int(np.round(2 * Ny / 3)) + 1) 
    # Only for the left boundry and inside the corrected interval:
    if j == 1 and i_min < i < i_max:
        BCs_Term_fuel[i, j] = 2*Sbxf*(delta_x*mfr_fuel/(v_x[i,j]*rho*Lx))  
    else:
        BCs_Term_fuel[i, j] = 0

#Source created by the discretization of species equation (Oxidizer)
BCs_Term_O = np.zeros((Ny, Nx))
for k in range(Ni):
    i = 1 + k // (Nx - 2)
    j = 1 + k % (Nx - 2)
    
    DOx = D_Ox* Dt / (2 * (delta_x)**2)
    COxij = -(v_x[i,j] * Dt) / (2 * delta_x)       
    SbxO = DOx - COxij 

    # Oxidizer injection Limits
    i_min = max(1, int(np.round(Ny / 3)) - 1)           
    i_max = min(Ny - 2, int(np.round(2 * Ny / 3)) + 1) 
    # Only for the left boundry and inside the corrected interval:
    if j == 1 and i_min < i < i_max:
        BCs_Term_O[i, j] =2*SbxO*(delta_x*mfr_O/(v_x[i,j]*rho*Lx))
    else:
        BCs_Term_O[i, j] = 0


#. 4 SOLVE THE COUPLED SYSTEM FOR SPECIES (Yf, Yo) AND ENERGY (T) IN THE SAME LOOP

mPex = max_V_magnitude*Lx/D_fuelx          #peclet number for mass diffusion
tPex = max_V_magnitude*Lx/alpha            #peclet number for termal diffusion

from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

Yfs = []  
YOs = []  
Ymists = []
Ts = []   

ignition_triggered = False 
combustion_started = False 

# Initializing lists to store T and Y arrays in the time step loop
T_mean_list  = []
Yf_mean_list = []
YO_mean_list = []
Ymist_mean_ignition_list = []


#SIMULATION TIME STEP LOOP
for t in range(0,Nt):
    step = t
    print(f'step ={step}')

    #4.1
    if step == 0:
        # Initialize fuel fraction and temperature matrices
        Yf = define_initial_conditions_Yf()
        YO = define_initial_conditions_YO()
        Ymist= Yf+YO
        Yffraction_reached_in_ignition_region = 0.0 
        YOfraction_reached_in_ignition_region = 0.0
        Yffraction_reached_total = 0.0      
        YOfraction_reached_total = 0.0      
        T_initial , T_ignition  = define_initial_conditions_T()

        T = T_initial
        filling_matrixes_temperature(v_x,v_y)
        BCs_Term, Q = modelate_source_term_Q(Yf, YO, T) 
    #4.2
    else:
        ## Solve: Fuel and oxidezer injection kept or stopped
        #4.2.1
        if not combustion_started:                                    
            #1A  Solve  A_f·Yf[n+1] = M_f·Yf[n] for Yf at time n+1
            Yf_n = Yf
            # converting matrixes into sparse form
            A_f = np.zeros((Ni, Ni))
            M_f = np.zeros((Ni, Ni))
            filling_matrixes_f(A_f,M_f,T, YO, v_x,v_y)
            A_f_sparse = csc_matrix(A_f)
            M_f_sparse = csc_matrix(M_f)
            # convert internal values to vector form
            Yfn_vect = Yf[1:-1, 1:-1].reshape(Ni)
            BCs_Term_fuel_n_vect =  BCs_Term_fuel[1:-1, 1:-1].reshape(Ni)
            # solve
            b_f = M_f_sparse @ Yfn_vect +  BCs_Term_fuel_n_vect
            Yfnext = spsolve(A_f_sparse, b_f)
            Yf[1:-1, 1:-1] = Yfnext.reshape(Ny-2, Nx-2)
            Yf = apply_boundries_conditions_Yf(Yf)           
            Yf = np.maximum(Yf, 0)
            Yf = np.minimum(Yf, 1)
            Yf_n1= Yf

            #1B Solve  A_O·Y[n+1] = M_O·Y[n] for YO at time n+1
            # converting matrixes into sparse form
            YO_n = YO
            A_O = np.zeros((Ni, Ni))
            M_O = np.zeros((Ni, Ni))
            filling_matrixes_O(A_O,M_O,T, Yf, v_x,v_y)
            A_O_sparse = csc_matrix(A_O)
            M_O_sparse = csc_matrix(M_O)
            # convert internal values to vector form
            YOn_vect = YO[1:-1, 1:-1].reshape(Ni)
            BCs_Term_O_n_vect =  BCs_Term_O[1:-1, 1:-1].reshape(Ni)
            # solve
            b_O = M_O_sparse @ YOn_vect +  BCs_Term_O_n_vect
            YOnext = spsolve(A_O_sparse, b_O)
            YO[1:-1, 1:-1] = YOnext.reshape(Ny-2, Nx-2)
            YO = apply_boundries_conditions_YO(YO)              
            YO = np.maximum(YO, 0)
            YO = np.minimum(YO, 1)
            YO_n1 = YO

            #1C Temperature
            # Temperature variaton due to hot gas injection
            a=235.6       
            Tmisture = T_initial + (800-a)*Yf*np.ones((Ny,Nx)) + (800-a)*YO*np.ones((Ny,Nx))
            T = Tmisture
        #4.2.2
        else:                                                     
            #2A Solve  A_f·Yf[n+1] = M_f·Yf[n] for Yf at time n+1
            # converting matrixes into sparse form
            Yf_n = Yf
            Yf_n = np.maximum(Yf_n, 0)
            Yf_n = np.minimum(Yf_n, 1)

            A_f = np.zeros((Ni, Ni))
            M_f = np.zeros((Ni, Ni))
 
            #convection manually decreased until 0 to #reactants as they are  consumed          
            if step<390:
               filling_matrixes_f(A_f,M_f,T, YO,v_x,v_y) 
            elif step >= 390 and step < 440:
               filling_matrixes_f(A_f,M_f,T, YO,v_x/20,v_y/4)
            elif step >= 440 and step < 480:
                filling_matrixes_f(A_f,M_f,T, YO,v_x/100,v_y/20)
            else:
                 filling_matrixes_f(A_f,M_f,T, YO,0*v_x,0*v_y) 
            A_f_sparse = csc_matrix(A_f)
            M_f_sparse = csc_matrix(M_f)
            # convert internal values to vector form
            Yfn_vect = Yf[1:-1, 1:-1].reshape(Ni)
            # solve:
            b_f = M_f_sparse @ Yfn_vect  
            Yfnext = spsolve(A_f_sparse, b_f)
            Yf[1:-1, 1:-1] = Yfnext.reshape(Ny-2, Nx-2)
            Yf = stop_injection_Yf(Yf,step)              
            Yf = np.maximum(Yf, 0)
            Yf = np.minimum(Yf, 1)
            Yf_n1= Yf
            print(f'Yf_{t}')

            #2B Solve  A_O·Y[n+1] = M_O·Y[n] for YO at time n+1
            #2.1 converting matrixes into sparse form
            YO_n = YO
            YO_n = np.maximum(YO_n, 0)
            YO_n = np.minimum(YO_n, 1)
            A_O = np.zeros((Ni, Ni))
            M_O = np.zeros((Ni, Ni))
            
            #convection is manually decreased until 0 to #reactants as they are  consumed
            if step<390:
               filling_matrixes_O(A_O,M_O,T, Yf,v_x,v_y) 
            elif step >= 390 and step < 440:
               filling_matrixes_O(A_O,M_O,T, Yf,v_x/20,v_y/4)
            elif step >= 440 and step < 480:
                filling_matrixes_O(A_O,M_O,T, Yf,v_x/100,v_y/20)
            else:
                 filling_matrixes_O(A_O,M_O,T, Yf,0*v_x,0*v_y)  
            
            A_O_sparse = csc_matrix(A_O)
            M_O_sparse = csc_matrix(M_O)
            #2.2 convert internal values to vector form
            YOn_vect = YO[1:-1, 1:-1].reshape(Ni)
            #2.3 solve:
            b_O = M_O_sparse @ YOn_vect 
            YOnext = spsolve(A_O_sparse, b_O)
            YO[1:-1, 1:-1] = YOnext.reshape(Ny-2, Nx-2) 
            YO = stop_injection_YO(YO,step)                         
            YO = np.maximum(YO, 0)
            YO = np.minimum(YO, 1)
            YO = S_ratio*YO
            YO_n1 = YO

        Ymist= Yf+YO                                        #fraction of the mixture defined as the sum of fuel and oxidizer fractions 

        #4.2.3 Location, Parameters and conditions to start combustion
        # IGNITION CONDITION: Critical fraction of fuel and oxidizer in the ignition region to activate the resolution of energy equation
        Yf_critical = 0.062    
        YO_critical = 0.222
        # Ignation Region
        j_igmin = int(np.round(Nx/9)-1)
        j_igmax = int(np.round(2*Nx/9)+1)
        i_igmin = int(np.round(4*Ny/9)-2)
        i_igmax = int(np.round(5*Ny/9)+2)
        ignition_region = (slice(i_igmin, i_igmax), slice(j_igmin, j_igmax))  
        Yf_in_ignition = Yf[ignition_region]
        YO_in_ignition = YO[ignition_region]

        #4.2.4 Fuel  and oxidizer Avarage fraction in the ignition region and inside the chamber
        Yffraction_reached_in_ignition_region = np.mean(Yf_in_ignition)  
        YOfraction_reached_in_ignition_region = np.mean(YO_in_ignition)
        Yffraction_reached_total = np.mean(Yf)        
        YOfraction_reached_total = np.mean(YO)

        #4.2.5 conditions ignite the spark: "combustion_started" defines when combustion begins.
        if t>=2:
            if (Yffraction_reached_in_ignition_region >= Yf_critical and YOfraction_reached_in_ignition_region >= YO_critical) or combustion_started:   
                combustion_started = True
                print(f"Time step {t}: Yf_fraction reached(%) = {np.round(Yffraction_reached_in_ignition_region,4)}, YO_fraction reached(%) = {np.round(YOfraction_reached_in_ignition_region,4)}, Solve energy equation?= {combustion_started}")
                print(f"Time step {t}: Yf_fraction reached total(%) = {np.round(Yffraction_reached_total,4)}, YO_fraction reached total(%) = {np.round(YOfraction_reached_total,4)}, Solve energy equation?= {combustion_started}")       
            else:
                combustion_started = False
                print(f"Time step {t}: Yf_fraction reached(%) = {np.round(Yffraction_reached_in_ignition_region,4)}, YO_fraction reached(%) = {np.round(YOfraction_reached_in_ignition_region,4)}, Solve energy equation? = {combustion_started}")
                print(f"Time step {t}: Yf_fraction reached total(%) = {np.round(Yffraction_reached_total,4)}, YO_fraction reached total(%) = {np.round(YOfraction_reached_total,4)}, Solve energy equation?= {combustion_started}")
        
        #4.2.6 Ignite, Update the source term Q and Solve  A·T[n+1] = M·T[n] + Q for T at time n+1.
        if combustion_started:
            # Trigger ignition
            if not ignition_triggered:
                #Modeling the spark inition
                T = Tmisture + T_ignition       
                print(f'T_ignition{t}')
                ignition_triggered = True
            else:
                # Update the source term Q
                BCs_Term, Q = modelate_source_term_Q(Yf_n, YO_n, T)
                print(f'Q_{t}')

                # Updadte T solving energy equation
                # converting matrixes into sparse form
                A_sparse  = csc_matrix(A)      
                M_sparse  = csc_matrix(M)
                Tn_vect = T[1:-1, 1:-1].reshape(Ni)
                Qn_vect = Q[1:-1, 1:-1].reshape(Ni)
                #solving
                b_T = M_sparse @ Tn_vect + Qn_vect
                Tnext = spsolve(A_sparse, b_T)
                # convert internal values to vector form
                T[1:-1, 1:-1] = Tnext.reshape(Ny-2, Nx-2)
                T = apply_boundries_conditions_T(T)
                T = np.maximum(T, 100)
                T = np.minimum(T, 3000)
                print(f'T_{t}')

    #4.3 STORE VALUES TO PLOT LATER
    Yfs.append(np.copy(Yf))
    YOs.append(np.copy(YO))
    Ymists.append(np.copy(Ymist))
    Ts.append(np.copy(T))

    # time step description
    time_steps = np.arange(Nt)  # from 0 to 1405 steps

    # Average of T, Yf, YO and Ymist
    T_mean  =  np.mean(T)
    print(f'Mean temperature: {T_mean}')
    Yf_mean = Yffraction_reached_total       
    print(f'Mean Fuel: {Yf_mean}')
    YO_mean = YOfraction_reached_total
    print(f'Mean oxidizer: {YO_mean}')
    Ymist_mean_ignition = Yffraction_reached_in_ignition_region+YOfraction_reached_in_ignition_region 
    print(f'Mean_mist: {Ymist_mean_ignition}')

    # Storing the values in lists
    T_mean_list.append(T_mean)
    Yf_mean_list.append(Yf_mean)
    YO_mean_list.append(YO_mean)
    Ymist_mean_ignition_list.append(Ymist_mean_ignition)


    #4.4 BUILDING A PERSONALISED COLORMAP
    from matplotlib import pyplot as plt
    import matplotlib.colors as mcolors 

    if not combustion_started: 
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

        #Image resolution definations
        dpi = 100
        width_inches =1180/dpi
        height_inches =580/dpi
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
        plt.title('Temperature distribution')

        plt.subplot(1, 2, 2)                #MIXTURE FRACTION
        plt.imshow(Ymist, cmap=custom_cmapf, extent=[0, Lx, 0, Ly], vmin=Ymist_min, vmax=Ymist_max)
        plt.colorbar(label='')
        plt.xticks([0.0, Lx])              
        plt.yticks([0.0, Ly])
        plt.title('Air-fuel mixture distribution')

        # show time steps on screen
        plt.suptitle(f"Time step: {t}")


    else: 
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
        custom_cmapf = mcolors.LinearSegmentedColormap.from_list("custom_fuel", colorsf)

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

        #Image resolution definations
        dpi = 100
        width_inches =1180/dpi
        height_inches =580/dpi
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
        plt.title('Temperature distribution')

        plt.subplot(1, 2, 2)                #MIXTURE FRACTION
        plt.imshow(Ymist, cmap=custom_cmapf, extent=[0, Lx, 0, Ly], vmin=Ymist_min, vmax=Ymist_max)
        plt.colorbar(label='')
        plt.xticks([0.0, Lx])              
        plt.yticks([0.0, Ly])
        plt.title('Air-fuel mixture distributio')

        # show time step on screen
        plt.suptitle(f"Time step: {step}")


    #4.6 DISPAYING AND STORING IMAGES (To further use)
    if t%5==0 or t in [324,325, 326, 327, 328]:
        os.makedirs(f"Simulation3_images",exist_ok=True)
        plt.savefig(f"Simulation3_images/simulation_step_{step:04d}.png",dpi=150,bbox_inches='tight')
        plt.pause(5)
        plt.close()
        plt.show()
    plt.close()  

    #4.7 CLEANING MEMORY
    #clean RAM space to each 50th time step (if necessary)
    if t%50==0:
        gc.collect()
    if step >0:
        del A_f
        del M_f
        del A_O
        del M_O


# 5 SAVE THE T, Yf AND Yo DATA ARRAYS AS CSV (for further analysis)
np.savetxt("simulation3_data.csv",
        np.column_stack([time_steps, T_mean_list, Yf_mean_list, YO_mean_list, Ymist_mean_ignition_list]),
        delimiter=",",
        header="Time,T_mean,Yf_mean,YO_mean,Ymist_mean_ignition",
        comments='')

np.savez("simulacao3_data.npz", Yf_list=Yfs, YO_list=YOs, T_list=Ts,Ymist_list=Ymists)


