# ===========================================================================
# Initial and Boundary conditions (Species and Energy)
# ===========================================================================
def define_initial_conditions_Yf():
    Yf_initial= np.zeros((Ny, Nx))
    return Yf_initial

def define_initial_conditions_YO():
    YO_initial= np.zeros((Ny, Nx))
    return YO_initial

def define_initial_conditions_T():
    xign = (Lx/9+2*Lx/9)/2        # center of the ignation region
    yign = 0.5*Ly                 
    sigma = 0.01
    delta_T = 1200
    T_initial = np.ones((Ny, Nx))*320
    #euler function to modelate the ignition
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
                Yf[i, 0] = Yf[i, 1] + delta_x*(mfr_fuel/(v_x[i,j]*rho*Lx))  
            else:
                Yf[0:i_min, 0] = Yf[i_max:-1, 0]= 0.0
        else:
            Yf[0,:]  = 0.0              # No fuel on the bottom boundary
            Yf[:,-1] = 0.0              # No fuel at the right boundary 
            Yf[-1,:] = 0.0              # No fuel at the upper boundary
    return Yf

def stop_injection_Yf(Yf,step):         #Fuel injection stopped after ignition
    for k in range(Ni):
        i = 1 + k // (Nx-2)
        j = 1 + k % (Nx-2)
        Yf[i, 0] = 0.0                  # injection stopped
        Yf[0,:]  = 0.0                  # No fuel at the Top side
        Yf[:,-1] = 0.0                  # No fuel at the right side
        Yf[-1,:] = 0.0                  # No fuel at the Bottom side
    return Yf

def apply_boundries_conditions_YO(YO): 
    for k in range(Ni):
        i = 1 + k // (Nx-2)
        j = 1 + k % (Nx-2)
        i_min = max(1, int(np.round(Ny / 3)) - 1)           
        i_max = min(Ny - 2, int(np.round(2 * Ny / 3)) + 1) 
        if j == 1:
            if i_min < i < i_max:     # Oxidizer injected from the left boundary 
                YO[i, 0] = YO[i, 1] + delta_x*(mfr_O/(v_x[i,j]*rho*Lx))
            else:
                YO[0:i_min, 0] = YO[i_max:-1, 0]= 0.0
        else:
            YO[0,:]  = 0.0            # No fuel at the Top boundary
            YO[:,-1] = 0.0            # No fuel at the right boundarvy
            YO[-1,:] = 0.0            # No fuel at the Bottom boundary
    return YO

def stop_injection_YO(YO, step):      #Oxidizer injection stopped after ignition
    for k in range(Ni):
        i = 1 + k // (Nx-2)
        j = 1 + k % (Nx-2)
        YO[i, 0] = 0.0                # injection stopped
        YO[0,:]  = 0.0                # No fuel at the upper boundary 
        YO[:,-1] = 0.0                # No fuel at the right boundary 
        YO[-1,:] = 0.0                # No fuel at the Bottom boundary
    return Yf

def apply_boundries_conditions_T(T):
    for k in range(Ni):
        i = 1 + k // (Nx-2)
        j = 1 + k % (Nx-2)
        i_min = max(1, int(np.round(Ny / 3)) - 1)
        i_max = min(Ny - 2, int(np.round(2 * Ny / 3)) + 1) 
        if j == 1:
            if i_min < i < i_max:
                T[i,0] = Te                    #high temperature of the gas flowing to inside the chamber
            else:
                T[i,0] = Ke + (1/Ce)*T[i,1]
        else:
            T[0,:]  = Ki + (1/Ci)*T[1, :]     # Temperature at the bottom boundary
            T[:,-1] = Kd + (1/Cd)*T[i,-2]     # Temperature at the Right boundary
            T[-1,:] = Ks + (1/Cs)*T[-2,:]     # Temperature at the upper boundary
    return T