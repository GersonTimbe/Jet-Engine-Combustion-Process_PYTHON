# ===========================================================================
# Modeling the source term (Heat generation)
# ===========================================================================
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

    Q = np.zeros((Ny, Nx))                    # Initialize source term matrix
    for k in range(Ni):
        i = 1 + k // (Nx - 2)
        j = 1 + k % (Nx - 2)
        q = Dt*Delta_H / (rho*MW_fuel*Cp)
        r_ij_energy = - (q)*(rho**2/(MW_fuel*MW_O2))*Ae*Yf[i,j]*YO[i,j]*(np.exp(-E_a / (R * Temp[i, j])))
        Q[i, j] = r_ij_energy + BCs_Term[i,j]
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