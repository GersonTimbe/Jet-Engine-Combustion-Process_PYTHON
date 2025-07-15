# ===========================================================================
# Matrix Composition (constant velocity field) 
# ===========================================================================
# Initializing matrices A and M for the temperature equation
A = np.zeros((Ni, Ni))
M = np.zeros((Ni, Ni))

# Function to fill matrices A and M for the temperature equation
def filling_matrixes_temperature(v_x,v_y):
    for k in range(Ni):
        i = 1 + k // (Nx-2)
        j = 1 + k % (Nx-2)
        # Diffusion and convection coefficients(energy equation)
        Dx   = (Dt * alpha) / (2 * (delta_x)**2)         # Diffusion coefficient (x)
        Dy   = (Dt * 0.4*alpha) / (2 * (delta_y)**2)     # Diffusion coefficient (y)
        Cxij = -(v_x[i,j] * Dt) / (2 * delta_x)          # Convection coefficient (x)
        Cyij = -(v_y[i,j] * Dt) / (4 * delta_y)          # Convection coefficient (y)

        # Coefficients for the discretized equation(energy equation)
        Sbx = Dx - Cxij                                  # T(i,j-1) coefficient
        Smx = Dx + 0*Cxij                                # T(i,j+1) coefficient
        Sby = Dy - Cyij                                  # T(i-1,j) coefficient
        Smy = Dy + Cyij                                  # T(i+1,j) coefficient

        a = 1 + 2 * Dx + 2 * Dy - Cxij                   # T(i,j) coefficient for t=n+1
        b = 1 - 2 * Dx - 2 * Dy + Cxij                   # T(i,j) coefficient for t=n
        #Main Diagonal Modified by BCi:
        if i == 1:                            # bottom boundary
            if j == 1:                              # bottom left corner canto inferior esquerdo
                A[k, k] =  a - Sby / Ci - Sbx/Ce
                M[k, k] =  b + Sby / Ci + Sbx/Ce
            elif j == Nx - 2:                       # bottom right corner canto inferior direito
                A[k, k] =  a - Sby / Ci - Smx / Cd
                M[k, k] =  b + Sby / Ci + Smx / Cd
            else:                                   # other nodes of the bottom boundary
                A[k, k] =  a - Sby / Ci
                M[k, k] =  b + Sby / Ci

        elif i == Ny - 2:                     # upper boundary
            if j == 1:                             # upper left corner
                A[k, k] = a - Smy / Cs - Sbx/Ce
                M[k, k] = b + Smy / Cs + Sbx/Ce
            elif j == Nx - 2:                      # upper right corner
                A[k, k] = a - Smy / Cs - Smx / Cd
                M[k, k] = b + Smy / Cs + Smx / Cd
            else:                                  # other nodes of the upper boundary
                A[k, k] = a - Smy / Cs
                M[k, k] = b + Smy / Cs

        elif j == 1 and i!=1 and i!=Ny-2:    # left boundary (except corners)
            i_min = max(1, int(round(Ny / 3))-1)            
            i_max = min(Ny - 2, int(round(2 * Ny / 3))+1)  
            if i_min < i < i_max:
                A[k, k] = a
                M[k, k] = b
            else:
                A[k, k] = a - Sbx/Ce
                M[k, k] = b + Sbx/Ce
        elif j == Nx - 2 and i!=1 and i!=Ny-2:  # right boundary (except corners)
             A[k, k] = a - Smx / Cd
             M[k, k] = b + Smx / Cd

        else:                                  # internal nodes
             A[k, k] = a
             #A[k, k] = 9
             M[k, k] = b

        #Non-Diagonal Terms    
        if i != 1:                             # Lower lone diagonal
            A[k, (i-2)*(Nx-2) + (j-1)] = -Sby
            M[k, (i-2)*(Nx-2) + (j-1)] =  Sby

        if i != Ny-2:                          # Upper lone diagonal
            A[k, i*(Nx-2) + (j-1)] = -Smy
            M[k, i*(Nx-2) + (j-1)] =  Smy

        if j != 1:                             # Lower main diagonal
            A[k, k-1] = -Sbx
            M[k, k-1] =  Sbx

        if j != Nx-2:                          # Upper main diagonal
            A[k, k+1] = -Smx
            #A[k, k+1] = 11
            M[k, k+1] =  Smx

# Initializing matrices for fuel and oxydizer transport (A_f, A_O, M_f and M_O)
A_f = np.zeros((Ni, Ni))
M_f = np.zeros((Ni, Ni))
A_O = np.zeros((Ni, Ni))
M_O = np.zeros((Ni, Ni))
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
        Cfxij = -(v_x[i,j] * Dt) / (2 * delta_x)        # Convection coeffiCsent (x)
        Cfyij = -(v_y[i,j] * Dt) / (4 * delta_y)        # Convection coeffiCsent (y)

        Sbxf = Dfx - Cfxij                              # Yf(i,j-1) coeffiCsent
        Smxf = Dfx + 0*Cfxij                            # Yf(i,j+1) coeffiCsent
        Sbyf = Dfy - Cfyij                              # Yf(i-1,j) coeffiCsent
        Smyf = Dfy + Cfyij
 
        #Coefficients for fuel equation(2)
        r_ij_fuel = ((Dt )/(2*rho))*K[i,j]*(rho**2/(MW_O2))*YO[i,j]  # reaction rate term
        a_fij = 1 + 2 * Dfx + 2 * Dfy - Cfxij +r_ij_fuel
        b_fij = 1 - 2 * Dfx - 2 * Dfy + Cfyij -r_ij_fuel

        # Fuel injection limits 
        i_min = max(1, int(round(Ny / 3)) -1)          
        i_max = min(Ny - 2, int(round(2 * Ny / 3))+1 )  

        # Only to the left boundary and inside the corrected range 
        if j == 1 and i_min < i < i_max:
            A_f[k, k] = a_fij - Sbxf 
            M_f[k, k] = b_fij + Sbxf
        else:
            A_f[k, k] = a_fij
            M_f[k, k] = b_fij

        #Non-Diagonal Terms    
        if i != 1:                                    # Lower lone diagonal
            A_f[k, (i-2)*(Nx-2) + (j-1)] = -Sbyf
            M_f[k, (i-2)*(Nx-2) + (j-1)] =  Sbyf

        if i != Ny-2:                                 # Upper lone diagonal
            A_f[k, i*(Nx-2) + (j-1)] = -Smyf
            M_f[k, i*(Nx-2) + (j-1)] = Smyf

        if j != 1:                                    # Lower main diagonal
            A_f[k, k-1] = -Sbxf
            M_f[k, k-1] =  Sbxf

        if j != Nx-2:                                 # Upper main diagonal
            A_f[k, k+1] = -Smxf
            M_f[k, k+1] =  Smxf

# Function to fill matrices for O_2 transport (A_O and M_O)
def filling_matrixes_O(A_O,M_O,T, Yf,v_x,v_y):
    for k in range(Ni):
        i = 1 + k // (Nx-2)
        j = 1 + k % (Nx-2) 
        K[i,j]=Ae*np.exp(-E_a/(R*T[i,j])) 

        #CoeffiCsents for the discretized equation(oxidizer equation)
        DOx = D_Ox* Dt / (2 * (delta_x)**2)
        DOy = D_Oy * Dt / (2 * (delta_y)**2)
        COxij = -(v_x[i,j] * Dt) / (2 * delta_x)        # Convection coefficient (x)
        COyij = -(v_y[i,j] * Dt) / (4 * delta_y)        # Convection coefficient (y)

        SbxO = DOx - COxij                              # O(i,j-1) coefficient
        SmxO = DOx + 0*COxij                            # O(i,j+1) coefficient
        SbyO = DOy - COyij                              # O(i-1,j) coefficient
        SmyO = DOy + COyij                              # O(i+1,j) coefficient

        r_ij_O = ((Dt )/(2*rho))*K[i,j]*S_ratio*(rho**2/(MW_fuel))*Yf[i,j]  # Reaction rate term
        a_Oij = 1 + 2 * DOx + 2 * DOy - COxij +r_ij_O
        b_Oij = 1 - 2 * DOx - 2 * DOy + COyij -r_ij_O

        # Cálculo ajustado dos índices (Fuel injection limits Limites de injeccao do Fuel)
        i_min = max(1, int(round(Ny / 3)) -1)           
        i_max = min(Ny - 2, int(round(2 * Ny / 3))+1 )  

        # Apenas para a borda esquerda e dentro do intervalo corrigido:
        if j == 1 and i_min < i < i_max:
            A_O[k, k] = a_Oij -SbxO 
            M_O[k, k] = b_Oij +SbxO
        else:
            A_O[k, k] = a_Oij
            M_O[k, k] = b_Oij

        #Non-Diagonal Terms    
        if i != 1:                                        # Lower lone diagonal
            A_O[k, (i-2)*(Nx-2) + (j-1)] = -SbyO
            M_O[k, (i-2)*(Nx-2) + (j-1)] =  SbyO

        if i != Ny-2:                                     # Upper lone diagonal
            A_O[k, i*(Nx-2) + (j-1)] = -SmyO
            M_O[k, i*(Nx-2) + (j-1)] = SmyO

        if j != 1:                                        # Lower main diagonal
            A_O[k, k-1] = -SbxO
            M_O[k, k-1] =  SbxO

        if j != Nx-2:                                     # Upper main diagonal
            A_O[k, k+1] = -SmxO
            M_O[k, k+1] =  SmxO



# ===========================================================================
# Matrix Composition (variable velocity field)  
# ===========================================================================
# Initializing matrices A and M for the  Energy equation
A = np.zeros((Ni, Ni))
M = np.zeros((Ni, Ni))

# Function to fill matrices A and M
def filling_matrixes_temperature(v_x,v_y):
    for k in range(Ni):
        i = 1 + k // (Nx-2)
        j = 1 + k % (Nx-2)
        # Diffusion and convection coefficients
        Dx = (Dt * alpha) / (2 * (delta_x)**2)  
        Dy = (Dt * 0.4*alpha) / (2 * (delta_y)**2)  
        Cxij = -(v_x[i,j] * Dt) / (2 * delta_x)  
        Cyij = -(v_y[i,j] * Dt) / (4 * delta_y) 

        # Coefficient for the discretized equation
        Sbx = Dx - Cxij   
        Smx = Dx + 0*Cxij 
        Sby = Dy - Cyij   
        Smy = Dy + Cyij  
        a = 1 + 2 * Dx + 2 * Dy - Cxij  
        b = 1 - 2 * Dx - 2 * Dy + Cxij 
        #Main diagonal:
        if i == 1:                                 # bottom boundary 
            if j == 1:                               #bottom-left boudary
                A[k, k] =  a - Sby / Ci - Sbx/Ce
                M[k, k] =  b + Sby / Ci + Sbx/Ce
            elif j == Nx - 2:                        # bottom-right boudary
                A[k, k] =  a - Sby / Ci - Smx / Cd
                M[k, k] =  b + Sby / Ci + Smx / Cd
            else:                                    # other nodes of the bottom boundory 
                A[k, k] =  a - Sby / Ci 
                M[k, k] =  b + Sby / Ci

        elif i == Ny - 2:                          # upper boundry
            if j == 1:                              # upper-left corner
                A[k, k] = a - Smy / Cs - Sbx/Ce
                M[k, k] = b + Smy / Cs + Sbx/Ce
            elif j == Nx - 2:                       # upper-right corner
                A[k, k] = a - Smy / Cs - Smx / Cd
                M[k, k] = b + Smy / Cs + Smx / Cd
            else:                                   # other nodes of the upper boundary 
                A[k, k] = a - Smy / Cs
                M[k, k] = b + Smy / Cs

        elif j == 1 and i!=1 and i!=Ny-2:        # left boundary (except corners)
            i_min = max(1, int(round(Ny / 3))-1)     
            i_max = min(Ny - 2, int(round(2 * Ny / 3))+1)  
            if i_min < i < i_max:
                A[k, k] = a
                M[k, k] = b
            else:
                A[k, k] = a - Sbx/Ce
                M[k, k] = b + Sbx/Ce
        elif j == Nx - 2 and i!=1 and i!=Ny-2:  # right boundary (except corners)
             A[k, k] = a - Smx / Cd
             M[k, k] = b + Smx / Cd

        else:                                  # internal nodes
             A[k, k] = a
             M[k, k] = b

        #Non-Diagonal Terms    
        if i != 1:                             # Lower lone diagonal
            A[k, (i-2)*(Nx-2) + (j-1)] = -Sby
            M[k, (i-2)*(Nx-2) + (j-1)] =  Sby

        if i != Ny-2:                          # Upper lone diagonal
            A[k, i*(Nx-2) + (j-1)] = -Smy
            M[k, i*(Nx-2) + (j-1)] =  Smy

        if j != 1:                             # Lower main diagonal
            A[k, k-1] = -Sbx
            M[k, k-1] =  Sbx

        if j != Nx-2:                          # Upper main diagonal
            A[k, k+1] = -Smx
            M[k, k+1] =  Smx

# Initializing matrices for F and Oxidizer transport (A_f, A_O, M_f and M_O)
A_f = np.zeros((Ni, Ni))
M_f = np.zeros((Ni, Ni))
A_O = np.zeros((Ni, Ni))
M_O = np.zeros((Ni, Ni))

# Function to fill matrices for fuel transport (A_f and M_f)
def filling_matrixes_f(A_f,M_f,T, YO,v_x,v_y):
    for k in range(Ni):
        i = 1 + k // (Nx-2)
        j = 1 + k % (Nx-2)
        #Arrenhius Model
        K[i,j]=Ae*np.exp(-E_a/(R*T[i,j]))

        #Coefficient for the discretized equation(fuel equation)
        Dfx = D_fuelx* Dt / (2 * (delta_x)**2)
        Dfy = D_fuely * Dt / (2 * (delta_y)**2)
        Cfxij = -(v_x[i,j] * Dt) / (2 * delta_x)        
        Cfyij = -(v_y[i,j] * Dt) / (4 * delta_y)        

        Sbxf = Dfx - Cfxij  
        Smxf = Dfx + 0*Cfxij  
        Sbyf = Dfy - Cfyij  
        Smyf = Dfy + Cfyij
 
        #Coefficients for fuel equation(2)
        #Reaction rate term
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

        #Non-Diagonal Terms    
        if i != 1:                                   # Lower lone diagonal
            A_f[k, (i-2)*(Nx-2) + (j-1)] = -Sbyf
            M_f[k, (i-2)*(Nx-2) + (j-1)] =  Sbyf

        if i != Ny-2:                                # Upper lone diagonal
            A_f[k, i*(Nx-2) + (j-1)] = -Smyf
            M_f[k, i*(Nx-2) + (j-1)] = Smyf

        if j != 1:                                   # Lower main diagonal
            A_f[k, k-1] = -Sbxf 
            M_f[k, k-1] =  Sbxf

        if j != Nx-2:                                # Upper main diagonal
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
        
        #Reaction rate term
        r_ij_O = ((Dt )/(2*rho))*K[i,j]*S_ratio*(rho**2/(MW_fuel))*Yf[i,j]  
        a_Oij = 1 + 2 * DOx + 2 * DOy - COxij +r_ij_O
        b_Oij = 1 - 2 * DOx - 2 * DOy + COyij -r_ij_O

        # Oxidizer injection (as part of the mixture) 
        i_min = max(1, int(round(Ny / 3)) -1)           
        i_max = min(Ny - 2, int(round(2 * Ny / 3))+1 ) 

        # only the left boundary :
        if j == 1 and i_min < i < i_max:
            A_O[k, k] = a_Oij -SbxO 
            M_O[k, k] = b_Oij +SbxO
        else:
            A_O[k, k] = a_Oij
            M_O[k, k] = b_Oij

        #Non-Diagonal Terms    
        if i != 1:                                   # Lower lone diagonal
            A_O[k, (i-2)*(Nx-2) + (j-1)] = -SbyO
            M_O[k, (i-2)*(Nx-2) + (j-1)] =  SbyO

        if i != Ny-2:                                # Upper lone diagonal
            A_O[k, i*(Nx-2) + (j-1)] = -SmyO
            M_O[k, i*(Nx-2) + (j-1)] = SmyO

        if j != 1:                                   # Lower main diagonal
            A_O[k, k-1] = -SbxO
            M_O[k, k-1] =  SbxO

        if j != Nx-2:                                # Upper main diagonal
            A_O[k, k+1] = -SmxO
            M_O[k, k+1] =  SmxO
