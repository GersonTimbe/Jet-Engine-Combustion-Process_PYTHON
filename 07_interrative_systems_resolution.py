# ===========================================================================
# Interrative system resolution
# ===========================================================================
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

Yfs = []                      # To store the evolution of fuel fraction
YOs = []                      # To store the evolution of oxidizer fraction
Ymists = []
Ts = []                       # To store the evolution of temperature

# Time-stepping loop
ignition_triggered = False    #To ignite the mixture
combustion_started = False    #To stop air+fuel injection after ignition

# Avarage of the parameters of interest
T_mean_list  = []
Yf_mean_list = []
YO_mean_list = []
Ymist_mean_ignition_list = []

for t in range(0,Nt):
    step = t

    if step == 0:
        # Initialize fuel fraction and temperature matrices
        Yf = define_initial_conditions_Yf()
        YO = define_initial_conditions_YO()
        Ymist= Yf+YO
        Yffraction_reached_in_ignition_region = 0.0 
        YOfraction_reached_in_ignition_region = 0.0
        Yffraction_reached_total = 0.0                         # avarage fuel fraction inside the chamber
        YOfraction_reached_total = 0.0                         # avarage oxidizer fraction  inside the chamber
        T_initial , T_ignition  = define_initial_conditions_T()

        T = T_initial
        filling_matrixes_temperature(v_x,v_y)
        BCs_Term, Q = modelate_source_term_Q(Yf, YO, T) 
 
    else:
        ## Solve: Fuel and oxidizer injection kept or stopped
        if not combustion_started:                  #INJECTION KEPT
            #1.a  Solve  A_f·Yf[n+1] = M_f·Yf[n] for Yf at time n+1
            Yf_n = Yf
            #1.1 converting matrixes into sparse form
            A_f = np.zeros((Ni, Ni))
            M_f = np.zeros((Ni, Ni))
            filling_matrixes_f(A_f,M_f,T, YO, v_x,v_y)      
            #filling_matrixes_f(T, YO,v_x,v_y)
            A_f_sparse = csc_matrix(A_f)
            M_f_sparse = csc_matrix(M_f)
            #1.2 convert internal values to vector form
            Yfn_vect = Yf[1:-1, 1:-1].reshape(Ni)
            BCs_Term_fuel_n_vect =  BCs_Term_fuel[1:-1, 1:-1].reshape(Ni)
            #1.3 solve
            b_f = M_f_sparse @ Yfn_vect +  BCs_Term_fuel_n_vect
            Yfnext = spsolve(A_f_sparse, b_f)
            Yf[1:-1, 1:-1] = Yfnext.reshape(Ny-2, Nx-2)
            Yf = apply_boundries_conditions_Yf(Yf)           
            Yf_n1= Yf

            #2.a Solve  A_O·Y[n+1] = M_O·Y[n] for YO at time n+1
            #2.1 converting matrixes into sparse form
            YO_n = YO
            A_O = np.zeros((Ni, Ni))
            M_O = np.zeros((Ni, Ni))
            filling_matrixes_O(A_O,M_O,T, Yf, v_x,v_y)
            A_O_sparse = csc_matrix(A_O)
            M_O_sparse = csc_matrix(M_O)
            #2.2 convert internal values to vector form
            YOn_vect = YO[1:-1, 1:-1].reshape(Ni)
            BCs_Term_O_n_vect =  BCs_Term_O[1:-1, 1:-1].reshape(Ni)
            #2.3 solve
            b_O = M_O_sparse @ YOn_vect +  BCs_Term_O_n_vect
            YOnext = spsolve(A_O_sparse, b_O)
            YO[1:-1, 1:-1] = YOnext.reshape(Ny-2, Nx-2)
            YO = apply_boundries_conditions_YO(YO)              
            YO_n1 = YO

            #3.a Temperature
            #3.1 Temperature variation due to hot gas injection
            a=235.6        #parameter to adjust the gases inlet temperature 
            Tmixture = T_initial+(800-a)*Yf*np.ones((Ny,Nx))+(800-a)*YO*np.ones((Ny,Nx))
            T = Tmixture
            

        else:                                  #INJECTION STOPPED
            #1b Solve  A_f·Yf[n+1] = M_f·Yf[n] for Yf at time n+1
            #1.1 converting matrixes into sparse form
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
            #1.2 convert internal values to vector form
            Yfn_vect = Yf[1:-1, 1:-1].reshape(Ni)
            #1.3 solve:
            b_f = M_f_sparse @ Yfn_vect  
            Yfnext = spsolve(A_f_sparse, b_f)
            Yf[1:-1, 1:-1] = Yfnext.reshape(Ny-2, Nx-2)
            Yf = stop_injection_Yf(Yf,step)              
            Yf_n1= Yf

            #2b Solve  A_O·Y[n+1] = M_O·Y[n] for YO at time n+1
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
            YO = S_ratio*YO
            YO_n1 = YO

        Ymist= Yf+YO
        Ymist = np.maximum(Ymist, 0)
        Ymist = np.minimum(Ymist, 1)

        #3b Location, Parameters and conditions to start combustion
        #3.1 Critical fraction of fuel and oxidizer in the ignition region to activate the combustion
        Yf_critical = 0.062    
        YO_critical = 0.222
        # 3.2 Ignition Region
        j_igmin = int(np.round(Nx/9)-1)
        j_igmax = int(np.round(2*Nx/9)+1)
        i_igmin = int(np.round(4*Ny/9)-2)
        i_igmax = int(np.round(5*Ny/9)+2)
        ignition_region = (slice(i_igmin, i_igmax), slice(j_igmin, j_igmax))  
        Yf_in_ignition = Yf[ignition_region]
        YO_in_ignition = YO[ignition_region]

        # 3.3 Fuel  and oxidizer Average fraction in the ignition region and inside the chamber
        Yffraction_reached_in_ignition_region = np.mean(Yf_in_ignition)  
        YOfraction_reached_in_ignition_region = np.mean(YO_in_ignition)
        Yffraction_reached_total = np.mean(Yf)    # avarage fuel fraction inside the chamber
        YOfraction_reached_total = np.mean(YO)    # avarage oxidizer fraction inside the chamber 

        #3.4 conditions to ignite the spark: "combustion_started" defines when combustion begins.
        if t>=2:
            if (Yffraction_reached_in_ignition_region >= Yf_critical and
            YOfraction_reached_in_ignition_region >= YO_critical) or combustion_started:   
                combustion_started = True       
            else:
                combustion_started = False
        
        #4.3 Ignite, Update Q and Solve A·T[n+1]=M·T[n] + Q for T at time n+1
        if combustion_started:
            #4.3.1 Trigger ignition
            if not ignition_triggered:
                #Modeling the spark inition
                T = Tmisture + T_ignition       
                ignition_triggered = True
            else:
                #4.3.2 Update the source term Q
                BCs_Term, Q = modelate_source_term_Q(Yf_n, YO_n, T)

                #4.3.3 Updadte T solving energy equation
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
          
    ##5 Store values to plotting later
    Yfs.append(np.copy(Yf))
    YOs.append(np.copy(YO))
    Ymists.append(np.copy(Ymist))
    Ts.append(np.copy(T))

    # Average of T, Yf, YO and Ymist
    T_mean  =  np.mean(T)
    Yf_mean = Yffraction_reached_total 
    YO_mean = YOfraction_reached_total
    Ymist_mean_in_ignition_region = Yffraction_reached_in_ignition_region + YOfraction_reached_in_ignition_region 

    # Storing the values in lists
    T_mean_list.append(T_mean)
    Yf_mean_list.append(Yf_mean)
    YO_mean_list.append(YO_mean)
    Ymist_mean_ignition_list.append(Ymist_mean_ignition)
