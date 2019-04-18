import numpy as np, warnings, seawater as sw, warnings, matplotlib.pyplot as plt

# main, time-stepping function
def model_timestep(T,S,U,V,z,I,L,E,P,tau_x,tau_y,dt,nabla_b=None,Ekman_Q_flux=None,use_static_stability=True,use_mixed_layer_stability=True,use_shear_stability=True,use_Ekman_flux=False,use_MLI=False,tracer=None,vert_diffusivity=None,verbose=False,I1=0.62,I2=None,lambda1=.6,lambda2=20,T0=17,S0=36,rho0=None,alpha=None,beta=None,f=sw.f(40),return_MLD=False):
    # define initial variables
    c = 4218 # heat capacity (J kg^-1 K^-1)
    if I2 is None:
        I2 = 1-I1
    if I1+I2 != 1:
        raise Exception('Shortwave insolation amplitudes need to sum to unity')
    if rho0 is None:
        rho0 = sw.dens(S0,T0,0)
    if alpha is None:
        alpha = -sw.alpha(S0,T0,0)*rho0 # multiply by rho to fit into nice equation of state (see get_rho function)
    if beta is None:
        beta = sw.beta(S0,T0,0)*rho0 # # multiply by rho to fit into nice equation of state (see get_rho function)
    
    dz = z[1]-z[0]
    
    if use_Ekman_flux and Ekman_Q_flux is None:
        raise Exception('Using Ekman-induced buoyacy flux but no buoyancy gradients were given.')
    if use_MLI and nabla_b is None:
        raise Exception('Using MLI but no horizontal buoyancy gradient given.')
    
    T = T.copy(); S = S.copy(); U = U.copy(); V = V.copy(); # so don't overwrite data
    if tracer is not None:
        tracer = tracer.copy();
    
    # make initial heat profile
    I_profile = -I/dz*(I1*np.exp(-z/lambda1)*(np.exp(-dz/2/lambda1)-np.exp(dz/2/lambda1))+I2*np.exp(-z/lambda2)*(np.exp(-dz/2/lambda2)-np.exp(dz/2/lambda2)))
    L_profile = np.zeros(len(z)); L_profile[0] = L/dz
    Q_profile = I_profile+L_profile
            
    if use_Ekman_flux:
        A = 0.1 # eddy viscosity m^2 s^-1
        z_Ek = np.sqrt(A/np.abs(f))
        if verbose:
            print('Using Ekman depth of %d m'%z_Ek)
        z_Ek_ind = np.where(z>z_Ek)[0][0]
        Q_Ek_profile = np.zeros(len(z))
        Q_Ek_profile[0:z_Ek_ind] = Ekman_Q_flux/z_Ek*dz; 
        Q_profile += Q_Ek_profile
    
    
    # update temperature
    dTdt = Q_profile/(c*rho0)
    if use_MLI:
        mld_ind = get_mld_ind(T,S,U,V,z)
        mld = z[mld_ind]
        C_e = 0.06; 
        g = 9.81; # gravitational acceleration (m s^-2)
        c = 4218; # heat capacity (J kg^-1 degC^-1
        MLI_dTdt = -C_e*nabla_b**2*mld**2*rho0/(np.abs(f)*alpha*g)
        vert_profile = 4/mld*(1-2*z/mld)*(16+10*(1-2*z/mld)**2)/21 # this is vertical derivative of mu(z)
        vert_profile[mld_ind::] = 0 
        vert_profile[0:mld_ind] -= np.mean(vert_profile[0:mld_ind]) # just to ensure that no heat added to system
        dTdt += MLI_dTdt*vert_profile
    T += dTdt * dt
    
    # update salinity
    dSdt = S[0]*(E-P)/dz/1000
    S[0] += dSdt * dt
    
    if use_MLI:
        rho = get_rho(T,S,T0,S0,rho0,alpha,beta)
        half_mld_ind = int(mld_ind/2)
        if np.any(np.diff(rho[half_mld_ind::])<0):
            if verbose:
                print('Need to homogenize discontinuity at base of previous mixed layer')
            # get rid of discontinuity at base of previous mixed layer
            # homogenize the region of water from mld/2 to z*
            # z* is the shallowest value (> mld/2) such that the homogenized rho <= rho(z*)
            zstar_ind = mld_ind.copy()
            while np.mean(rho[half_mld_ind:zstar_ind])>=rho[zstar_ind]:
                if verbose:
                    print('Deepening z*...')
                zstar_ind += 1
            T[half_mld_ind:zstar_ind] = np.mean(T[half_mld_ind:zstar_ind])
            S[half_mld_ind:zstar_ind] = np.mean(S[half_mld_ind:zstar_ind])
            if tracer is not None:
                tracer[:,half_mld_ind:zstar_ind] = np.mean(tracer[:,half_mld_ind:zstar_ind],axis=1)
        elif verbose:
            print('No need to homogenize base of previous mixed layer')
        
    
    # update momentum
    # first rotate momentum halfway
    angle = -f*dt/2 # currently assuming this is in rad
    U,V = rotate(angle,U,V)
    # then add wind stress
    mld_ind = get_mld_ind(T,S,U,V,z)
    mld = z[mld_ind]
    U[0:mld_ind] += tau_x/mld/rho0*dz * dt; 
    V[0:mld_ind] += tau_y/mld/rho0*dz * dt;
    # then rotate second half
    U,V = rotate(angle,U,V)
    
    if use_static_stability:
        T,S,U,V,tracer = static_stability(T,S,U,V,z,T0,S0,rho0,alpha,beta,tracer=tracer,verbose=verbose)
    if use_mixed_layer_stability:
        T,S,U,V,tracer = mixed_layer_stability(T,S,U,V,z,T0,S0,rho0,alpha,beta,tracer=tracer,verbose=verbose)
    if use_shear_stability:
        T,S,U,V,tracer = shear_stability(T,S,U,V,z,T0,S0,rho0,alpha,beta,tracer=tracer,verbose=verbose)
    
    if vert_diffusivity is not None:
        dTdt = np.zeros(len(T))
        dTdt[1:-1] = np.diff(np.diff(T))/dz**2
        T += vert_diffusivity*dTdt*dt
        
        dSdt = np.zeros(len(S))
        dSdt[1:-1] = np.diff(np.diff(S))/dz**2
        S += vert_diffusivity*dSdt*dt
        
        dUdt = np.zeros(len(U))
        dUdt[1:-1] = np.diff(np.diff(U))/dz**2
        U += vert_diffusivity*dUdt*dt
        
        dVdt = np.zeros(len(V))
        dVdt[1:-1] = np.diff(np.diff(V))/dz**2
        V += vert_diffusivity*dSdt*dt
        
        if tracers is not None:
            dtdt = np.zeros(shape=tracer.shape)
            dtdt[:,1:-1] = np.diff(np.diff(tracer,axis=0),axis=0)/dz**2
            tracer += vert_diffusivity*dtdt*dt
        
        
    return_variables = (T,S,U,V,)
    if tracer is not None:
        return_variables += (tracer,)
    if return_MLD:
        return_variables += (get_mld(T,S,U,V,z),)
    return return_variables    
   
    





# function to get density (rho)    
def get_rho(T,S,T0,S0,rho0,alpha,beta):
    return rho0 + alpha*(T-T0) + beta*(S-S0)

def get_mld(T,S,U,V,z):
    return z[get_mld_ind(T,S,U,V,z)]

# function to get mld index
def get_mld_ind(T,S,U,V,z):
    return np.where(np.logical_or(np.abs(T-T[0])>0.02,np.logical_or(np.abs(U-U[0])>0.01,np.abs(V-V[0])>0.01)))[0][0]

# function to get mld (simply calls get_mld_ind function)    
def get_mld(T,S,U,V,z):
    return z[get_mld_ind(T,S,U,V,z)]

# function to calculate mld from a density threshold using density at 10 m
def calc_mld_ind(T,S,z,T0,S0,rho0,alpha,beta):
    rho = get_rho(T,S,T0,S0,rho0,alpha,beta)
    z_10_ind = np.argmin(np.abs(z-10))
    rho_10 = rho[z_10_ind]
    if np.any(rho[z_10_ind::]-rho_10>0.03):
        return np.where(rho[z_10_ind::]-rho_10>0.03)[0][0]+z_10_ind
    else:
        return None
    
# function to get calculated mld (simply calls calc_mld_ind function)
def calc_mld(T,S,z,T0,S0,rho0,alpha,beta):
    ind = calc_mld_ind(T,S,z,T0,S0,rho0,alpha,beta)
    return np.nan if ind is None else z[ind]

# function to rotate zonal (U) and meridional (V) velocities by a given angle (in radians)    
def rotate(angle,U,V):
    temp = (U+V*1j)*np.exp(angle*1j)
    return np.real(temp),np.imag(temp)

# function to ensure static stability in the water column
def static_stability(T,S,U,V,z,T0,S0,rho0,alpha,beta,tracer=None,verbose=False):
    mld_ind = get_mld_ind(T,S,U,V,z)
    rho = get_rho(T,S,T0,S0,rho0,alpha,beta)
    if ~np.any(np.diff(rho)<0):
        if verbose:
            print('No need to deepen ML to ensure static stability')
        return T,S,U,V,tracer
    T = T.copy(); S = S.copy(); U = U.copy(); V = V.copy();
    if tracer is not None:
        tracer = tracer.copy()
    
    T[0:mld_ind] = np.mean(T[0:mld_ind]); S[0:mld_ind] = np.mean(S[0:mld_ind]); 
    U[0:mld_ind] = np.mean(U[0:mld_ind]); V[0:mld_ind] = np.mean(V[0:mld_ind]);
    if tracer is not None:
        tracer[:,0:mld_ind] = np.mean(tracer[:,0:mld_ind],axis=1);
    rho = get_rho(T,S,T0,S0,rho0,alpha,beta)
    while np.any(np.diff(rho)<0):
        if verbose:
            print('Deepening ML to ensure static stability...')
        mld_ind += 1
        T[0:mld_ind] = np.mean(T[0:mld_ind]); S[0:mld_ind] = np.mean(S[0:mld_ind]); 
        U[0:mld_ind] = np.mean(U[0:mld_ind]); V[0:mld_ind] = np.mean(V[0:mld_ind]);
        if tracer is not None:
            tracer[0:mld_ind] = np.mean(tracer[0:mld_ind]);
        rho = get_rho(T,S,T0,S0,rho0,alpha,beta)
    return T,S,U,V,tracer

# function to ensure mixed layer stability in the water column    
def mixed_layer_stability(T,S,U,V,z,T0,S0,rho0,alpha,beta,tracer=None,verbose=False):
    T = T.copy(); S = S.copy(); U = U.copy(); V = V.copy();
    if tracer is not None:
        tracer = tracer.copy()
    
    rho = get_rho(T,S,T0,S0,rho0,alpha,beta)
    h = get_mld(T,S,U,V,z)
    Rb = calculate_Rb(rho,h,U,V,z,rho0)
    if verbose:
        print('Rb: %.2f'%Rb)
    while Rb<0.65:
        h = get_mld(T,S,U,V,z)
        j = np.where(z>h)[0][0]
        U[0:j+1] = np.mean(U[0:j+1]); V[0:j+1] = np.mean(V[0:j+1]); 
        T[0:j+1] = np.mean(T[0:j+1]); S[0:j+1] = np.mean(S[0:j+1]);
        if tracer is not None:
            tracer[:,0:j+1] = np.mean(tracer[:,0:j+1],axis=1);
        rho = get_rho(T,S,T0,S0,rho0,alpha,beta)
        h = get_mld(T,S,U,V,z)
        Rb = calculate_Rb(rho,h,U,V,z,rho0)
        if verbose:
            print('Rb: %.2f'%Rb)
    return T,S,U,V,tracer

# function to ensure shear stability in the water column    
def shear_stability(T,S,U,V,z,T0,S0,rho0,alpha,beta,tracer=None,verbose=False,N_THRESHOLD=10000):
    T = T.copy(); S = S.copy(); U = U.copy(); V = V.copy();
    if tracer is not None:
        tracer = tracer.copy()
    
    Rgprime = 0.3
    rho = get_rho(T,S,T0,S0,rho0,alpha,beta)
    h = get_mld(T,S,U,V,z)
    Rg = calculate_Rg(rho,h,U,V,z,rho0)
    ind = np.nanargmin(Rg)
    if verbose:
        print('Rg: %.2f'%Rg[ind])
    h = get_mld(T,S,U,V,z)
    counter = 0
    while np.round(Rg[ind],4)<0.25:
        U[ind],U[ind+1] = U[ind]-(1-Rg[ind]/Rgprime)*(U[ind]-U[ind+1])/2,U[ind+1]+(1-Rg[ind]/Rgprime)*(U[ind]-U[ind+1])/2
        V[ind],V[ind+1] = V[ind]-(1-Rg[ind]/Rgprime)*(V[ind]-V[ind+1])/2,V[ind+1]+(1-Rg[ind]/Rgprime)*(V[ind]-V[ind+1])/2
        T[ind],T[ind+1] = T[ind]-(1-Rg[ind]/Rgprime)*(T[ind]-T[ind+1])/2,T[ind+1]+(1-Rg[ind]/Rgprime)*(T[ind]-T[ind+1])/2
        S[ind],S[ind+1] = S[ind]-(1-Rg[ind]/Rgprime)*(S[ind]-S[ind+1])/2,S[ind+1]+(1-Rg[ind]/Rgprime)*(S[ind]-S[ind+1])/2
        if tracer is not None:
            tracer[:,ind],tracer[:,ind+1] = tracer[:,ind]-(1-Rg[ind]/Rgprime)*(tracer[:,ind]-tracer[:,ind+1])/2,tracer[:,ind+1]+(1-Rg[ind]/Rgprime)*(tracer[:,ind]-tracer[:,ind+1])/2
        rho = get_rho(T,S,T0,S0,rho0,alpha,beta)
        Rg = calculate_Rg(rho,h,U,V,z,rho0)
        ind = np.nanargmin(Rg)
        if verbose:
            print('Rg: %.2f'%Rg[ind])
        T,S,U,V,tracer = mixed_layer_stability(T,S,U,V,z,T0,S0,rho0,alpha,beta,tracer=tracer,verbose=verbose) # remove any mixed layer instabilities
        counter += 1
        if counter>N_THRESHOLD:
            warnings.warn('Could not resolve shear instability (n>%d)'%N_THRESHOLD)
            break;
    return T,S,U,V,tracer

# calculate bulk Richardson number    
def calculate_Rb(rho,h,U,V,z,rho0):
    g = 9.81 # gravitational acceleration (m s^-2)
    j = np.where(z>h)[0][0]
    dU = U[j]-np.mean(U[0:j]); dV = V[j]-np.mean(V[0:j]); drho = rho[j]-np.mean(rho[0:j]);
    return g*drho*h/rho0/(dU**2+dV**2) # calculate balanced Richardson number

# calculate gradient Richardson number    
def calculate_Rg(rho,h,U,V,z,rho0):
    g = 9.81 # gravitational acceleration (m s^-2)
    Rg = g*np.diff(rho)*np.diff(z)/rho0/(np.diff(U)**2+np.diff(V)**2)
    Rg[z[0:-1]+np.diff(z)/2<h] = .25 # disregard these
    return Rg