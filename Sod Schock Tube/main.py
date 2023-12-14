def lax_friedrichs(U_inicial,nt,dt,nx):
    
    dx = Lx/(nx-1)
    
    U = U_inicial.copy()
    f = flux(U)
    
    for n in range(nt):
        
        Un = U.copy()
        fn = f.copy()

        Ux = 0.5*(Un[:,2:] + Un[:,:-2])
        fx = 0.5*(fn[:,2:] - fn[:,:-2])

        U[:,1:-1] = Ux - (dt/dx)*fx
        U[:,0] = U[:,1]
        U[:,-1] = U[:,-2]
        
        f = flux(U)
            
    return U,f
      
 #-------------------------   
   
def richtmyers(U_inicial,nt,dt,nx):
    
    dx = Lx/(nx-1)
    
    U = U_inicial.copy()
    f = flux(U)
    U1 = U.copy()
    U2 = U.copy()
    
    for n in range(nt):
        
        Un = U.copy()
        fn = f.copy()

        U1[:,1:-1] = 0.5*(Un[:,2:] + Un[:,1:-1]) - (dt/dx)*0.5*(fn[:,2:] - fn[:,1:-1])
        U2[:,1:-1] = 0.5*(Un[:,1:-1] + Un[:,:-2]) - (dt/dx)*0.5*(fn[:,1:-1] - fn[:,:-2])
        
        U[:,1:-1] = Un[:,1:-1] - (dt/dx)*(flux(U1)[:,1:-1]-flux(U2)[:,1:-1])
        U[:,0] = U[:,1]
        U[:,-1] = U[:,-2]
        
        f = flux(U)

    return U,f

    
#---------------------
    
def obtener_var(A):
    rho = A[0,:].copy()
    u = (A[1,:]/rho)
    E = A[2,:].copy()
    p = (gamma-1)*(E-0.5*rho*u**2)
    R = 287
    T = p/(rho*R)
    
    return rho,u,E,p,T

def flux(A):
    f = np.zeros((3,nx))
    rho,u,E,p,_ = obtener_var(A)
    
    f[0,:] = np.full(nx,u*rho)
    f[1,:] = np.full(nx,rho*u**2+p) 
    f[2,:] = np.full(nx,u*E+u*p)
    
    return f