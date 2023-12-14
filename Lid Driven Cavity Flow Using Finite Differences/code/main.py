import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import sparse
from scipy.sparse.linalg import spsolve
plt.rcParams["figure.figsize"] = (10,6)
import time
import sys
from matrix import Matrix


class NS_Solver(Matrix):

    def __init__(self,Lx,Ly,nx,ny,dt):
        self.Lx = Lx
        self.Ly = Ly
        self.nx = nx
        self.ny = ny
        self.dt = dt

    def get_params(self,u_0,mu,rho):
        self.u_0 = u_0
        self.mu = mu
        self.rho = rho



    def method(self,U_inicial,V_inicial,nt,A,u_0,flag):
        
        dx = self.Lx/(self.nx-1)
        dy = self.Ly/(self.ny-1)
        
        u = U_inicial.copy()
        v = V_inicial.copy()
        
        Res_u = np.zeros(nt)
        Res_v = np.zeros(nt)
        Res_div = np.zeros(nt)
        Steady = np.zeros((3,nt))
        
        u_st = np.zeros((self.ny,self.nx),dtype='float64')
        v_st = np.zeros((self.ny,self.nx), dtype='float64')
        P = np.zeros((self.ny,self.nx))
        
        for n in range(nt):
            
            un = u.copy()
            vn = v.copy()
            Pn = P.copy()

            u_1x = (un[1:-1,1:-1]-un[1:-1,:-2])/dx
            u_1y = (un[1:-1,1:-1]-un[:-2,1:-1])/dy
            u_2x = (un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,:-2])/(dx**2)
            u_2y = (un[2:,1:-1]-2*un[1:-1,1:-1]+un[:-2,1:-1])/(dy**2)
            
            v_1x = (vn[1:-1,1:-1]-vn[1:-1,:-2])/dx
            v_1y = (vn[1:-1,1:-1]-vn[:-2,1:-1])/dy
            v_2x = (vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,:-2])/(dx**2)
            v_2y = (vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[:-2,1:-1])/(dy**2)
            
            u_st[1:-1,1:-1] = un[1:-1,1:-1] + self.dt*(-un[1:-1,1:-1]*u_1x - vn[1:-1,1:-1]*u_1y) + self.dt*self.mu*(u_2x + u_2y)
            v_st[1:-1,1:-1] = vn[1:-1,1:-1] + self.dt*(-un[1:-1,1:-1]*v_1x - vn[1:-1,1:-1]*v_1y) + self.dt*self.mu*(v_2x + v_2y)

            u_st[0,:] = np.zeros(self.nx)
            v_st[0,:] = np.zeros(self.nx)
            u_st[-1,:] = np.full(self.nx,u_0)
            v_st[-1,:] = np.zeros(self.nx)
            
            u_st[:,0] = np.zeros(self.ny)
            v_st[:,0] = np.zeros(self.ny)
            u_st[:,-1] = np.zeros(self.ny)
            v_st[:,-1] = np.zeros(self.ny)
            

            u_st_x = (u_st[1:-1,2:]-u_st[1:-1,:-2])/(2*dx)
            v_st_y = (v_st[2:,1:-1]-v_st[:-2,1:-1])/(2*dy)
            div_ji = (self.rho/self.dt)*(u_st_x + v_st_y)
            
            b = -1*div_ji.flatten()
            
            P_1d = spsolve(A,b)
            
            P[1:-1,1:-1] = np.reshape(P_1d,(-1,self.nx-2))
            
            P[0,:] = P[1,:]
            P[-1,:] = P[-2,:]
            P[:,0] = P[:,1]
            P[:,-1] = P[:,-2]
            
            P[0,0] = (1/(1/dx+1/dy))* (P[0,1]/dx + P[1,0]/dy)
            P[-1,-1] = (1/(1/dx+1/dy))* (P[-1,-2]/dx + P[-2,-1]/dy)
   
            P_x = (P[1:-1,2:]-P[1:-1,:-2])/(2*dx)
            P_y = (P[2:,1:-1]-P[:-2,1:-1])/(2*dy)
            
            u[1:-1,1:-1] = u_st[1:-1,1:-1] - (self.dt/self.rho)*P_x
            v[1:-1,1:-1] = v_st[1:-1,1:-1] - (self.dt/self.rho)*P_y
            
            u[0,:] = np.zeros(self.nx)
            v[0,:] = np.zeros(self.nx)
            u[-1,:] = np.full(self.nx,u_0)
            v[-1,:] = np.zeros(self.nx)
            
            u[:,0] = np.zeros(self.ny)
            v[:,0] = np.zeros(self.ny)
            u[:,-1] = np.zeros(self.ny)
            v[:,-1] = np.zeros(self.ny)
                  
            if flag:
                u_t = (u[1:-1,1:-1]-un[1:-1,1:-1])/self.dt

                u_x = (u[1:-1,1:-1]-u[1:-1,:-2])/dx
                u_y = (u[1:-1,1:-1]-u[:-2,1:-1])/dy
                u_xx = (u[1:-1,2:]-2*u[1:-1,1:-1]+u[1:-1,:-2])/(dx**2)
                u_yy = (u[2:,1:-1]-2*u[1:-1,1:-1]+u[:-2,1:-1])/(dy**2)

                P_x = (P[1:-1,2:]-P[1:-1,:-2])/(2*dx)

                conv_u = u_t + u[1:-1,1:-1]*u_x + v[1:-1,1:-1]*u_y 
                diff_u = -P_x/self.rho + self.mu*(u_xx + u_yy)

                Res_u[n] = np.sqrt(1/((self.nx-2)*(self.ny-2))*np.sum((conv_u-diff_u)**2))

                v_t = (v[1:-1,1:-1]-vn[1:-1,1:-1])/self.dt

                v_x = (v[1:-1,1:-1]-v[1:-1,:-2])/dx
                v_y = (v[1:-1,1:-1]-v[:-2,1:-1])/dy
                v_xx = (v[1:-1,2:]-2*v[1:-1,1:-1]+v[1:-1,:-2])/(dx**2)
                v_yy = (v[2:,1:-1]-2*v[1:-1,1:-1]+v[:-2,1:-1])/(dy**2)

                P_y = (P[2:,1:-1]-P[:-2,1:-1])/(2*dy)

                conv_v = v_t + u[1:-1,1:-1]*v_x + v[1:-1,1:-1]*v_y 
                diff_v = -P_y/self.rho + self.mu*(v_xx + v_yy)

                Res_v[n] = np.sqrt(1/((self.nx-2)*(self.ny-2))*np.sum((conv_v-diff_v)**2))

                ux = (u[1:-1,2:]-u[1:-1,:-2])/(2*dx)
                vy = (v[2:,1:-1]-v[:-2,1:-1])/(2*dy)

                Res_div[n] = np.sqrt(1/((self.nx-2)*(self.ny-2))*np.sum((ux+vy)**2)) 

                Steady[0,n] = np.sqrt(1/((self.nx-2)*(self.ny-2))*np.sum((u_t)**2))
                Steady[1,n] = np.sqrt(1/((self.nx-2)*(self.ny-2))*np.sum((v_t)**2))
                Steady[2,n] = np.sqrt(1/((self.nx-2)*(self.ny-2))*np.sum(((P-Pn)/self.dt)**2))
                
            if np.max(u)>10**4:
                break
            
        return u,v,P,Res_u,Res_v,Res_div,Steady


    def Adim_numbers(self):
        Re = int(self.u_0*self.Lx/self.mu)
        print(f'Número de Reynolds: {Re}')

        C = self.u_0*self.dt/self.dx

        CFL_diff = self.mu*self.dt/self.dx**2 + self.mu*self.dt/self.dy**2

        print('Número de Courant: ',np.max(C))
        print('Número de Difusión: ',CFL_diff)

        Pe = C/CFL_diff

        print('Número de Peclet: ', np.max(Pe))

        return Re,C,CFL_diff,Pe


    def vorticity(self,u,v,dx,dy):
        v_x = (v[1:-1,2:]-v[1:-1,:-2])/(2*dx)
        u_y = (u[2:,1:-1]-u[:-2,1:-1])/(2*dy)

        w = np.zeros((ny,nx))
        w[1:-1,1:-1] = v_x - u_y 

        w[0,1:-1] = (v[0,2:]-v[0,1:-1])/(dx) - (u[1,1:-1]-u[0,1:-1])/(dy)
        w[-1,1:-1] = (v[-1,2:]-v[-1,1:-1])/(dx) - (u[-1,1:-1]-u[-2,1:-1])/(dy)
        w[1:-1,0] = (v[1:-1,1]-v[1:-1,0])/(dx) - (u[2:,0]-u[1:-1,0])/(dy)
        w[1:-1,-1] = (v[1:-1,-1]-v[1:-1,0])/(dx) - (u[2:,-1]-u[1:-1,-1])/(dy)

        return w


    def stream_lines(self,u,nx,ny,dy):
        psi = np.zeros((ny,nx))

        for j in range(ny-1):
            psi[j+1,:] = psi[j,:] + (dy/2)*(u[j,:]+u[j+1,:])

        return psi


if __name__=='__main__':
    Lx = 2                           
    Ly = 2                           
    rho = 1.                         
    mu = 0.1                         
    u_0 = 1                          
    
    dt = 0.001                       
    nx = 41                          
    ny = 41                          
          
    flag= True
    nt =10000
    U_inicial = np.zeros((ny,nx))
    U_inicial[-1,:] = np.full(nx,u_0)
    V_inicial = np.zeros((ny,nx))
    P = np.zeros((ny,nx))

    NS = NS_Solver(Lx,Ly,nx,ny,dt)
    A = NS.generate_Matrix(NS.nx,NS.ny,NS.dx,NS.dy)

    NS.get_params(u_0,mu,rho)
    Re,C,C_diff,Pe = NS.Adim_numbers()

    u,v,P,Res_u,Res_v,Res_div,Steady = NS.metodo(U_inicial,V_inicial,nt,A,u_0,flag)

    print(f'Estado después de {nt*dt} segundos')

    x = np.linspace(0,Lx,nx)
    y = np.linspace(0,Ly,ny)

    X,Y = np.meshgrid(x,y)

    plt.streamplot(X, Y, u, v, color='k', density=1)

    plt.contourf(X,Y,P/np.linalg.norm(P), cmap="rainbow")
    plt.colorbar()

    plt.title(f'Campo de Velocidades y Presión Re = {Re}');
