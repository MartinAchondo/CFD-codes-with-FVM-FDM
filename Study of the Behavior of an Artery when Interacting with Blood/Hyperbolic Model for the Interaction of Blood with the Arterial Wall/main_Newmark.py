from f90compiler import f90_compiler
from f90compiler import get_types
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import cm
types = get_types()


class Newmark():

    def __init__(self,problem=1):
        
        types = get_types()
        function_name = 'newmark'
        file_name = 'newmark.so'
        arg_types = [
            types['pointer_d'],
            types['double'],
            types['double'],
            types['double'],
            types['int'],
            types['int'],
            types['int'],
            types['double'],
            types['double'],
            types['pointer_d']
        ]  

        self.newmark_class = f90_compiler(function_name, file_name, arg_types)

        self.prob = problem
        self.set_params()
        


    def set_params(self):
        rho_w = 10**3
        H = 3*10**-4
        sigma = 1
        self.L = 5*10**-2
        self.T = 1
        self.gamma = np.sqrt(sigma/(rho_w*H))
        if self.prob == 1:
            self.gamma = 0.9
        self.f_source = self.prob


    def calculate_uk_NW(self,k):

        if self.prob==1:
            tf = 1
            xmax = 1
            dx,dt = self.get_mesh_1(k)
        elif self.prob==2:
            tf = self.T
            xmax = self.L
            dx,dt = self.get_mesh_2(k)

        Nt = int(tf/dt) +1
        Nx = int(xmax/dx) +1
        x = np.linspace(0,xmax,Nx)

        u_0 = np.zeros((2,Nx), order='F')
        if self.prob==1:
            u_0[0,:] = np.sin(np.pi*x)
            u_0[1,:] = -np.sin(np.pi*x)

        beta = 0.25
        theta = 0.5

        u_out = np.zeros((Nt,Nx), order='F')
        u_out = self.newmark_class.newmark_func(u_out,dt,dx,self.gamma,Nt,Nx,self.f_source,beta,theta,u_0)


        if self.prob==2:
            return u_out,x,np.linspace(0,self.T,Nt)

        t = np.linspace(0.1,tf,10)
        n_p = int(tf/(10*dt))
        u_j = np.zeros((10,Nx))

        for j in range(1,11):
            u_j[j-1,:] = u_out[j*n_p,:]

        return u_j,x,t

    def analytic(self,x,t):
        an = np.exp(-t)*np.sin(np.pi*x)
        return an


    def get_mesh_1(self,k):
        dx = 1/(2**k*10)
        dt = 1/(2**k*10)
        return dx,dt

    
    def get_mesh_2(self,k):
        dx = self.L/10
        dt = self.T/(100*(k)**2)
        return dx,dt


    def get_error(self,u_an,u_ap):
        return np.max(np.abs(u_an-u_ap))


    def p_convergence(self,k):
        ujk,xk,t = self.calculate_uk_NW(k)
        uj0,x0,t2 = self.calculate_uk_NW(0)

        ejk = np.zeros(len(t))
        ej0 = np.zeros(len(t))

        for j in range(len(t)):
            u_an = self.analytic(xk,t[j])
            ejk[j] = self.get_error(u_an,ujk[j,:])

        for j in range(len(t2)):
            u_an = self.analytic(x0,t2[j])
            ej0[j] = self.get_error(u_an,uj0[j,:])

        pj = np.log(ej0/ejk)/np.log(2**k)

        return pj


def plot_p1(k):
    nm = Newmark(problem=1)

    u_out,x,t = nm.calculate_uk_NW(k)

    for j in range(10):
        plt.plot(x,u_out[j,:]);
    plt.show();

    p = nm.p_convergence(k)
    print(p)


def plot_p2(k):
    nm = Newmark(problem=2)

    u_out,x,t = nm.calculate_uk_NW(k)

    # for j in range(len(t)):
    #     plt.plot(x,u_out[j,:]);
    # plt.show();

    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})

    X,T = np.meshgrid(x,t)
    surf = ax.plot_surface(X,T,u_out, cmap = cm.coolwarm)
    ax.view_init(20,20)
    plt.show();


def plot_an():
    newm = Newmark()
    x = np.linspace(0,1,100)
    t = np.linspace(0.1,1,10)
    for j in range(10):
        u_an = newm.analytic(x,t[j])
        plt.plot(u_an)
    plt.show()


if __name__=='__main__':

    plot_p1(k=1)
    # plot_an()