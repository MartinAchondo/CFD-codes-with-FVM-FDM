
from f90compiler import f90_compiler
from f90compiler import get_types
import numpy as np
import matplotlib.pyplot as plt

class Theta_Method():

    def __init__(self):

        types = get_types()
        function_name = 'theta_method'
        file_name = 'sub_theta.so'
        arg_types = [
            types['double'],
            types['pointer_d'],
            types['double'],
            types['double'],
            types['double'],
            types['pointer_d'],
            types['pointer_d'],
            types['pointer_d'],
            types['int'],
            types['double'],
            types['double'],
            types['double']
        ]  

        self.theta_class = f90_compiler(function_name, file_name, arg_types)


    def set_params(self):
        L = 5*10**-2
        R0 = 5*10**-3
        rho_w = 10**3
        H = 3*10**-4
        E = 9*10**5

        self.alpha = E/(rho_w*R0**2)
        self.beta = 0
        self.gamma = 1.0/(rho_w*H)
    
        self.x = L/2


    def solve_system(self,theta,tf,h):

        y0 = np.array([[0],[0]])
        t = 0.0
        
        N = int(tf/h)

        y1 = np.zeros((N,1), order='F')
        y2 = np.zeros((N,1), order='F')
        t = np.zeros((N,1), order='F')

        y1,y2,t = self.theta_class.theta_func(theta,y0,h,self.x,tf,t,y1,y2,N,self.alpha,self.beta,self.gamma)

        return y1,y2,t

    
    def eigen_values(self):
        lambda_1 = (-self.beta + np.sqrt(self.beta**2-4*self.alpha + 0j))/2
        lambda_2 = (-self.beta - np.sqrt(self.beta**2-4*self.alpha + 0j))/2
        return lambda_1,lambda_2

if __name__=='__main__':
    
    th = Theta_Method()

    th.set_params()
    th.beta = np.sqrt(th.alpha)
    #th.beta = th.alpha

    theta = 0.5
    tf = 2.5*10**-3
    h = 10**-4

    y1,y2,t = th.solve_system(theta,tf,h)

    print(th.eigen_values())
    print(th.alpha,th.beta )
    print(th.x)

    print(th.gamma)
