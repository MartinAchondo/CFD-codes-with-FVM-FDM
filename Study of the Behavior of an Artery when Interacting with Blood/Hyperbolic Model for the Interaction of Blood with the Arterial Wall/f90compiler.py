
import os
import ctypes as ct



class f90_compiler():

    def __init__(self,name,lib,vars):
        root = os.path.dirname(os.path.realpath(__file__))
        path = os.path.join(root,lib)
        fortlib = ct.CDLL(path) 
        f = fortlib[name]
        f.argtypes=[*vars] 
        self.vars = vars
        self.f = f


    def leap_frog_func(self,u_out,dt,dx,gamma,Nt,Nx,f,u0):

        self.u = u_out

        L_args = [
            self.u.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.c_double(dt),
            ct.c_double(dx),
            ct.c_double(gamma),
            ct.c_int(Nt),
            ct.c_int(Nx),
            ct.c_int(f),
            u0.ctypes.data_as(ct.POINTER(ct.c_double))
        ]

        self.f(*L_args)

        return self.u

        
    def newmark_func(self,u_out,dt,dx,gamma,Nt,Nx,f,beta,theta,u0):

        self.uu = u_out

        L_args = [
            self.uu.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.c_double(dt),
            ct.c_double(dx),
            ct.c_double(gamma),
            ct.c_int(Nt),
            ct.c_int(Nx),
            ct.c_int(f),
            ct.c_double(beta),
            ct.c_double(theta),
            u0.ctypes.data_as(ct.POINTER(ct.c_double))
        ]

        self.f(*L_args)

        return self.uu




def get_types():
    types = {
        'pointer_d': ct.POINTER(ct.c_double),
        'pointer_i': ct.POINTER(ct.c_int),
        'double': ct.c_double,
        'int': ct.c_int
    }
    return types
