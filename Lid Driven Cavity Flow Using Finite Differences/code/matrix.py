import numpy as np
import scipy as sp
from scipy import sparse

class Matrix():

    def __init__(self):
        pass

    def generate_Matrix(self,nx,ny,dx,dy):
        nxn = nx-2
        nyn = ny-2
        N2 = nxn*nyn
        
        rows = np.array([])
        cols = np.array([])
        vals = np.array([])
        
        dx2 = dx**2
        dy2 = dy**2
        
        flag = True
        for k in range(nyn):  #num de sector
            for i in range(nxn):  #fila en sect
                
                ii = i + k*nxn
                #primer sector:
                if k==0:
                    if i==0:
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii, 1/dx2 + 1/dy2)
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii+1, -1/dx2)
                    elif i == nxn-1:
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii, 1/dx2 + 1/dy2)
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii-1, -1/dx2)
                    else:
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii, 2/dx2 + 1/dy2)
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii-1, -1/dx2)
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii+1, -1/dx2)
                    rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii+nxn, -1/dy2)
                    
                elif k==nyn-1:
                    if i==0:
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii, 1/dx2 + 1/dy2)
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii+1, -1/dx2)
                    elif i == nxn-1:
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii, 1/dx2 + 1/dy2)
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii-1, -1/dx2)
                    else:
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii, 2/dx2 + 1/dy2)
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii-1, -1/dx2)
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii+1, -1/dx2)
                    rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii-nxn, -1/dy2)
                        
                else:
                    if i==0:
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii, 1/dx2 + 2/dy2)
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii+1, -1/dx2)
                    elif i == nxn-1:
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii, 1/dx2 + 2/dy2)
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii-1, -1/dx2)
                    else:
                        if flag:
                            if k==int(nyn/2):
                                if i == int(nxn/2):
                                    rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii, 5/dx2 + 5/dy2)
                                    flag = False
                                else:
                                    rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii, 2/dx2 + 2/dy2)
                            else:
                                rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii, 2/dx2 + 2/dy2) 
                        else:
                            rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii, 2/dx2 + 2/dy2)
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii-1, -1/dx2)
                        rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii+1, -1/dx2)
                    rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii+nxn, -1/dy2)
                    rows,cols,vals = self.coordinates(rows,cols,vals,ii,ii-nxn, -1/dy2)
        
        A = sparse.csr_matrix((vals,(rows,cols)),shape = (N2,N2))
        return A

    def coordinates(self,rows,cols,vals,r,c,v):
        rows = np.append(rows,r)
        cols = np.append(cols,c)
        vals = np.append(vals,v)
        return rows,cols,vals