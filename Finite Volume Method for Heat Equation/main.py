import numpy
from matplotlib import pyplot

class cell():

    def __init__(self,nx,ny):
        self.ctr    = numpy.zeros((ny,nx,2))
        self.vSW = numpy.zeros((ny,nx,2))
        self.vSE = numpy.zeros((ny,nx,2))
        self.vNE = numpy.zeros((ny,nx,2))
        self.vNW = numpy.zeros((ny,nx,2))
        self.nN  = numpy.zeros((ny,nx,2))
        self.nS  = numpy.zeros((ny,nx,2))
        self.nE  = numpy.zeros((ny,nx,2))
        self.nW  = numpy.zeros((ny,nx,2))
        self.area = numpy.zeros((ny,nx))

        self.T  = numpy.zeros((ny,nx))
        self.Tn = numpy.zeros((ny,nx))

def generate_mesh(nx, ny, Lx, Ly, theta, wedge):

    x = numpy.linspace(0, Lx, nx+1) # +1 porque nx es numero de celdas,
    X = numpy.zeros((ny+1,nx+1))    # no vertices
    Y = numpy.zeros((ny+1,nx+1))
    dy = Ly/ny
    dx = Lx/nx

    for j in range(ny+1):

        X[j,:] = x
        i_before = numpy.where(x<=wedge)[0]
        Y[j,i_before] = dy*j
        
        theta_j = theta*(Ly-j*dy)/Ly
        i_after = numpy.where(x>wedge)[0]
        Y[j,i_after] = dy*j + (dx*i_after-wedge)*numpy.tan(theta_j)

    
    vol_test = numpy.zeros((ny,nx,2))
    volumes = cell(nx,ny)
    volumes.vSW = numpy.dstack((X[:-1,:-1],Y[:-1,:-1]))
    volumes.vSE = numpy.dstack((X[:-1,1:],Y[:-1,1:]))
    volumes.vNE = numpy.dstack((X[1:,1:],Y[1:,1:]))
    volumes.vNW = numpy.dstack((X[1:,:-1],Y[1:,:-1]))

    volumes.ctr = ( volumes.vSW + volumes.vSE + volumes.vNE + volumes.vNW )/4

    volumes.area = compute_area(volumes.vSW,\
                                volumes.vSE,\
                                volumes.vNE,\
                                volumes.vNW)

    lS = numpy.linalg.norm(volumes.vSE-volumes.vSW, axis=2)
    volumes.nS  = numpy.dstack((volumes.vSE[:,:,1]-volumes.vSW[:,:,1], \
                              -(volumes.vSE[:,:,0]-volumes.vSW[:,:,0]))/lS)

    lE = numpy.linalg.norm(volumes.vSE-volumes.vNE, axis=2)
    volumes.nE  = numpy.dstack((volumes.vNE[:,:,1]-volumes.vSE[:,:,1], \
                              -(volumes.vNE[:,:,0]-volumes.vSE[:,:,0]))/lE)

    lN = numpy.linalg.norm(volumes.vNW-volumes.vNE, axis=2)
    volumes.nN  = numpy.dstack((volumes.vNW[:,:,1]-volumes.vNE[:,:,1], \
                              -(volumes.vNW[:,:,0]-volumes.vNE[:,:,0]))/lN)

    lW = numpy.linalg.norm(volumes.vNW-volumes.vSW, axis=2)
    volumes.nW  = numpy.dstack((volumes.vSW[:,:,1]-volumes.vNW[:,:,1], \
                              -(volumes.vSW[:,:,0]-volumes.vNW[:,:,0]))/lW)

    return volumes

def compute_area(A, B, C, D):

    p = A - C
    q = D - B

    pq = numpy.cross(p,q)
    
    normpq = numpy.abs(pq)

    area = 0.5*normpq
    
    return area

def initial_conditions(volumes, T, Tb):

    ny, nx = numpy.shape(volumes.T)

    volumes.T[:,:]  = T
    volumes.Tn[:,:] = T


def node_g(cSW, cSE, cNE, cNW):

    area = compute_area(cSW[1], cSE[1], cNE[1], cNW[1])

    g = -1/(2*area) * ( (cSW[0]+cSE[0])*(cSE[1][:,:,0]-cSW[1][:,:,0]) \
                      + (cSE[0]+cNE[0])*(cNE[1][:,:,0]-cSE[1][:,:,0]) \
                      + (cNE[0]+cNW[0])*(cNW[1][:,:,0]-cNE[1][:,:,0]) \
                      + (cNW[0]+cSW[0])*(cSW[1][:,:,0]-cNW[1][:,:,0]) )
    return g


def node_f(cSW, cSE, cNE, cNW):

    area = compute_area(cSW[1], cSE[1], cNE[1], cNW[1])
    #print(numpy.shape(cSW[1]),numpy.shape(cNE[1]),numpy.shape(cSE[1]),numpy.shape(cNW[1]))
    #print(numpy.shape(cSW[0]),numpy.shape(cNE[0]),numpy.shape(cSE[0]),numpy.shape(cNW[0]))
    f = 1/(2*area) * ( (cSW[0]+cSE[0])*(cSE[1][:,:,1]-cSW[1][:,:,1]) \
                     + (cSE[0]+cNE[0])*(cNE[1][:,:,1]-cSE[1][:,:,1]) \
                     + (cNE[0]+cNW[0])*(cNW[1][:,:,1]-cNE[1][:,:,1]) \
                     + (cNW[0]+cSW[0])*(cSW[1][:,:,1]-cNW[1][:,:,1]) )
    return f

def computeF(cell, cN, cS, cE, cW, cSW, cSE, cNW, cNE):

#   order:  D---C
#           |   |
#           |   |
#           A---B

    fA = node_f(cSW, cS, cell, cW) 
    fB = node_f(cS, cSE, cE, cell) 
    fC = node_f(cell, cE, cNE, cN) 
    fD = node_f(cW, cell, cN, cNW) 

    flux_bot = (fA+fB)/2
    flux_rgt = (fB+fC)/2
    flux_top = (fC+fD)/2
    flux_lft = (fD+fA)/2

    return flux_bot, flux_rgt, flux_top, flux_lft


def computeG(cell, cN, cS, cE, cW, cSW, cSE, cNW, cNE):

#   order:  D---C
#           |   |
#           |   |
#           A---B

    gA = node_g(cSW, cS, cell, cW) 
    gB = node_g(cS, cSE, cE, cell) 
    gC = node_g(cell, cE, cNE, cN) 
    gD = node_g(cW, cell, cN, cNW) 

    flux_bot = (gA+gB)/2
    flux_rgt = (gB+gC)/2
    flux_top = (gC+gD)/2
    flux_lft = (gD+gA)/2

    return flux_bot, flux_rgt, flux_top, flux_lft


def euler_step(mesh, dt, nt, nx, ny, alpha, qbot, qtop, qrgt, qlft):

    f_bot = numpy.zeros((ny,nx))
    f_rgt = numpy.zeros((ny,nx))
    f_top = numpy.zeros((ny,nx))
    f_lft = numpy.zeros((ny,nx))
    g_bot = numpy.zeros((ny,nx))
    g_rgt = numpy.zeros((ny,nx))
    g_top = numpy.zeros((ny,nx))
    g_lft = numpy.zeros((ny,nx))

    for t in range(nt):

        mesh.Tn[:,:] = mesh.T[:,:]

        cell= [mesh.Tn[1:-1,1:-1],mesh.ctr[1:-1,1:-1]]
        cN  = [mesh.Tn[2:,1:-1],mesh.ctr[2:,1:-1]]
        cS  = [mesh.Tn[:-2,1:-1],mesh.ctr[:-2,1:-1]]
        cE  = [mesh.Tn[1:-1,2:],mesh.ctr[1:-1,2:]]
        cW  = [mesh.Tn[1:-1,:-2],mesh.ctr[1:-1,:-2]]
        cSW = [mesh.Tn[:-2,:-2],mesh.ctr[:-2,:-2]]
        cSE = [mesh.Tn[:-2,2:],mesh.ctr[:-2,2:]]
        cNE = [mesh.Tn[2:,2:],mesh.ctr[2:,2:]]
        cNW = [mesh.Tn[2:,:-2],mesh.ctr[2:,:-2]]

#       Caso general
        f_bot[1:-1,1:-1], f_rgt[1:-1,1:-1], f_top[1:-1,1:-1], f_lft[1:-1,1:-1] = \
                            computeF(cell, cN, cS, cE, cW, cSW, cSE, cNW, cNE)
        g_bot[1:-1,1:-1], g_rgt[1:-1,1:-1], g_top[1:-1,1:-1], g_lft[1:-1,1:-1] = \
                            computeG(cell, cN, cS, cE, cW, cSW, cSE, cNW, cNE)

#       Casos particulares
#       Borde inferior 
        cell= [mesh.Tn[[0],1:-1],mesh.ctr[[0],1:-1]]
        cN  = [mesh.Tn[[1],1:-1],mesh.ctr[[1],1:-1]]
        cE  = [mesh.Tn[[0],2:],mesh.ctr[[0],2:]]
        cW  = [mesh.Tn[[0],:-2],mesh.ctr[[0],:-2]]
        cNE = [mesh.Tn[[1],2:],mesh.ctr[[1],2:]]
        cNW = [mesh.Tn[[1],:-2],mesh.ctr[[1],:-2]]

        fC = node_f(cell, cE, cNE, cN) 
        gC = node_g(cell, cE, cNE, cN) 
        fD = node_f(cW,cell, cN, cNW) 
        gD = node_g(cW,cell, cN, cNW) 

        f_bot[0,1:-1] = qbot*mesh.nS[0,1:-1,0]
        f_lft[0,1:-1] = (fD+qbot*mesh.nS[0,1:-1,0])/2
        f_rgt[0,1:-1] = (fC+qbot*mesh.nS[0,1:-1,0])/2
        f_top[0,1:-1] = (fC+fD)/2

        g_bot[0,1:-1] = qbot*mesh.nS[0,1:-1,1]
        g_lft[0,1:-1] = (gD+qbot*mesh.nS[0,1:-1,1])/2
        g_rgt[0,1:-1] = (gC+qbot*mesh.nS[0,1:-1,1])/2
        g_top[0,1:-1] = (gC+gD)/2

#       Borde derecho
        cell= [mesh.Tn[1:-1,[-1]],mesh.ctr[1:-1,[-1]]]
        cN  = [mesh.Tn[2:,[-1]],mesh.ctr[2:,[-1]]]
        cS  = [mesh.Tn[:-2,[-1]],mesh.ctr[:-2,[-1]]]
        cW  = [mesh.Tn[1:-1,[-2]],mesh.ctr[1:-1,[-2]]]
        cSW = [mesh.Tn[:-2,[-2]],mesh.ctr[:-2,[-2]]]
        cNW = [mesh.Tn[2:,[-2]],mesh.ctr[2:,[-2]]]

        fA = node_f(cSW,cS,cell,cW) 
        gA = node_g(cSW,cS,cell,cW) 
        fD = node_f(cW,cell, cN, cNW) 
        gD = node_g(cW,cell, cN, cNW) 

       # print(numpy.shape(fA), numpy.shape(mesh.nE[1:-1,-1,0]))
       # print(numpy.shape(qrgt),numpy.shape(mesh.nE[1:-1,-1,0]))
        f_bot[1:-1,-1] = (fA[:,0]+qrgt*mesh.nE[1:-1,-1,0])/2
        f_rgt[1:-1,-1] = qrgt*mesh.nE[1:-1,-1,0]
        f_lft[1:-1,-1] = (fA[:,0]+fD[:,0])/2
        f_top[1:-1,-1] = (fD[:,0]+qrgt*mesh.nE[1:-1,-1,0])/2

        g_bot[1:-1,-1] = (gA[:,0]+qrgt*mesh.nE[1:-1,-1,1])/2
        g_rgt[1:-1,-1] = qrgt*mesh.nE[1:-1,-1,1]
        g_lft[1:-1,-1] = (gA[:,0]+gD[:,0])/2
        g_top[1:-1,-1] = (gD[:,0]+qrgt*mesh.nE[1:-1,-1,1])/2

#       Borde superior
        cell= [mesh.Tn[[-1],1:-1],mesh.ctr[[-1],1:-1]]
        cS  = [mesh.Tn[[-2],1:-1],mesh.ctr[[-2],1:-1]]
        cE  = [mesh.Tn[[-1],2:],mesh.ctr[[-1],2:]]
        cW  = [mesh.Tn[[-1],:-2],mesh.ctr[[-1],:-2]]
        cSW = [mesh.Tn[[-2],:-2],mesh.ctr[[-2],:-2]]
        cSE = [mesh.Tn[[-2],2:],mesh.ctr[[-2],2:]]

        fA = node_f(cSW,cS,cell,cW) 
        gA = node_g(cSW,cS,cell,cW) 
        fB = node_f(cS,cSE,cE,cell) 
        gB = node_g(cS,cSE,cE,cell) 

        f_bot[-1,1:-1] = (fA+fB)/2
        f_rgt[-1,1:-1] = (fB+qtop*mesh.nN[-1,1:-1,0])/2
        f_lft[-1,1:-1] = (fA+qtop*mesh.nN[-1,1:-1,0])/2
        f_top[-1,1:-1] = qtop*mesh.nN[-1,1:-1,0]

        g_bot[-1,1:-1] = (gA+gB)/2
        g_rgt[-1,1:-1] = (gB+qtop*mesh.nN[-1,1:-1,1])/2
        g_lft[-1,1:-1] = (gA+qtop*mesh.nN[-1,1:-1,1])/2
        g_top[-1,1:-1] = qtop*mesh.nN[-1,1:-1,1]

#       Borde izquierdo
        cell= [mesh.Tn[1:-1,[0]],mesh.ctr[1:-1,[0]]]
        cN  = [mesh.Tn[2:,[0]],mesh.ctr[2:,[0]]]
        cS  = [mesh.Tn[:-2,[0]],mesh.ctr[:-2,[0]]]
        cE  = [mesh.Tn[1:-1,[1]],mesh.ctr[1:-1,[1]]]
        cSE = [mesh.Tn[:-2,[1]],mesh.ctr[:-2,[1]]]
        cNE = [mesh.Tn[2:,[1]],mesh.ctr[2:,[1]]]

        fC = node_f(cell, cE, cNE, cN) 
        gC = node_g(cell, cE, cNE, cN) 
        fB = node_f(cS,cSE,cE,cell) 
        gB = node_g(cS,cSE,cE,cell) 

        f_bot[1:-1,0] = (fB[:,0]+qlft*mesh.nW[1:-1,0,0])/2
        f_rgt[1:-1,0] = (fB[:,0]+fC[:,0])/2
        f_top[1:-1,0] = (fC[:,0]+qlft*mesh.nW[1:-1,0,0])/2
        f_lft[1:-1,0] = qlft*mesh.nW[1:-1,0,0]

        g_bot[1:-1,0] = (gB[:,0]+qlft*mesh.nW[1:-1,0,1])/2
        g_rgt[1:-1,0] = (gB[:,0]+fC[:,0])/2
        g_top[1:-1,0] = (gC[:,0]+qlft*mesh.nW[1:-1,0,1])/2
        g_lft[1:-1,0] = qlft*mesh.nW[1:-1,0,1]

#       Esquina inferior izquierda
        # cell= [mesh.Tn[0,0],mesh.ctr[0,0]]
        # cN  = [mesh.Tn[1,0],mesh.ctr[1,0]]
        # cE  = [mesh.Tn[0,1],mesh.ctr[0,1]]
        # cNE = [mesh.Tn[1,1],mesh.ctr[1,1]]

        # fC = node_f(cell, cE, cNE, cN) 
        # gC = node_g(cell, cE, cNE, cN) 

        # f_bot[0,0] = qbot*mesh.nS[0,0,0]
        # f_lft[0,0] = qlft*mesh.nW[0,0,0]
        # f_rgt[0,0] = (fC+qbot*mesh.nS[0,0,0])/2
        # f_top[0,0] = (fC+qlft*mesh.nW[0,0,0])/2

        # g_bot[0,0] = qbot*mesh.nS[0,0,1]
        # g_lft[0,0] = qlft*mesh.nW[0,0,1]
        # g_rgt[0,0] = (gC+qbot*mesh.nS[0,0,1])/2
        # g_top[0,0] = (gC+qlft*mesh.nW[0,0,1])/2



        # for j in range(ny):
        #     for i in range(nx):


        #         cell = mesh[j,i]

        #         if i==0 and j==0:  # esquina inferior izquierda
        #             cN  = mesh[j+1,i]
        #             cE  = mesh[j,i+1]
        #             cNE = mesh[j+1,i+1]
        #             fC = node_f(cell, cE, cNE, cN) 
        #             gC = node_g(cell, cE, cNE, cN) 

        #             f_bot = qbot*cell.nS[0]
        #             f_lft = qlft*cell.nW[0]
        #             f_rgt = (fC+qbot*cell.nS[0])/2
        #             f_top = (fC+qlft*cell.nW[0])/2

        #             g_bot = qbot*cell.nS[1]
        #             g_lft = qlft*cell.nW[1]
        #             g_rgt = (gC+qbot*cell.nS[1])/2
        #             g_top = (gC+qlft*cell.nW[1])/2

        #         elif i==nx-1 and j==0: # esquina inferior derecha
        #             cN  = mesh[j+1,i]
        #             cW  = mesh[j,i-1]
        #             cNW = mesh[j+1,i-1]
        #             fD = node_f(cW,cell, cN, cNW) 
        #             gD = node_g(cW,cell, cN, cNW) 

        #             f_bot = qbot*cell.nS[0]
        #             f_rgt = qrgt*cell.nE[0]
        #             f_lft = (fD+qbot*cell.nS[0])/2
        #             f_top = (fD+qrgt*cell.nE[0])/2

        #             g_bot = qbot*cell.nS[1]
        #             g_rgt = qrgt*cell.nE[1]
        #             g_lft = (gD+qbot*cell.nS[1])/2
        #             g_top = (gD+qrgt*cell.nE[1])/2


        #         elif i==0 and j==ny-1: # esquina superior izquierda
        #             cS  = mesh[j-1,i]
        #             cE  = mesh[j,i+1]
        #             cSE = mesh[j-1,i+1]
        #             fB = node_f(cS,cSE,cE,cell) 
        #             gB = node_g(cS,cSE,cE,cell) 

        #             f_bot = (fB+qlft*cell.nW[0])/2
        #             f_rgt = (fB+qtop*cell.nN[0])/2
        #             f_lft = qlft*cell.nW[0]
        #             f_top = qtop*cell.nN[0]

        #             g_bot = (gB+qlft*cell.nW[1])/2
        #             g_rgt = (gB+qtop*cell.nN[1])/2
        #             g_lft = qlft*cell.nW[1]
        #             g_top = qtop*cell.nN[1]

        #         elif i==nx-1 and j==ny-1: # esquina superior derecha
        #             cS  = mesh[j-1,i]
        #             cW  = mesh[j,i-1]
        #             cSW = mesh[j-1,i-1]
        #             fA = node_f(cSW,cS,cell,cW) 
        #             gA = node_g(cSW,cS,cell,cW) 

        #             f_bot = (fA+qrgt*cell.nE[0])/2
        #             f_rgt = qrgt*cell.nE[0]
        #             f_lft = (fA+qtop*cell.nN[0])/2
        #             f_top = qtop*cell.nN[0]

        #             g_bot = (gA+qrgt*cell.nE[1])/2
        #             g_rgt = qrgt*cell.nE[1]
        #             g_lft = (gA+qtop*cell.nN[1])/2
        #             g_top = qtop*cell.nN[1]

                                      
    mesh.T[:,:] = mesh.Tn[:,:] + alpha*dt/mesh.area[:,:] * \
                  ( f_bot * (mesh.vSE[:,:,1]-mesh.vSW[:,:,1]) \
                  - g_bot * (mesh.vSE[:,:,0]-mesh.vSW[:,:,0]) \
                  + f_rgt * (mesh.vNE[:,:,1]-mesh.vSE[:,:,1]) \
                  - g_rgt * (mesh.vNE[:,:,0]-mesh.vSE[:,:,0]) \
                  + f_top * (mesh.vNW[:,:,1]-mesh.vNE[:,:,1]) \
                  - g_top * (mesh.vNW[:,:,0]-mesh.vNE[:,:,0]) \
                  + f_lft * (mesh.vSW[:,:,1]-mesh.vNW[:,:,1]) \
                  - g_lft * (mesh.vSW[:,:,0]-mesh.vNW[:,:,0]) )
                   

alpha = 0.1
T0 = 20.
Tb = 300.

qbot = 10.
qrgt = 0.
qtop = -10.
qlft = 0.

Lx = 1.5
Ly = 1.
wedge = 0.5
theta = 15*numpy.pi/180.
nx = 15
ny = 10

dt = 0.001
nt = 500

sigma = alpha*dt/(Lx/nx)**2
print('CFL = %1.4f'%sigma)

mesh = generate_mesh(nx, ny, Lx, Ly, theta, wedge)

initial_conditions(mesh, T0, Tb)

euler_step(mesh, dt, nt, nx, ny, alpha, qbot, qtop, qrgt, qlft)

# center = numpy.zeros((ny,nx,2))
# T = numpy.zeros((ny,nx))
# for j in range(ny):
#     for i in range(nx):
#         center[j,i,:] = mesh.ctr[j,i]
#         T[j,i] = mesh.ctr[j,i].T

#pyplot.pcolor(center[:,:,0],center[:,:,1],T)
pyplot.contourf(mesh.ctr[:,:,0],mesh.ctr[:,:,1],mesh.T)
# a = center[:,nx/2,0]
# b= center[:,nx/2,1]
# c = T[:,nx/2]
#pyplot.plot(b,c)
pyplot.colorbar()
pyplot.show()
