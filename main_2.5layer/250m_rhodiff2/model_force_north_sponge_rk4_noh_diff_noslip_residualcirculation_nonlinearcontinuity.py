#In this file we have grid functions (e.g. boundary conditions and taking derivatives) and model functions (the actual equations implemented)
#Realized because of the staggered grid that model_force_north_sponge_rk4_noh_diff_noslip was implementing another condition on top of no slip 
#Here we implement no slip as in https://extranet.gfdl.noaa.gov/~aja/papers/adcroft_marshall_TELLUS_1998.pdf

import numpy as np
import copy
import sys
import os

# Get parent directory then add it
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)


#use this function to import in everything that we defined in the notebook we execute
def import_initialvalues(dto, to, _u1o,_v1o,_h1o,_h2o,_u2o, _v2o,nyo,dxo,dyo,tau_mino,tau_maxo,H1o,H2o,nuo,ro,tau_xo,f0o,betao,rho1o,rho2o,rho3o,uyo,vyo,m1o=0,m2o=0,extento=0,start_to=0,vol_transfero=0,exponential_decayro=5): #the o after each variable name implies original (so from the file where we run)
#Here we make the mass input terms (m1 for top layer, m2 for bottom layer) optional by setting them as =0 by default (e.g. no mas input  by default). Extent is the number of grid cells that the mass is input over
    global dt,t,u1, v1, h1, u2, v2, h2, _u1, _v1,_h1,_h2,_u2, _v2,ny,dx,dy,tau_min,tau_max,H1,H2,nu,r,f0,beta,rho1,rho2,rho3,g1,g2,uy,vy,tau_x,m1,m2,extent, start_t,vol_transfer,exponential_decayr #make all these global so they can be used
    
    dt=dto
    t=to
    _u1=_u1o
    _v1=_v1o
    _h1=_h1o
    _h2=_h2o
    _u2=_u2o
    _v2=_v2o
    ny=nyo
    dx=dxo
    dy=dyo
    u1 = _u1[1:-1, 1:-1]               # (nx+1, ny)
    v1 = _v1[1:-1, 1:-1]               # (nx, ny+1)
    h1 = _h1[1:-1, 1:-1]               # (nx, ny)
    u2 = _u2[1:-1, 1:-1]               # (nx+1, ny)
    v2 = _v2[1:-1, 1:-1]               # (nx, ny+1)
    h2 = _h2[1:-1, 1:-1]               # (nx, ny)
    tau_min=tau_mino
    tau_max=tau_maxo
    H1=H1o
    H2=H2o
    nu=nuo
    tau_x=tau_xo
    r=ro
    f0=f0o
    beta=betao
    rho1=rho1o
    rho2=rho2o
    rho3=rho3o
    g1= 9.8*(rho2-rho1)/rho1 #reduced gravity 1
    g2= 9.8*(rho3-rho2)/rho1 #reduced gravity 1
    uy=uyo
    vy=vyo
    m1=m1o
    m2=m2o
    extent=extento
    start_t=start_to #this is the pickup t
    vol_transfer=vol_transfero
    exponential_decayr=exponential_decayro
    
    

## GRID FUNCTIONS
def update_boundaries():
    global _u1,_v1,_h1,_h2,_u2, _v2, u1, u2,h1,h2
    #Solid walls left and right
    #    - No zonal (u) flow through the left and right walls
    #    - Zero x-derivative in v and h
    # No flow through the boundary at x=0
    
   # _u[0, :] = 0
   # _u[1, :] = 0
   # _u[-1, :] = 0
   # _u[-2, :] = 0
    
    #PERIODIC WHERE THERE IS WIND
    _u1[0, tau_min+1:tau_max+1] = _u1[-3, tau_min+1:tau_max+1]
    _u1[1, tau_min+1:tau_max+1] = _u1[-2, tau_min+1:tau_max+1]
    _u1[-1, tau_min+1:tau_max+1] = _u1[2, tau_min+1:tau_max+1]
    
    _h1[0, tau_min+1:tau_max+1] = _h1[-2, tau_min+1:tau_max+1]
    _h1[-1, tau_min+1:tau_max+1] = _h1[1, tau_min+1:tau_max+1]

    
    _v1[0, tau_min:tau_max+2] = _v1[-2, tau_min:tau_max+2] #CHECK THIS IS RIGHT FOR V
    _v1[-1, tau_min:tau_max+2] = _v1[1, tau_min:tau_max+2]
    
    _u2[0, tau_min+1:tau_max+1] = _u2[-3, tau_min+1:tau_max+1]
    _u2[1, tau_min+1:tau_max+1] = _u2[-2, tau_min+1:tau_max+1]
    _u2[-1, tau_min+1:tau_max+1] = _u2[2, tau_min+1:tau_max+1]
    
    _h2[0, tau_min+1:tau_max+1] = _h2[-2, tau_min+1:tau_max+1]
    _h2[-1, tau_min+1:tau_max+1] = _h2[1, tau_min+1:tau_max+1]

    
    _v2[0, tau_min:tau_max+2] = _v2[-2, tau_min:tau_max+2]
    _v2[-1, tau_min:tau_max+2] = _v2[1, tau_min:tau_max+2]   

    #vs on east and west walls are 0 (no slip). This is adcroft/marshall implementation
    _v1[0, 0:tau_min+1] = - _v1[1, 0:tau_min+1]
    _v1[-1, 0:tau_min+1] = -_v1[-2, 0:tau_min+1]
    
    _v2[0, 0:tau_min+1] = -_v2[1, 0:tau_min+1]
    _v2[-1, 0:tau_min+1] = -_v2[-2, 0:tau_min+1]
    
    _v1[0, tau_max:-1] = -_v1[1, tau_max:-1]
    _v1[-1, tau_max:-1] = -_v1[-2, tau_max:-1]
    
    _v2[0, tau_max:-1] = -_v2[1, tau_max:-1]
    _v2[-1, tau_max:-1] = -_v2[-2, tau_max:-1]
    
    #hs on east/west walls: zero-derivative (free slip). Again this doesn't matter because we dont actually form 2nd derivs for h
    
    _h1[0, 0:tau_min+1] = _h1[1, 0:tau_min+1]
    _h1[-1, 0:tau_min+1] = _h1[-2, 0:tau_min+1]
    
    _h2[0, 0:tau_min+1] = _h2[1, 0:tau_min+1]
    _h2[-1, 0:tau_min+1] = _h2[-2, 0:tau_min+1]
    
    _h1[0, tau_max:-1] = _h1[1, tau_max:-1]
    _h1[-1, tau_max:-1] = _h1[-2, tau_max:-1]
    
    _h2[0, tau_max:-1] = _h2[1, tau_max:-1]
    _h2[-1, tau_max:-1] = _h2[-2, tau_max:-1]

    
    #us on east and west walls are 0
    
    _u1[0, 0:tau_min+1] = 0
    _u1[1, 0:tau_min+1] = 0
    _u1[-1, 0:tau_min+1] = 0
    _u1[-2, 0:tau_min+1] = 0
    
    _u2[0, 0:tau_min+1] = 0
    _u2[1, 0:tau_min+1] = 0
    _u2[-1, 0:tau_min+1] = 0
    _u2[-2, 0:tau_min+1] = 0
    
    _u1[0, tau_max:-1] = 0
    _u1[1, tau_max:-1] = 0
    _u1[-1, tau_max:-1] = 0
    _u1[-2, tau_max:-1] = 0
    
    _u2[0, tau_max:-1] = 0
    _u2[1, tau_max:-1] = 0
    _u2[-1, tau_max:-1] = 0
    _u2[-2, tau_max:-1] = 0
    
    # No flow through at top and bottom and no slip (us also 0)
    _v1[:, 0] = 0
    _v1[:, 1] = 0
    _v1[:, -1] = 0
    _v1[:, -2] = 0

    _v2[:, 0] = 0
    _v2[:, 1] = 0
    _v2[:, -1] = 0
    _v2[:, -2] = 0
    
    _u1[:, 0] = -_u1[:, 1]
    _u1[:, -1]=-_u1[:, -2]

    _u2[:, 0] = -_u2[:, 1]
    _u2[:, -1]=-_u2[:, -2]
    
    #h free slip. again does not matter 
    _h1[:, 0] = _h1[:, 1]
    _h1[:, -1]=_h1[:, -2]

    _h2[:, 0] = _h2[:, 1]
    _h2[:, -1]= _h2[:, -2]

    state1=np.array([u1,v1, h1, u2, v2,h2],dtype=object)
    # This applied for both boundary cases above
    for field in state1:

    	# fix corners to be average of neighbours
        field[0, 0] =  0.5*(field[1, 0] + field[0, 1])
        field[-1, 0] = 0.5*(field[-2, 0] + field[-1, 1])
        field[0, -1] = 0.5*(field[1, -1] + field[0, -2])
        field[-1, -1] = 0.5*(field[-1, -2] + field[-2, -1])
        
        
def update_boundaries_intermediate(bdd_u1,bdd_v1,bdd_h1,bdd_u2,bdd_v2,bdd_h2,u1_next, h1_next, u2_next, h2_next): #update boundaries so the intermediate rk steps obey
    #PERIODIC WHERE THERE IS WIND
    bdd_u1[0, tau_min+1:tau_max+1] = bdd_u1[-3, tau_min+1:tau_max+1]
    bdd_u1[1, tau_min+1:tau_max+1] = bdd_u1[-2, tau_min+1:tau_max+1]
    bdd_u1[-1, tau_min+1:tau_max+1] = bdd_u1[2, tau_min+1:tau_max+1]
    
    bdd_h1[0, tau_min+1:tau_max+1] = bdd_h1[-2, tau_min+1:tau_max+1]
    bdd_h1[-1, tau_min+1:tau_max+1] = bdd_h1[1, tau_min+1:tau_max+1]

    
    bdd_v1[0, tau_min:tau_max+2] = bdd_v1[-2, tau_min:tau_max+2] #CHECK THIS IS RIGHT FOR V
    bdd_v1[-1, tau_min:tau_max+2] = bdd_v1[1, tau_min:tau_max+2]
    
    bdd_u2[0, tau_min+1:tau_max+1] = bdd_u2[-3, tau_min+1:tau_max+1]
    bdd_u2[1, tau_min+1:tau_max+1] = bdd_u2[-2, tau_min+1:tau_max+1]
    bdd_u2[-1, tau_min+1:tau_max+1] = bdd_u2[2, tau_min+1:tau_max+1]
    
    bdd_h2[0, tau_min+1:tau_max+1] = bdd_h2[-2, tau_min+1:tau_max+1]
    bdd_h2[-1, tau_min+1:tau_max+1] = bdd_h2[1, tau_min+1:tau_max+1]

    
    bdd_v2[0, tau_min:tau_max+2] = bdd_v2[-2, tau_min:tau_max+2]
    bdd_v2[-1, tau_min:tau_max+2] = bdd_v2[1, tau_min:tau_max+2]   

    #vs on east and west walls are no slip
    bdd_v1[0, 0:tau_min+1] = - bdd_v1[1, 0:tau_min+1]
    bdd_v1[-1, 0:tau_min+1] = - bdd_v1[-2, 0:tau_min+1]
    
    bdd_v2[0, 0:tau_min+1] = - bdd_v2[1, 0:tau_min+1]
    bdd_v2[-1, 0:tau_min+1] = - bdd_v2[-2, 0:tau_min+1]
    
    bdd_v1[0, tau_max:-1] = -bdd_v1[1, tau_max:-1]
    bdd_v1[-1, tau_max:-1] = -bdd_v1[-2, tau_max:-1]
    
    bdd_v2[0, tau_max:-1] = -bdd_v2[1, tau_max:-1]
    bdd_v2[-1, tau_max:-1] = -bdd_v2[-2, tau_max:-1]
    
    #hs on east/west walls: zero-derivative (free slip)
    
    bdd_h1[0, 0:tau_min+1] = bdd_h1[1, 0:tau_min+1]
    bdd_h1[-1, 0:tau_min+1] = bdd_h1[-2, 0:tau_min+1]
    
    bdd_h2[0, 0:tau_min+1] = bdd_h2[1, 0:tau_min+1]
    bdd_h2[-1, 0:tau_min+1] = bdd_h2[-2, 0:tau_min+1]
    
    bdd_h1[0, tau_max:-1] = bdd_h1[1, tau_max:-1]
    bdd_h1[-1, tau_max:-1] = bdd_h1[-2, tau_max:-1]
    
    bdd_h2[0, tau_max:-1] = bdd_h2[1, tau_max:-1]
    bdd_h2[-1, tau_max:-1] = bdd_h2[-2, tau_max:-1]

    
    #us on east and west walls are 0
    
    bdd_u1[0, 0:tau_min+1] = 0
    bdd_u1[1, 0:tau_min+1] = 0
    bdd_u1[-1, 0:tau_min+1] = 0
    bdd_u1[-2, 0:tau_min+1] = 0
    
    bdd_u2[0, 0:tau_min+1] = 0
    bdd_u2[1, 0:tau_min+1] = 0
    bdd_u2[-1, 0:tau_min+1] = 0
    bdd_u2[-2, 0:tau_min+1] = 0
    
    bdd_u1[0, tau_max:-1] = 0
    bdd_u1[1, tau_max:-1] = 0
    bdd_u1[-1, tau_max:-1] = 0
    bdd_u1[-2, tau_max:-1] = 0
    
    bdd_u2[0, tau_max:-1] = 0
    bdd_u2[1, tau_max:-1] = 0
    bdd_u2[-1, tau_max:-1] = 0
    bdd_u2[-2, tau_max:-1] = 0
    
    # No flow through at top and bottom and no slip us and h free slip
    bdd_v1[:, 0] = 0
    bdd_v1[:, 1] = 0
    bdd_v1[:, -1] = 0
    bdd_v1[:, -2] = 0

    bdd_v2[:, 0] = 0
    bdd_v2[:, 1] = 0
    bdd_v2[:, -1] = 0
    bdd_v2[:, -2] = 0
    
    bdd_u1[:, 0] = - bdd_u1[:, 1]
    bdd_u1[:, -1] = - bdd_u1[:, -2]

    bdd_u2[:, 0] = -bdd_u2[:, 1]
    bdd_u2[:, -1] = -bdd_u2[:, -2]
    
    bdd_h1[:, 0] = bdd_h1[:, 1]
    bdd_h1[:, -1] = bdd_h1[:, -2]

    bdd_h2[:, 0] = bdd_h2[:, 1]
    bdd_h2[:, -1] = bdd_h2[:, -2]
    
    state2=np.array([bdd_u1[1:-1, 1:-1],bdd_v1[1:-1, 1:-1],bdd_h1[1:-1, 1:-1], bdd_u2[1:-1, 1:-1], bdd_v2[1:-1, 1:-1], bdd_h2[1:-1, 1:-1]],dtype=object) #this needs to be changed maybe
    # This applied for both boundary cases above
    for field in state2:
    	# fix corners to be average of neighbours
        field[0, 0] =  0.5*(field[1, 0] + field[0, 1])
        field[-1, 0] = 0.5*(field[-2, 0] + field[-1, 1])
        field[0, -1] = 0.5*(field[1, -1] + field[0, -2])
        field[-1, -1] = 0.5*(field[-1, -2] + field[-2, -1])



def diffx(psi):
   # """Calculate ∂/∂x[psi] over a single grid square.

    #i.e. d/dx(psi)[i,j] = (psi[i+1/2, j] - psi[i-1/2, j]) / dx

    #The derivative is returned at x points at the midpoint between
    #x points of the input array."""
    global dx
    return (psi[1:,:] - psi[:-1,:]) / dx

def diff2x(psi):
    #"""Calculate ∂2/∂x2[psi] over a single grid square.

    #i.e. d2/dx2(psi)[i,j] = (psi[i+1, j] - psi[i, j] + psi[i-1, j]) / dx^2

    #The derivative is returned at the same x points as the
    #x points of the input array, with dimension (nx-2, ny)."""
    global dx
    return (psi[:-2, :] - 2*psi[1:-1, :] + psi[2:, :]) / dx**2

def diff2y(psi):
   # """Calculate ∂2/∂y2[psi] over a single grid square.

   # i.e. d2/dy2(psi)[i,j] = (psi[i, j+1] - psi[i, j] + psi[i, j-1]) / dy^2

    #The derivative is returned at the same y points as the
   # y points of the input array, with dimension (nx, ny-2)."""
    global dy
    return (psi[:, :-2] - 2*psi[:, 1:-1] + psi[:, 2:]) / dy**2

def diffy(psi):
   # """Calculate ∂/∂y[psi] over a single grid square.

   # i.e. d/dy(psi)[i,j] = (psi[i, j+1/2] - psi[i, j-1/2]) / dy

   # The derivative is returned at y points at the midpoint between
   # y points of the input array."""
    global dy
    return (psi[:, 1:] - psi[:,:-1]) / dy

def centre_average(phi):
    #"""Returns the four-point average at the centres between grid points."""
    return 0.25*(phi[:-1,:-1] + phi[:-1,1:] + phi[1:, :-1] + phi[1:,1:])

def y_average(phi):
   # """Average adjacent values in the y dimension.
   # If phi has shape (nx, ny), returns an array of shape (nx, ny - 1)."""
    return 0.5*(phi[:,:-1] + phi[:,1:])

def x_average(phi):
   # """Average adjacent values in the x dimension.
   # If phi has shape (nx, ny), returns an array of shape (nx - 1, ny)."""
    return 0.5*(phi[:-1,:] + phi[1:,:])

def divergence1():
   # """Returns the top layer horizontal divergence at h1 points."""
    return diffx(u1) + diffy(v1)

def divergence2():
   # """Returns the bottom layer horizontal divergence at h2 points."""
    return diffx(u2) + diffy(v2)

def del2(phi):
   # """Returns the Laplacian of phi."""
    return diff2x(phi)[:, 1:-1] + diff2y(phi)[1:-1, :]

def uvatuv():
   # """Calculate the value of us at vs and vs at us."""
    global _u1, _v1, _u2, _v2
    u1bar = centre_average(_u1)[1:-1, :]
    v1bar = centre_average(_v1)[:, 1:-1]
    u2bar = centre_average(_u2)[1:-1, :]
    v2bar = centre_average(_v2)[:, 1:-1]
    return u1bar, v1bar, u2bar, v2bar

def uvath():
    global u1, v1, u2, v2
    u1bar = x_average(u1)
    v1bar = y_average(v1)
    u2bar = x_average(u2)
    v2bar = y_average(v2)
    return u1bar, v1bar

def absmax(psi):
    return np.max(np.abs(psi))

## Model functions
def damping(var):
    # sponges are active at the top and bottom of the domain by applying Rayleigh friction
    # with exponential decay towards the centre of the domain
    #relax is the number we are going to relax to
    #This relaxes for both ends of the domain
    global sponge, sponge_ny
    sponge_ny = ny//7
    sponge = np.exp(-np.linspace(0, exponential_decayr, sponge_ny))
    var_sponge = np.zeros_like(var)
    var_sponge[:, :sponge_ny] = sponge[np.newaxis, :]
    var_sponge[:, -sponge_ny:] = sponge[np.newaxis, ::-1]
    return var_sponge*var

def rhs():
    global u1, v1, h1, h2, u2, v2, _u1,_v1,_h1,_h2,_u2, _v2, m1, m2, extent,dt,t,start_t
   # Calculate the right hand side of the equations (u,v,h) in both layers
    u1_at_v, v1_at_u, u2_at_v, v2_at_u = uvatuv()   # (nx, ny+1), (nx+1, ny)

    
    
    # the u equations
    dh1dx = diffx(_h1)[:, 1:-1]  # (nx+1, ny)
    dh2dx = diffx(_h2)[:, 1:-1]  # (nx+1, ny)
    dh1dy  = diffy(_h1)[1:-1, :]   # (nx, ny+1)
    dh2dy  = diffy(_h2)[1:-1, :]   # (nx, ny+1)
    # the height equations
    h1_rhs = -h1*divergence1()-x_average(u1)*x_average(dh1dx)-y_average(v1)*y_average(dh1dy) #-u1*dh1dx #*u1 #this is part of the nonlinear expression
    h2_rhs = -h2*divergence2()-x_average(u2)*x_average(dh2dx)-y_average(v2)*y_average(dh2dy)
    
    if extent>0:
        if np.size(m1)==1:
            h1_rhs[:,ny//7:ny//7+extent] += m1 #adding here an optional mass input term
            h2_rhs[:,ny//7:ny//7+extent] += m2 #adding here an optional mass input term
        if np.size(m1)>1:
            h1_rhs[:,ny//7:ny//7+extent] += m1[int((t-start_t)/dt)]
            #print(m1[int((t-start_t)/dt)])
            #print(start_t)
        if np.size(m2)>1:
            h2_rhs[:,ny//7:ny//7+extent] += m2[int((t-start_t)/dt)]
    
    h1_rhs+=vol_transfer
    h2_rhs+=-vol_transfer
    
    

    
    u1_rhs = (f0 + beta*uy)*v1_at_u - g1*dh1dx -g2*(dh1dx+dh2dx) + nu*del2(_u1) - r*damping(u1) + tau_x/(1024*H1) 
    u2_rhs = (f0 + beta*uy)*v2_at_u - (rho1/rho2)*g2*(dh1dx+dh2dx)+ nu*del2(_u2) - r*damping(u2)
    
    # the v equations
  
    
    v1_rhs = -(f0 + beta*vy)*u1_at_v - g1*dh1dy -g2*(dh1dy+dh2dy) + nu*del2(_v1) - r*damping(v1) 
    v2_rhs = -(f0 + beta*vy)*u2_at_v - (rho1/rho2)*g2*(dh1dy+dh2dy) + nu*del2(_v2) - r*damping(v2)
    
    k1=np.array([u1_rhs, v1_rhs, h1_rhs, u2_rhs, v2_rhs, h2_rhs],dtype=object)
    
    def uvatuv_given(bdd_u1,bdd_v1,bdd_u2,bdd_v2):
        u1bar = centre_average(bdd_u1)[1:-1, :]
        v1bar = centre_average(bdd_v1)[:, 1:-1]
        u2bar = centre_average(bdd_u2)[1:-1, :]
        v2bar = centre_average(bdd_v2)[:, 1:-1]
        return u1bar, v1bar, u2bar, v2bar
    
    #form the next k functions
    def k(k_prev,step):
        u1_next=u1+k_prev[0]*step
        v1_next=v1+k_prev[1]*step
        h1_next=h1+k_prev[2]*step
        u2_next=u2+k_prev[3]*step
        v2_next=v2+k_prev[4]*step
        h2_next=h2+k_prev[5]*step
        #now form the updated bdd terms
        bdd_u1_next=copy.deepcopy(_u1)
        bdd_u1_next[1:-1, 1:-1]=u1_next
        bdd_v1_next=copy.deepcopy(_v1)
        bdd_v1_next[1:-1, 1:-1]=v1_next
        bdd_h1_next=copy.deepcopy(_h1)
        bdd_h1_next[1:-1, 1:-1]=h1_next
        bdd_u2_next=copy.deepcopy(_u2)
        bdd_u2_next[1:-1, 1:-1]=u2_next
        bdd_v2_next=copy.deepcopy(_v2)
        bdd_v2_next[1:-1, 1:-1]=v2_next
        bdd_h2_next=copy.deepcopy(_h2)
        bdd_h2_next[1:-1, 1:-1]=h2_next
        update_boundaries_intermediate(bdd_u1_next,bdd_v1_next,bdd_h1_next,bdd_u2_next,bdd_v2_next,bdd_h2_next,u1_next,h1_next,u2_next,h2_next)
        #now form the derivatives 
        dh1dx_next = diffx(bdd_h1_next)[:, 1:-1]  # (nx+1, ny)
        dh2dx_next = diffx(bdd_h2_next)[:, 1:-1]  # (nx+1, ny)
        dh1dy_next = diffy(bdd_h1_next)[1:-1, :]   # (nx, ny+1)
        dh2dy_next = diffy(bdd_h2_next)[1:-1, :]   # (nx, ny+1)
        u1_at_v_next, v1_at_u_next, u2_at_v_next, v2_at_u_next = uvatuv_given(bdd_u1_next,bdd_v1_next,bdd_u2_next,bdd_v2_next)   # (nx, ny+1), (nx+1, ny)
        #now form the next h equations
        h1_rhs_next = -h1_next*(diffx(u1_next) + diffy(v1_next)) -x_average(u1_next)*x_average(dh1dx_next)-y_average(v1_next)*y_average(dh1dy_next)
        h2_rhs_next = -h2_next*(diffx(u2_next) + diffy(v2_next))-x_average(u2_next)*x_average(dh2dx_next)-y_average(v2_next)*y_average(dh2dy_next)
        if extent>0:
            h1_rhs_next[:,ny//7:ny//7+extent] += m1 #adding here an optional mass input term
            h2_rhs_next[:,ny//7:ny//7+extent] += m2 #adding here an optional mass input term
        h1_rhs_next+=vol_transfer
        h2_rhs_next+=-vol_transfer
        #now form the u equations
        u1_rhs_next = (f0 + beta*uy)*v1_at_u_next - g1*dh1dx_next -g2*(dh1dx_next+dh2dx_next) + nu*del2(bdd_u1_next) - r*damping(u1_next) + tau_x/(1024*H1) 
        u2_rhs_next = (f0 + beta*uy)*v2_at_u_next - (rho1/rho2)*g2*(dh1dx_next+dh2dx_next)+ nu*del2(bdd_u2_next) - r*damping(u2_next)
        #now form the v equations
        v1_rhs_next = -(f0 + beta*vy)*u1_at_v_next - g1*dh1dy_next -g2*(dh1dy_next+dh2dy_next) + nu*del2(bdd_v1_next) - r*damping(v1_next) 
        v2_rhs_next = -(f0 + beta*vy)*u2_at_v_next - (rho1/rho2)*g2*(dh1dy_next+dh2dy_next) + nu*del2(bdd_v2_next) - r*damping(v2_next)
        
        return np.array([u1_rhs_next, v1_rhs_next, h1_rhs_next, u2_rhs_next, v2_rhs_next, h2_rhs_next],dtype=object)
    
    k2=k(k1,dt/2)
    k3=k(k2,dt/2)
    k4=k(k3,dt)

    return (k1+2*k2+2*k3+k4)/6


def step():
    global dt, t, _u1,_v1,_h1,_h2,_u2, _v2

    update_boundaries()

    dstate = rhs()

    state = np.array([u1, v1, h1, u2, v2, h2],dtype=object)
    # timestep forward
    newstate = state + dt*dstate
    u1[:], v1[:], h1[:], u2[:], v2[:], h2[:] = newstate
    

    t  += dt
    return u1,v1,h1,u2,v2,h2,t