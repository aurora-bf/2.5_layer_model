import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import copy
import sys
import os

# Get parent directory then add it
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

## CONFIGURATION
# Domain

#use this option to choose if want to start from 0 or load in a previous run and continue from there
start_from0=1 #0 if start from 0 and 1 if load in previous file

nx = 128
ny = 129

H1 = 250.0          # [m]  Average depth of the fluid
H2 = 250.0          # [m]  Average depth of the fluid
H3 = 4000.0

L_x_h =3795751.714012305        # [m]  Zonal width of domain. #see compute_domainsize.ipynb in 1 degree folder for how we got this
L_y_h = 13785809.584401697        # [m]  Meridional height of domain


L_x =(ny/(ny-1))*L_x_h     # [m]  Zonal width of domain. #see compute_domainsize.ipynb in 1 degree folder for how we got this
L_y = (nx/(nx-1))*L_y_h        # [m]  Meridional height of domain

dx = L_x / nx            # [m]
dy = L_y / ny            # [m]


#Set coriolis values
omega=7.2921E-5
ref_lat=0
f0=2*omega*np.sin(np.deg2rad(ref_lat))
beta=2*omega*np.cos(np.deg2rad(ref_lat))/(6.637E6)

# Diffusion and Friction
nu = 8e3        # [m^2.s^-1] Coefficient of diffusion
r = 1.0e-4          # Rayleigh damping at top and bottom of domain
dt = 125.0         # Timestep [s]


# Set densities and reduced gravity
rho0 = 1024.0 #reference density of fluid
rho3=1030 #density abyssal fluid
rho2=1028 #density bottom active layer
rho1=1026 #density top active layer

g1= 9.8*(rho2-rho1)/rho1 #reduced gravity 1
g2=9.8*(rho3-rho2)/rho1

#Set wind intensity and channel through which it is forced
tau_0=0 #amplitude wind stress
tau_min=10 #smallest grid point where there will be wind in channel
tau_max=10 #largest grid point where there will be wind in channel

#Set the strength of the mass input step function
#m1=0 #the strength of mass input in the upper layer
#m2=0 #the strength of mass input in the lower layer
#extent=0 #the number of grid cells that the mass input is placed over

#Initialize lists to save to
output_interval = 2100          # how often to save variables to lists
h1_list = list()
u1_list = list()
v1_list = list()
h2_list = list()
u2_list = list()
v2_list = list()
h3_list = list()
u3_list = list()
v3_list = list()
t_list=list()

## GRID
# Arakawa-C Grid:
#
# +-- v --+
# |       |    * (nx, ny)   h points at grid centres
# u   h   u    * (nx+1, ny) u points on vertical edges  (u[0] and u[nx] are boundary values)
# |       |    * (nx, ny+1) v points on horizontal edges
# +-- v --+
#
# Variables preceeded with underscore  (_u, _v, _h) include the boundary values,
# variables without (u, v, h) are only the values defined within the domain
_u1 = np.zeros((nx+3, ny+2))
_v1 = np.zeros((nx+2, ny+3))
_h1 = H1*np.ones((nx+2, ny+2))

_u2 = np.zeros((nx+3, ny+2))
_v2 = np.zeros((nx+2, ny+3))
_h2 = H2*np.ones((nx+2, ny+2))

_u3 = np.zeros((nx+3, ny+2))
_v3 = np.zeros((nx+2, ny+3))
_h3 = H3*np.ones((nx+2, ny+2))


u1 = _u1[1:-1, 1:-1]               # (nx+1, ny)
v1 = _v1[1:-1, 1:-1]               # (nx, ny+1)
h1 = _h1[1:-1, 1:-1]               # (nx, ny)

u2 = _u2[1:-1, 1:-1]               # (nx+1, ny)
v2 = _v2[1:-1, 1:-1]               # (nx, ny+1)
h2 = _h2[1:-1, 1:-1]               # (nx, ny)

u3 = _u3[1:-1, 1:-1]               # (nx+1, ny)
v3 = _v3[1:-1, 1:-1]               # (nx, ny+1)
h3 = _h3[1:-1, 1:-1]               # (nx, ny)

state = np.array([u1, v1, h1, u2, v2, h2, u3, v3, h3])


dx = L_x / nx            # [m]
dy = L_y / ny            # [m]

# positions of the value points in [m]
ux = (-L_x/2 + np.arange(nx+1)*dx)[:, np.newaxis]
vx = (-L_x/2 + dx/2.0 + np.arange(nx)*dx)[:, np.newaxis]

vy = (-L_y/2 + np.arange(ny+1)*dy)[np.newaxis, :]
uy = (-L_y/2 + dy/2.0 + np.arange(ny)*dy)[np.newaxis, :]

hx = vx
hy = uy

t = 0 #Initial time


#Define the zonal wind

tau_x=np.zeros(ny)
for i in range(tau_min,tau_max):
    tau_x[i] = tau_0*np.sin(np.pi*uy[0,i]/(L_y/((ny)/(15-5))))
    
plt.plot(tau_x,uy[0,:])
plt.xlabel('Magnitude wind stress')
plt.ylabel('y distance')
plt.title('Wind')


if start_from0==0:
    ## INITIAL CONDITIONS
    # Set the initial state of the model here 
    u0 = u1 * 0.0
    v0 = v1 * 0.0
    # set the variable fields to the initial conditions
    u1[:] = u0
    v1[:] = v0
    h1[:] = H1  
    u2[:] = u0
    v2[:] = v0
    h2[:] = H2  
    u3[:] = u0
    v3[:] = v0
    h3[:] = H3  

    #Initialize some lists to store output

    u1_list=[]
    h1_list=[]
    v1_list=[]
    u2_list=[]
    h2_list=[]
    v2_list=[]
    u3_list=[]
    h3_list=[]
    v3_list=[]

    t_list=[]

if start_from0==1: #if starting from a previous time step
    #Initialize some lists to store output

    u1_list=[]
    h1_list=[]
    v1_list=[]
    u2_list=[]
    h2_list=[]
    v2_list=[]
    u3_list=[]
    h3_list=[]
    v3_list=[]

    t_list=[]

    h1_pickup=xr.open_dataset("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/h1_xr_control_years0to10_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
    h2_pickup=xr.open_dataset("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/h2_xr_control_years0to10_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
    h3_pickup=xr.open_dataset("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/h3_xr_control_years0to10_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
    u1_pickup=xr.open_dataset("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/u1_xr_control_years0to10_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
    u2_pickup=xr.open_dataset("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/u2_xr_control_years0to10_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
    u3_pickup=xr.open_dataset("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/u3_xr_control_years0to10_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
    v1_pickup=xr.open_dataset("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/v1_xr_control_years0to10_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
    v2_pickup=xr.open_dataset("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/v2_xr_control_years0to10_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
    v3_pickup=xr.open_dataset("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/v3_xr_control_years0to10_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
    
    u1[:] = (np.array((u1_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1])
    v1[:] = (np.array((v1_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1])
    h1[:] = (np.array((h1_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1])  
    u2[:] = (np.array((u2_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1])
    v2[:] = (np.array((v2_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1])
    h2[:] = (np.array((h2_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1]) 
    u3[:] = (np.array((u3_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1])
    v3[:] = (np.array((v3_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1])
    h3[:] = (np.array((h3_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1]) 
    
    
    _u1[1:-1, 1:-1]=u1             # (nx+1, ny)
    _v1[1:-1, 1:-1]=v1               # (nx, ny+1)
    _h1[1:-1, 1:-1]=h1            # (nx, ny)

    _u2[1:-1, 1:-1]=u2              # (nx+1, ny)
    _v2[1:-1, 1:-1]=v2             # (nx, ny+1)
    _h2[1:-1, 1:-1]=h2
    
    _u3[1:-1, 1:-1]=u3              # (nx+1, ny)
    _v3[1:-1, 1:-1]=v3             # (nx, ny+1)
    _h3[1:-1, 1:-1]=h3
    t= ((h1_pickup.time)[-1])*(24*60*60)
    
    h1_list.append(np.array((h1_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1]) #append the list to include tlast time point of the previou step
    h2_list.append(np.array((h2_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1]) 
    h3_list.append(np.array((h3_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1]) 
    u1_list.append(np.array((u1_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1]) 
    u2_list.append(np.array((u2_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1]) 
    u3_list.append(np.array((u3_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1]) 
    v1_list.append(np.array((v1_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1]) 
    v2_list.append(np.array((v2_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1]) 
    v3_list.append(np.array((v3_pickup.to_array().squeeze()).transpose('x','y','time'))[:,:,np.shape(h1_pickup.to_array())[3]-1]) 

    
from model_true3layerboussinesq import *
import_initialvalues(dt, t, _u1,_v1,_h1,_h2,_u2, _v2,_u3,_v3,_h3,ny,dx,dy,tau_min,tau_max,H1,H2,H3,nu,r,tau_x,f0,beta,rho1,rho2,rho3,uy,vy,m1o=0,m2o=4.930484906118382e-08,extento=5)
nsteps = 8087040
for i in range(nsteps):
    u1_new,v1_new,h1_new,u2_new,v2_new,h2_new,u3_new,v3_new,h3_new,t=step()
    if i % output_interval == 0:
        u1_list.append(copy.deepcopy(np.array(u1_new)))
        h1_list.append(copy.deepcopy(np.array(h1_new)))
        v1_list.append(copy.deepcopy(np.array(v1_new)))
        u2_list.append(copy.deepcopy(np.array(u2_new)))
        h2_list.append(copy.deepcopy(np.array(h2_new)))
        v2_list.append(copy.deepcopy(np.array(v2_new)))
        u3_list.append(copy.deepcopy(np.array(u3_new)))
        h3_list.append(copy.deepcopy(np.array(h3_new)))
        v3_list.append(copy.deepcopy(np.array(v3_new)))
        t_list.append(copy.deepcopy(t))    
        
        
##Turn lists into xarrays of each variable with time as a dimension. Then plot different time slices

#Stack into 3d numpy array with time as dimension
h1_array=np.dstack(h1_list) 
u1_array=np.dstack(u1_list)
v1_array=np.dstack(v1_list)
h2_array=np.dstack(h2_list)
u2_array=np.dstack(u2_list)
v2_array=np.dstack(v2_list)
h3_array=np.dstack(h3_list)
u3_array=np.dstack(u3_list)
v3_array=np.dstack(v3_list)

#Convert numpy arrays into xarrays
if start_from0==1:

    h1_xr = xr.DataArray(h1_array[:,:,1:], coords={'x': hx[:,0],'y': hy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    h1_xr=h1_xr.transpose('y','x','time')

    u1_xr = xr.DataArray(u1_array[:,:,1:], coords={'x': ux[:,0],'y': uy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    u1_xr=u1_xr.transpose('y','x','time')


    v1_xr = xr.DataArray(v1_array[:,:,1:], coords={'x': vx[:,0],'y': vy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    v1_xr=v1_xr.transpose('y','x','time')

    h2_xr = xr.DataArray(h2_array[:,:,1:], coords={'x': hx[:,0],'y': hy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    h2_xr=h2_xr.transpose('y','x','time')

    u2_xr = xr.DataArray(u2_array[:,:,1:], coords={'x': ux[:,0],'y': uy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    u2_xr=u2_xr.transpose('y','x','time')


    v2_xr = xr.DataArray(v2_array[:,:,1:], coords={'x': vx[:,0],'y': vy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    v2_xr=v2_xr.transpose('y','x','time')
    
    
    h3_xr = xr.DataArray(h3_array[:,:,1:], coords={'x': hx[:,0],'y': hy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    h3_xr=h3_xr.transpose('y','x','time')

    u3_xr = xr.DataArray(u3_array[:,:,1:], coords={'x': ux[:,0],'y': uy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    u3_xr=u3_xr.transpose('y','x','time')


    v3_xr = xr.DataArray(v3_array[:,:,1:], coords={'x': vx[:,0],'y': vy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    v3_xr=v3_xr.transpose('y','x','time')
    
if start_from0==0:
    h1_xr = xr.DataArray(h1_array, coords={'x': hx[:,0],'y': hy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    h1_xr=h1_xr.transpose('y','x','time')

    u1_xr = xr.DataArray(u1_array, coords={'x': ux[:,0],'y': uy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    u1_xr=u1_xr.transpose('y','x','time')


    v1_xr = xr.DataArray(v1_array, coords={'x': vx[:,0],'y': vy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    v1_xr=v1_xr.transpose('y','x','time')

    h2_xr = xr.DataArray(h2_array, coords={'x': hx[:,0],'y': hy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    h2_xr=h2_xr.transpose('y','x','time')

    u2_xr = xr.DataArray(u2_array, coords={'x': ux[:,0],'y': uy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    u2_xr=u2_xr.transpose('y','x','time')


    v2_xr = xr.DataArray(v2_array, coords={'x': vx[:,0],'y': vy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    v2_xr=v2_xr.transpose('y','x','time')

    
    h3_xr = xr.DataArray(h3_array, coords={'x': hx[:,0],'y': hy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    h3_xr=h3_xr.transpose('y','x','time')

    u3_xr = xr.DataArray(u3_array, coords={'x': ux[:,0],'y': uy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    u3_xr=u3_xr.transpose('y','x','time')


    v3_xr = xr.DataArray(v3_array, coords={'x': vx[:,0],'y': vy[0,:],'time': np.array(t_list)/(24*60*60)}, dims=["x", "y", "time"])
    v3_xr=v3_xr.transpose('y','x','time')

h1_xr.to_netcdf("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/h1_xr_control_years10to42pt5_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
h2_xr.to_netcdf("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/h2_xr_control_years10to42pt5_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
h3_xr.to_netcdf("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/h3_xr_control_years10to42pt5_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
u1_xr.to_netcdf("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/u1_xr_control_years10to42pt5_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
u2_xr.to_netcdf("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/u2_xr_control_years10to42pt5_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
u3_xr.to_netcdf("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/u3_xr_control_years10to42pt5_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
v1_xr.to_netcdf("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/v1_xr_control_years10to42pt5_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
v2_xr.to_netcdf("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/v2_xr_control_years10to42pt5_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")
v3_xr.to_netcdf("/scratch/abf376/2.5layer_reducedgravitymodel/saved_runs/v3_xr_control_years10to42pt5_bottomlayerpert_fix_fromstationary_northofsponge_nohdamp_sealevel_tunedtomitgcm_nu8e3_dt125_rk4_nohdiffusion_noslip_properimplement_update_onlycorners_nonlinear_H1_250_H2_250_correcteddomainsizetruethreelayer_H34000.nc")