   # -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 10:26:14 2023

@author: matte
"""

import matplotlib
matplotlib.use('TkAgg')

import numpy as np
import sympy
from sympy import init_printing
init_printing(use_latex = True)
import math
import matplotlib.pyplot as plt

#Plots the 2D velocity field
from mpl_toolkits import mplot3d
from matplotlib import cm

plt.close('all')

'12 steps to NS equations'
#%% -----------------------------'Step 1 - 1D linear convection'-------------------------------------------------------------

'Step 1 - 1D linear convection'

xmax = 2
tmax = 0.5
nx   = 51               # Number of grid points
nt   = 151              # Number of time steps
dt   = tmax/(nt-1)      # Lenght of each time step 
dx   = xmax/(nx-1)      # Distance between two grid points
c    = 0.5               #Constant transoport velocity
u    = np.ones((nx,nt))
x    = np.linspace(0,xmax,nx)

'Boundary conditions'
u[0,:]    = 1
u[nx-1,:] = 1

'Initial conditions'

u[int(.5/dx):int(1/dx+1)]= 2

'''
# The same thing can be done, more awkwardly, this way:
for i in range(0,nx-1):
    if 0.5 <= x[i-1] <= 1:
        u[i,0] = 2
    else:
        u[i,0] = 1
'''
for n in range(0,nt-1):
    un = u.copy()
    for j in range(0,nx-1):
        u[j,n+1]=un[j,n]-c*dt/dx*(un[j,n]-un[j-1,n])

plt.figure(1)
plt_step_size = 10
for i in range(0,nt,plt_step_size):
      plt.plot(x, u[:,i], 'r')
      plt.xlabel('x')
      plt.ylabel('u')
      plt.ylim([0.8,2.2])
      plt.pause(0.1)
      plt.show()
      plt.figure(1,clear = True)
      # print(int(i/plt_step_size))
plt.close('all')

#%% -----------------------------'Step 2 - 1D convection'---------------------------------------------------------

'Step 2 - 1D convection'

p    = 4
xmax = p*2
tmax = p*0.5
nx   = p*51               # Number of grid points
nt   = p*151              # Number of time steps
dt   = tmax/(nt-1)        # Lenght of each time step 
dx   = xmax/(nx-1)        # Distance between two grid points
u    = np.ones((nx,nt))
x    = np.linspace(0,xmax,nx)

'Boundary conditions'
u[0,:]    = 1
u[nx-1,:] = 1

'Initial conditions'

u[int(.5/dx):int(1/dx+1)]= 2
# u[int(1/dx)::,0]= 0


'Upwind'
for n in range(0,nt-1):
    un = u.copy()
    for j in range(0,nx-1):
        u[j,n+1]=un[j,n]-un[j,n]*dt/dx*(un[j,n]-un[j-1,n])
'''
        
'Lax-Friederichs'
for n in range(0,nt-1):
    un = u.copy()
    for j in range(0,nx-1):
        u[j,n+1]=(un[j-1,n]+un[j+1,n])/2-un[j,n]*dt/2/dx*(un[j+1,n]-un[j-1,n])
'''

plt.figure(1)
plt_step_size = 10
for i in range(0,nt,plt_step_size):
      plt.plot(x, u[:,i], 'r')
      plt.xlabel('x')
      plt.ylabel('u')
      plt.ylim([0.8,2.5])
      plt.xlim([0,5])
      plt.pause(0.1)
      plt.show()
      plt.figure(1, clear=True)
      # print(int(i/plt_step_size))
plt.close('all')

#%% -----------------------------'Step 3 - 1D diffusion'--------------------------------------------------------

'Step 3 - 1D diffusion'
 
xmax = 2
tmax = 2.5
nx   = 21               # Number of grid points
nt   = 101              # Number of time steps
dt   = tmax/(nt-1)      # Lenght of each time step 
dx   = xmax/(nx-1)      # Distance between two grid points
vis  = 0.05             #Flow viscosity
u    = np.ones((nx,nt))
x    = np.linspace(0,xmax,nx)

'Boundary conditions'

u[0,:]    = 1
u[nx-1,:] = 1

'Initial conditions'

u[int(.5/dx):int(1/dx+1)]= 2


for n in range(0,nt-1):
    un = u
    for j in range(1,nx-1):
        u[j,n+1]=un[j,n]+vis*dt/dx/dx*(un[j+1,n]-2*un[j,n]+un[j-1,n])

import matplotlib.pyplot as plt
plt.figure(1)
plt_step_size = 10
for i in range(0,nt,plt_step_size):
      plt.plot(x, u[:,i], 'r')
      plt.xlabel('x')
      plt.ylabel('u')
      plt.ylim([0.8,2.2])
      plt.pause(0.1)
      plt.show()
      plt.figure(1, clear=True)
      # print(int(i/plt_step_size))
plt.close('all')

#%% -----------------------------'Step 4 - 1D Burgers equation'-------------------------------------------------------
'''
'Step 4 "- 1D Burgers" equation'

xmax = 2
tmax = 0.5
nx = 51               # Number of grid points
nt = 151              # Number of time steps
dt = tmax/(nt-1)      # Lenght of each time step 
dx = xmax/(nx-1)      # Distance between two grid points
vis = 0.1             # Flow viscosity

x,nu,t=sympy.symbols('x nu t')  # symbolic variables
phi=(sympy.exp(-(x-4*t)**2)/(4*nu*(t+1)))+sympy.exp(-(x-4*t-2*sympy.pi)**2/(4*nu*(t+1)))


for k in range(0,nx) 
    x[k]=dx*k

'Boundary conditions'

u[0,:] = 1
u[nx-1,:] = 1

'Initial conditions'

for i in range(0,nx-1):
    if 0.5 <= x[i-1] <= 1:
        u[i,0] = 2
    else:
        u[i,0] = 1
        
for i in range(0,nx):
    ip1[i] = i+2
    im1[i] = i
    x[i] = (i-1)*dx
    
ip1[nx-1] = 1
im1[0] = nx

for i in range(0,nx-1):
    phi = math.exp(-x[i]**2/4/vis)+math.exp(-(x[i]-2*math.pi)**2/4/vis)
    dphi = -0.5/vis*x[i]*math.exp(-x[i]**2/4/vis)-0.5/vis*(x[i]-2*math.pi)*math.exp(-(x[i]-2*math.pi)**2/4/vis)
    u[i,0] = -2*vis*dphi/phi+4

'Soluzione analitica esatta'
for it in range(1,nt):
    t = (it-1)*dt
    ua = u
    for i in range(1,nx):
        phi = math.exp(-(x[i]-4*t)**2/4/vis/(t+1))+math.exp(-(x[i]-4*t-2*math.pi)**2/4/vis/(t+1))
        dphi = -0.5/vis/(t+1)*(x[i]-4*t)*math.exp(-(x[i]-4*t)**2/4/vis/(t+1))-0.5/vis/(t+1)*(x[i]-4*t-2*math.pi)*math.exp(-(x[i]-4*t-2*math.pi)**2/4/vis/(t+1))
        ua[i] = -2*vis*dphi/phi+4
   
'Soluzione numerica'   
for n in range(0,nt-1):
    un = u
    for i in range(45,nx-1):
        u[i,n+1] = un[i,n]-un[i,n]*dt/dx*(un[i,n]-un[int(im1[i])-1,n])+vis*dt/dx**2*(un[int(ip1[i])-1,n]-2*un[i,n]+un[int(im1[i])-1,n])          
     
import matplotlib.pyplot as plt
plt.figure()
for i in range(0,nt,10):
      plt.plot(x, u[:,i], 'r')
      plt.xlabel('x (m)')
      plt.ylabel('u (m/s)')
      #plt.ylim([0,2.2])
      plt.show()
'''

'Step 4 "- 1D Burgers" equation'

p    = 4
xmax = p*2
tmax = p*0.5
nx   = p*51                    # Number of grid points
nt   = p*151                   # Number of time steps
dt   = tmax/(nt-1)             # Lenght of each time step 
dx   = xmax/(nx-1)             # Distance between two grid points
u    = np.ones((nx,nt))
x    = np.linspace(0,xmax,nx)
vis  = 0.2                     # Flow viscosity  (if vis = 0 we get back to the step 3)

'Initial conditions'

# u[int(.5/dx):int(1/dx+1)]= 2
u[int(1/dx)::,0] = 0
        
for n in range(0,nt-1):
    un = u.copy()
    for j in range(0,nx-1):
        u[j,n+1]=un[j,n]-un[j,n]*dt/dx*(un[j,n]-un[j-1,n])+vis*dt/dx/dx*(un[j+1,n]-2*un[j,n]+un[j-1,n])

plt.figure(1)
plt_step_size = 10
for i in range(0,nt,plt_step_size):
    plt.plot(x, u[:,i], 'r')
    plt.xlabel('x')
    plt.ylabel('u')
    plt.ylim([-0.5,1.5])
    plt.xlim([0,5])
    plt.pause(0.1)
    plt.show()
    plt.figure(1, clear=True)
    # print(int(i/plt_step_size))
plt.close('all')


#%% -----------------------------'Step 5 - 2D linear convection'--------------------------------------------------

'Step 5 - 2D linear convection'

xmax = 2
ymax = 2
tmax = 2
nx   = 31               
ny   = 31
nt   = 101              
dt   = tmax/(nt-1)      
dx   = xmax/(nx-1)     
dy   = ymax/(ny-1)     
c    = 1
u    = np.zeros((nx,ny,nt))
v    = np.zeros((nx,ny,nt))
x    = np.zeros(nx)
y    = np.zeros(ny)

for k in range(0,nx):
    x[k] = dx*k

for p in range(0,ny):
    y[p] = dy*p


u[0,:,:]    = 0
u[nx-1,:,:] = 0
u[:,0,:]    = 0
u[:,ny-1,:] = 0

for i in range(0,nx):
    for j in range(0,ny):
        if 0.5 <= x[i] <= 1 and 0.5 <= y[j] <= 1:
            u[i,j]=2
        else:
            u[i,j]=1
            
# 2D upwind method
for n in range(0,nt-1):
    un=u
    for i in range(0,nx-1):
        for j in range(0,ny-1):      
            u[i,j,n+1]  = un[i,j,n]-c*dt/dx*(un[i,j,n]-un[i-1,j,n])-c*dt/dy*(un[i,j,n]-un[i,j-1,n])

   
#Plots the 2D velocity field
from mpl_toolkits import mplot3d
from matplotlib import cm
import matplotlib.pyplot as plt

plt_step_size = 1
plt.figure(1,figsize=(14, 9))
for i in range(0,nt,plt_step_size):
    ax = plt.axes(projection ='3d')
    # Make data.
    X = np.arange(-5, 5, 10/nx)
    Y = np.arange(-5, 5, 10/ny)
    X, Y = np.meshgrid(X, Y)
    
    # Plot the surface.
    surf = ax.plot_surface(X, Y, u[:,:,i], cmap=cm.viridis, linewidth=0, antialiased=False)
    plt.pause(0.1)
    plt.show()
    plt.figure(1, clear=True, figsize=(14, 9))
    #print(int(i / plt_step_size))
plt.close('all')

#%% -----------------------------'Step 6 - 2D convection'-----------------------------------------------
'Step 6 - 2D convection'

xmax = 2
ymax = 2
tmax = 2
nx   = 31               
ny   = 31
nt   = 101              
dt   = tmax/(nt-1)      
dx   = xmax/(nx-1)     
dy   = ymax/(ny-1)     
c    = 1
u    = np.zeros((nx,ny,nt))
v    = np.zeros((nx,ny,nt))
x    = np.zeros(nx)
y    = np.zeros(ny)

for k in range(0,nx):
    x[k] = dx*k

for p in range(0,ny):
    y[p] = dy*p


u[0,:,:]    = 0
u[nx-1,:,:] = 0
u[:,0,:]    = 0
u[:,ny-1,:] = 0

v[0,:,:]    = 0
v[nx-1,:,:] = 0
v[:,0,:]    = 0
v[:,ny-1,:] = 0

for i in range(0,nx):
    for j in range(0,ny):
        if 0.5 <= x[i] <= 1 and 0.5 <= y[j] <= 1:
            u[i,j]=2
            v[i,j]=2
        else:
            u[i,j]=1
            v[i,j]=1
            
# 2D upwind method
for n in range(0,nt-1):
    un=u
    vn=v
    for i in range(0,nx-1):
        for j in range(0,ny-1):      
            u[i,j,n+1]  = un[i,j,n]-un[i,j,n]*dt/dx*(un[i,j,n]-un[i-1,j,n])-v[i,j,n+1]*dt/dy*(un[i,j,n]-un[i,j-1,n])
            v[i,j,n+1]  = vn[i,j,n]-v[i,j,n+1]*dt/dx*(vn[i,j,n]-vn[i-1,j,n])-un[i,j,n]*dt/dy*(vn[i,j,n]-vn[i,j-1,n])
     
#Plots the 2D velocity field
from mpl_toolkits import mplot3d
from matplotlib import cm
import matplotlib.pyplot as plt

plt_step_size = 1
plt.figure(1,figsize=(14, 9))
for i in range(0,nt,plt_step_size):
    # Make data.
    X = np.arange(-5, 5, 10/nx)
    Y = np.arange(-5, 5, 10/ny)
    X, Y = np.meshgrid(X, Y)
    # Plot the surface.
    ax = plt.axes(projection='3d')
    surf = ax.plot_surface(X, Y,u[:,:,i], cmap=cm.viridis, linewidth=0, antialiased=False)
    plt.pause(0.1)
    plt.show()
    plt.figure(1,clear = True,figsize=(14, 9))
    print(int(i/plt_step_size))
plt.close('all')


#%% -----------------------------'Step 7 - 2D Diffusion'-----------------------------------------------------------
'Step 7 - 2D Diffusion'

xmax = 2
ymax = 2
tmax = 2
nx   = 21               
ny   = 21
nt   = 101              
dt   = tmax/(nt-1)      
dx   = xmax/(nx-1)     
dy   = ymax/(ny-1)     
vis  = 0.1
u    = np.zeros((nx,ny,nt))
x    = np.zeros(nx)
y    = np.zeros(ny)

for k in range(0,nx):
    x[k] = dx*k

for p in range(0,ny):
    y[p] = dy*p

u[0,:,:]    = 0
u[nx-1,:,:] = 0
u[:,0,:]    = 0
u[:,ny-1,:] = 0

for i in range(0,nx):
    for j in range(0,ny):
        if 0.5 <= x[i] <= 1 and 0.5 <= y[j] <= 1:
            u[i,j]=2
        else:
            u[i,j]=1

# 2D 
for n in range(0,nt-1):
    un=u
    for i in range(0,nx-1):
        for j in range(0,ny-1):      
            u[i,j,n+1]  = un[i,j,n]+vis*dt/dx/dx*(un[i+1,j,n]-2*un[i,j,n]+un[i-1,j,n])+vis*dt/dy/dy*(un[i,j+1,n]-2*un[i,j,n]+un[i,j-1,n])
            
#Plots the 2D velocity field
from mpl_toolkits import mplot3d
from matplotlib import cm
import matplotlib.pyplot as plt

plt_step_size = 1
plt.figure(1, figsize =(14, 9))
for i in range(0,nt,plt_step_size):
    ax = plt.axes(projection ='3d')
    # Make data.
    X = np.arange(-5, 5, 10/nx)
    Y = np.arange(-5, 5, 10/ny)
    X, Y = np.meshgrid(X, Y)
    ax.set_zlim([1, 2])   # Cambia questi limiti se il plot esce dalle Z
    # Plot the surface.
    surf = ax.plot_surface(X, Y,u[:,:,i], cmap=cm.viridis, linewidth=0, antialiased=False)
    plt.pause(0.1)
    plt.show()
    plt.figure(1, clear=True, figsize=(14, 9))
    #print(int(i / plt_step_size))
plt.close('all')


#%% -----------------------------'Step 8 - 2D Burgers Equation'-----------------------------------------------------------
#%% -----------------------------'Step 9 - Laplace Equation'-----------------------------------------------------
'Step 9 - Laplace Equation'
    
xmax = 2
ymax = 1
tmax = 2
nx   = 21               
ny   = 21
nt   = 101              
dt   = tmax/(nt-1)      
dx   = xmax/(nx-1)     
dy   = ymax/(ny-1)     
vis  = 0.1
p    = np.zeros((ny,nx))
x    = np.linspace(0,xmax,nx)
y    = np.linspace(0,ymax,ny)

# Dirichlet omogeneo

p[:,0]  = 0
p[:,-1] = y

# Neumann omogeneo

p[0,:] = p[1,:]
p[-1,:]= p[-2,:]

plt.figure(1,figsize=(14, 9), dpi=100)
for n in range(0,nt-1):
    pn=p
    for i in range(0,nx-1):
        for j in range(0,ny-1):     
            p[i,j]  = ((pn[i+1,j]+pn[i-1,j])*dy**2+
                       (pn[i,j+1]+pn[i,j-1])*dx**2)/((dx**2+dy**2)*2)

    p[:,0]  = 0
    p[:,-1] = y    
    p[0,:]  = p[1,:]
    p[-1,:] = p[-2,:]

    ax = plt.axes(projection ='3d')
    
    # Make data.
    X, Y = np.meshgrid(x, y)
    ax.set_zlim([0, 1])   # Cambia questi limiti se il plot esce dalle Z
    
    # Plot the surface.
    surf = ax.plot_surface(X, Y,p[:], cmap=cm.viridis, linewidth=0, antialiased=False)
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 1)
    ax.view_init(30, 225)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')

    plt.figure(1, clear=True, figsize=(14, 9), dpi=100)
    plt.pause(0.1)
    plt.show()
    # print(n)
plt.close('all')

#%% -----------------------------'Step 10 - Poisson Equation'---------------------------------------------------
'Step 10 - Poisson Equation'
    
xmax = 2
ymax = 1
tmax = 2
nx   = 21               
ny   = 21
nt   = 101              
dt   = tmax/(nt-1)      
dx   = xmax/(nx-1)     
dy   = ymax/(ny-1)     
vis  = 0.1

p    = np.zeros((ny,nx))
b    = np.zeros((ny,nx))

x    = np.linspace(0,xmax,nx)
y    = np.linspace(0,ymax,ny)

b[int((nx-1)/4),int((ny-1)/4)]     = 100
b[int((nx-1)*3/4),int((ny-1)*3/4)] = -100

for n in range(0,nt-1):
    un=u
    for i in range(0,nx-1):
        for j in range(0,ny-1):     
            p[i,j]  = ((pn[i+1,j]+pn[i-1,j])*dy**2+
                       (pn[i,j+1]+pn[i,j-1])*dx**2-b[i,j]*dx**2*dy**2)/((dx**2+dy**2)*2)
    fig = plt.figure(figsize =(14, 9),dpi=100)
    ax  = plt.axes(projection ='3d')
    
    # Make data.
    X, Y = np.meshgrid(x, y)
    ax.set_zlim([0, 1])   # Cambia questi limiti se il plot esce dalle Z
    
    # Plot the surface.
    surf = ax.plot_surface(X, Y,p[:], cmap=cm.viridis, linewidth=0, antialiased=False)
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 1)
    ax.view_init(30, 225)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    

#%% -----------------------------'Step 11 - Cavity Flow'------------------------------------------------------
'Step 11 - Cavity Flow'
    
xmax = 2
ymax = 2
tmax = 2
nx   = 21               
ny   = 21
nt   = 501
nit  = 51  # number of iterations          
dt   = tmax/(nt-1)      
dx   = xmax/(nx-1)     
dy   = ymax/(ny-1)  
   
vis = 0.1
rho = 1

u = np.zeros((ny,nx))
v = np.zeros((ny,nx))
p = np.zeros((ny,nx))
b = np.zeros((ny,nx))

x = np.linspace(0,xmax,nx)
y = np.linspace(0,ymax,ny)

# Make data.
X, Y = np.meshgrid(x, y)

def sBrackets(b,rho,dt,dx,dy,nx,ny,u,v):

    for i in range(1,nx-1):
        for j in range(1,ny-1):
            b[i,j]=rho*(1/dt)*((u[i,j+1]-u[i,j-1])/(2*dx)+(v[i+1,j]-v[i-1,j])/(2*dy))-((u[i,j+1]-u[i,j-1])/(2*dx))**2-2*(u[i+1,j]-u[i-1,j])/(2*dy)*(v[i,j+1]-v[i,j-1])/(2*dx)-((v[i+1,j]-v[i-1,j])/(2*dy))**2

    return b

def Pressure(p,b,dx,dy,nx,ny,nit):
    
    for iit in range(0,nit):
        pn=p.copy()
        for i in range(1,nx-1):
            for j in range(1,ny-1):
                p[i,j]=((pn[i,j+1]+pn[i,j-1])*dy**2+(pn[i+1,j]+pn[i-1,j])*dx**2-b[i,j]*dx**2*dy**2)/(dx**2+dy**2)/2
                
        p[0,:]  = p[1,:];   # dp/dy = 0 @ x = 2
        p[:,0]  = p[:,1];   # dp/dx = 0 @ x = 0
        p[-1,:] = 0;        # p = 0 @ y = 2
        p[:,-1] = p[:,-2];  # dp/dx = 0 @ x = 2
 
    return p

for it in range(0,nt):
    
    b=sBrackets(b,rho,dt,dx,dy,nx,ny,u,v)
    
    p=Pressure(p,b,dx,dy,nx,ny,nit)
    raise ValueError
    un=u.copy()
    vn=v.copy()
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            u[i,j] = un[i,j]-un[i,j]*dt/dx*(un[i,j]-un[i,j-1])-vn[i,j]*dt/dy*(un[i,j]-un[i-1,j])-1/rho*(p[i,j+1]-p[i,j-1])*dt/2/dx+vis*dt/dx**2*(un[i,j-1]-2*un[i,j]+un[i,j+1])+vis*dt/dy**2*(un[i-1,j]-2*un[i,j]+un[i+1,j])
            v[i,j] = vn[i,j]-un[i,j]*dt/dx*(vn[i,j]-vn[i,j-1])-vn[i,j]*dt/dy*(vn[i,j]-vn[i-1,j])-1/rho*(p[i+1,j]-p[i-1,j])*dt/2/dy+vis*dt/dx**2*(vn[i,j-1]-2*vn[i,j]+vn[i,j+1])+vis*dt/dy**2*(vn[i-1,j]-2*vn[i,j]+vn[i+1,j])
    
    u[0,:]=0; u[-1,:]=1; u[:,0]=0; u[:,-1]=0;
    v[0,:]=0; v[-1,:]=0; v[:,0]=0; v[:,-1]=0;
     
    # Plot the surface.
    fig = plt.figure(figsize=(11, 7), dpi=100)
    plt.contourf(X, Y, p, alpha=0.5, cmap=cm.viridis)
    plt.colorbar()
    plt.contour(X, Y, p, cmap=cm.viridis)
    plt.streamplot(X, Y, u, v)
    plt.xlabel('X')
    plt.ylabel('Y')
    
#%% -----------------------------'Step 12 - Channel Flow'--------------------------------------------------------------
'Step 12 - Channel Flow'
    
xmax = 2
ymax = 2
tmax = 2
nx = 41               
ny = 41
nt = 100
nit=50  # number of iterations          
dt = tmax/(nt-1)      
dx = xmax/(nx-1)     
dy = ymax/(ny-1)  
   
vis = 0.1
rho=1
F=1

u = np.zeros((ny,nx))
v = np.zeros((ny,nx))
p = np.zeros((ny,nx))
b = np.zeros((ny,nx))

x = np.linspace(0,xmax,nx)
y = np.linspace(0,ymax,ny)

# Make data.
X, Y = np.meshgrid(x, y)

def sBrackets(b,rho,dt,dx,dy,nx,ny,u,v):

    for i in range(1,nx-1):
        for j in range(1,ny-1):
            b[i,j]=rho*(1/dt)*((u[i,j+1]-u[i,j-1])/(2*dx)+(v[i+1,j]-v[i-1,j])/(2*dy))-((u[i,j+1]-u[i,j-1])/(2*dx))**2-2*(u[i+1,j]-u[i-1,j])/(2*dy)*(v[i,j+1]-v[i,j-1])/(2*dx)-((v[i+1,j]-v[i-1,j])/(2*dy))**2
            
            # Periodic BC Pressure @ x = 2
            
            b[i,-1]=rho*(1/dt)*((u[i,0]-u[i,-2])/(2*dx)+(v[i+1,-1]-v[i-1,-1])/(2*dy))-((u[i,0]-u[i,-2])/(2*dx))**2-2*(u[i+1,-1]-u[i-1,-1])/(2*dy)*(v[i,0]-v[i,-2])/(2*dx)-((v[i+1,-1]-v[i-1,-1])/(2*dy))**2
            
            # Periodic BC Pressure @ x = 0
    
            b[i,0]=rho*(1/dt)*((u[i,1]-u[i,-1])/(2*dx)+(v[i+1,0]-v[i-1,0])/(2*dy))-((u[i,1]-u[i,-1])/(2*dx))**2-2*(u[i+1,0]-u[i-1,0])/(2*dy)*(v[i,1]-v[i,-1])/(2*dx)-((v[i+1,0]-v[i-1,0])/(2*dy))**2
    
    return b

def Pressure(p,b,dx,dy,nx,ny,nit):
    
    for iit in range(0,nit):
        pn=p.copy()
        for i in range(1,nx-1):
            for j in range(1,ny-1):
                p[i,j]=((pn[i,j+1]+pn[i,j-1])*dy**2+(pn[i+1,j]+pn[i-1,j])*dx**2-b[i,j]*dx**2*dy**2)/(dx**2+dy**2)/2
                
                # Periodic BC Pressure @ x = 2
                
                p[i,-1]=((pn[i,0]+pn[i,-2])*dy**2+(pn[i+1,-1]+pn[i-1,-1])*dx**2-b[i,-1]*dx**2*dy**2)/(dx**2+dy**2)/2
                
                # Periodic BC Pressure @ x = 0
                
                p[i,0]=((pn[i,1]+pn[i,-1])*dy**2+(pn[i+1,0]+pn[i-1,0])*dx**2-b[i,0]*dx**2*dy**2)/(dx**2+dy**2)/2
                
        p[0,:]=p[1,:]; # dp/dy = 0 @ y = 0
        p[-1,:]=p[-2,:];   # dp/dy = 0 @ y = 2
        
        return p

for it in range(0,nt):
    
    b=sBrackets(b,rho,dt,dx,dy,nx,ny,u,v)
    
    p=Pressure(p,b,dx,dy,nx,ny,nit)
    
    un=u.copy()
    vn=v.copy()
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            
            u[i,j] = un[i,j]-un[i,j]*dt/dx*(un[i,j]-un[i,j-1])-vn[i,j]*dt/dy*(un[i,j]-un[i-1,j])-1/rho*(p[i,j+1]-p[i,j-1])*dt/2/dx+vis*dt/dx**2*(un[i,j-1]-2*un[i,j]+un[i,j+1])+vis*dt/dy**2*(un[i-1,j]-2*un[i,j]+un[i+1,j])+F*dt
            v[i,j] = vn[i,j]-un[i,j]*dt/dx*(vn[i,j]-vn[i,j-1])-vn[i,j]*dt/dy*(vn[i,j]-vn[i-1,j])-1/rho*(p[i+1,j]-p[i-1,j])*dt/2/dy+vis*dt/dx**2*(vn[i,j-1]-2*vn[i,j]+vn[i,j+1])+vis*dt/dy**2*(vn[i-1,j]-2*vn[i,j]+vn[i+1,j])
            
            # Periodic BC u @ x = 2  
            
            u[i,-1] = un[i,-1]-un[i,-1]*dt/dx*(un[i,-1]-un[i,-2])-vn[i,-1]*dt/dy*(un[i,-1]-un[i-1,-1])-1/rho*(p[i,0]-p[i,-2])*dt/2/dx+vis*dt/dx**2*(un[i,-2]-2*un[i,-1]+un[i,0])+vis*dt/dy**2*(un[i-1,-1]-2*un[i,-1]+un[i+1,-1])+F*dt
            
            # Periodic BC u @ x = 0  
            
            u[i,0] = un[i,0]-un[i,0]*d8t/dx*(un[i,0]-un[i,-1])-vn[i,0]*dt/dy*(un[i,0]-un[i-1,0])-1/rho*(p[i,1]-p[i,-1])*dt/2/dx+vis*dt/dx**2*(un[i,-1]-2*un[i,0]+un[i,1])+vis*dt/dy**2*(un[i-1,0]-2*un[i,0]+un[i+1,0])+F*dt
            
            # Periodic BC u @ x = 2  
            
            v[i,-1] = vn[i,-1]-un[i,-1]*dt/dx*(vn[i,-1]-vn[i,-2])-vn[i,-1]*dt/dy*(vn[i,-1]-vn[i-1,-1])-1/rho*(p[i+1,-1]-p[i-1,-1])*dt/2/dy+vis*dt/dx**2*(vn[i,-2]-2*vn[i,-1]+vn[i,0])+vis*dt/dy**2*(vn[i-1,-1]-2*vn[i,-1]+vn[i+1,-1])
            
            # Periodic BC u @ x = 0  
            
            v[i,0] = vn[i,0]-un[i,0]*dt/dx*(vn[i,0]-vn[i,-1])-vn[i,0]*dt/dy*(vn[i,0]-vn[i-1,0])-1/rho*(p[i+1,0]-p[i-1,0])*dt/2/dy+vis*dt/dx**2*(vn[i,-1]-2*vn[i,0]+vn[i,1])+vis*dt/dy**2*(vn[i-1,0]-2*vn[i,0]+vn[i+1,0])
            
    u[0,:]=0; u[-1,:]=0;
    v[0,:]=0; v[-1,:]=0;
     
    # Plot the surface.
    fig = plt.figure(figsize=(11, 7), dpi=100)
    plt.quiver(X[::3, ::3], Y[::3, ::3], u[::3, ::3], v[::3, ::3])
    plt.xlabel('X')
    plt.ylabel('Y')


