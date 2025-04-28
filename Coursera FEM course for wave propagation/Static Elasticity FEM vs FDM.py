import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#FEM solution of static elasticity problem
#------------------------------------------------------------------
# Initialization of setup 
#------------------------------------------------------------------
nx = 20  # Number of element boundary points
u = np.zeros(nx) # Solution vector
f = np.zeros(nx) # source term vector
mu= 1 # Constant sheer force vector 

#------------------------------------------------------------------
#Element boundary points
#------------------------------------------------------------------
x = np.linspace(0, 1, nx) # x in [0,1]
h= x[2]-x[1] # Element size, constant for uniform mesh
#------------------------------------------------------------------
#Assembe stiffness matrix K_ij
#------------------------------------------------------------------
K = np.zeros((nx,nx)) # Stiffness matrix
for i in range(1,nx-1):
    for j in range(1,nx-1):
        if i==j:
            K[i,j]=2*mu/h
        elif j==i+1 or j== i-1:
            K[i,j]= -mu/h
        else:
            K[i,j]=0

#------------------------------------------------------------------
#Source term is a spike at i = 3*nx/4
f[int(3*nx/4)] = 1.0 

#Boundary at condition x = 0
u[0] = 0.15; f[1] = u[0]/h 

#Boundary at condition x = 1
u[nx-1] = 0.05 ; f[nx-2]= u[nx-1]/h

#------------------------------------------------------------------
#FEM solution
#------------------------------------------------------------------
u[1:nx-1] = np.linalg.inv(K[1:nx-1, 1:nx-1]) @ np.transpose(f[1:nx-1])

#------------------------------------------------------------------
#Plotting
#------------------------------------------------------------------
xfe = u
plt.plot(x, xfe, color= 'r', lw=2, label= 'Finite elements')
plt.title('Static Elasticity', size= 16)
plt.ylabel('Displacement $u(x)$', size= 16)
plt.xlabel('Position $x$', size= 16)    
plt.axis ([0, 1, 0.04, .28])
plt.legend()
plt.grid(True)
plt.show()

#FDM solution of static elasticity problem
#Poisson's equation with relaxation method
#------------------------------------------------------------------

nt= 500 #Number of iterations
iplot = 2 #Snapshot frequency

# non- zero boundary conditions
u = np.zeros(nx) # set u to zero
du = np.zeros(nx) # du/dx
f = np.zeros(nx) #forcing 

f[int(3*nx/4)] = 1./h
xfd = x.copy()  # Ensure xfd is used for FDM plotting
u = np.zeros_like(xfd)  # Initialize u with the same shape as xfd

# ---------------------------------------------------------------
# Initialize animated plot
# ---------------------------------------------------------------
plt.figure(figsize=(8,6))

line1 = plt.plot(x, xfe, color='r', lw=2, label='FE') 
line2 = plt.plot(xfd, u, color='k', ls='-.', label='FD relaxation')
plt.title('Static Elasticity with relaxation method', size=16)
plt.ylabel('Displacement, $u$', size=12)
plt.xlabel('Position, $x$', size=12)
plt.legend(loc=4)
plt.grid(True)

plt.ion()   # set interective mode
plt.show()
# ---------------------------------------------------------------
for it in range(nt):
    # Calculate the average of u (omit boundaries)
    for i in range(1, nx-1):
        du[i] =u [i+1] + u[i-1]
    u = 0.5*( f*h**2/mu + du )
    u[0] = 0.15    # Boundary condition at x=0
    u[nx-1] = 0.05 # Boundary condition at x=1
    fd = u.copy()
    
    # --------------------------------------   
    # Animation plot. Display both solutions
    if not it % iplot:
        for l in line2:
            l.remove()
        line1 = plt.plot(x, xfe, color='r', lw=2, label='FE')
        line2 = plt.plot(xfd, u, color='k', ls='-.', label='FD relaxation')  # Use updated u for FDM plot
        plt.legend(loc=4)
        plt.gcf().canvas.draw()
        plt.pause(0.01)  # Add a small pause to allow the plot to refresh
        plt.gcf().canvas.draw()  

    