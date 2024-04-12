import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import pdist
def minimum_distance_among_robots(x, y, z):
    # Concatenate x, y, and z into a single array
    positions = np.column_stack((x, y, z))
    # Compute pairwise distances between all combinations of robots
    distances = pdist(positions)
    # Find the minimum distance among all combinations of robots
    min_distance = np.min(distances)
    return min_distance

def circular_undirected_adjacency_matrix(N):
    adjacency_matrix = np.zeros((N, N))
    for i in range(N):
        adjacency_matrix[i, (i + 1) % N] = 1
        adjacency_matrix[(i + 1) % N, i] = 1
    return adjacency_matrix

def kuramoto(x, K, N, Omega, Ad):
    # Your implementation of the Kuramoto model here
    x = x.reshape(-1,1)
    f = Omega + (K ) * np.sum(np.sin(x * Ad - Ad.T * x.T), axis=1)
    return f

def fun(X, alpha):
    X = alpha * X
    t = np.mod(X, 2 * np.pi)
    x = 20 * np.cos(9*t)
    y = 20 * np.sin(5*t)#*np.sin(5*t)
    z = 7 + 2 * np.sin(1*t)
    return x,y,z

# Parameters
N = 14
p = 5
Ad = circular_undirected_adjacency_matrix(N) 
A = 20
B = 20
neigh = [np.where(Ad[i] == 1)[0] for i in range(N)]
iter = 2000
alpha = 1
rs = np.sin(np.pi/N)*np.sqrt(A**2+B**2)

# Initializations
theta = np.zeros((N, iter+1))
K = 4* np.ones((1, N))
Omega = (0.04/ alpha ) * np.ones((N,iter+1))+ np.random.normal(loc=0, scale=0.01, size=(N, iter+1))
ff = np.zeros((N, iter+1))
x, y, z = np.zeros((iter+1, N)), np.zeros((iter+1, N)), np.zeros((iter+1, N))
xnext, ynext, znext = np.zeros((iter+1, N)), np.zeros((iter+1, N)), np.zeros((iter+1, N))

# Initialize theta
for j in range(N):
    theta[j, 0] = np.mod((j-1) * 2 * np.pi * p / N,2*np.pi) #+ np.random.rand

#theta[1,0]= theta[1,0] - 0.05
#theta[9,0]= theta[9,0] + 0.05


# Time integration
dt = 0.1

plt.ion()
fig, ax= plt.subplots()
ax = fig.add_subplot(111, projection='3d')

co = 1
xX, yY, zZ = [], [], []
for c in np.arange(0, 2 * np.pi, 0.001):
    xX.append(fun(c, 1)[0])
    yY.append(fun(c, 1)[1])
    zZ.append(fun(c, 1)[2])
        

xb = [-A, A, A, -A, -A]
yb = [-B, -B, B, B, -B]

# Plot the square

for kk in range(iter):
    
    
    
    #if kk == 30:
    #    Omega = (0.0/ alpha ) * np.ones((N,iter+1))+ np.random.normal(loc=0, scale=0.01, size=(N, iter+1))

    for i in range(N):
        counte = 0
        for j in range(N):
            if Ad[i, j] == 1:
                counte += 1
                neigh[i][counte-1] = j

    for j in range(N):
        x[kk, j], y[kk, j], z[kk, j] = fun(theta[j, kk], alpha)
    #print(x[kk, :] , y[kk, :], z[kk,:], "theta", theta[:,kk])

    k1 = kuramoto(theta[:, kk], K, N, Omega[:,kk], Ad)
    k2 = kuramoto(theta[:, kk] + 0.5 * dt * k1, K, N, Omega[:,kk], Ad)
    k3 = kuramoto(theta[:, kk] + 0.5 * dt * k2, K, N, Omega[:,kk], Ad)
    k4 = kuramoto(theta[:, kk] + dt * k3, K, N, Omega[:,kk], Ad)
    theta[:, kk+1] = theta[:, kk] + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    #theta[j,kk+1] = theta[j,kk]

    
    #for j in range(0,50):
    #    theta[j,kk+1] = theta[j,kk]
       

    ff[:, kk+1] = (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    for j in range(N):

        xnext[kk, j], ynext[kk, j], znext[kk, j] = fun(theta[j, kk+1], alpha) 

        vel = np.sqrt((xnext[kk,j]-x[kk,j])**2 +(ynext[kk,j]-y[kk,j])**2 +(ynext[kk,j]-y[kk,j])**2)
        #if j <6:
    #circle = plt.Circle([x[kk, j] , y[kk, j]], rs, color='b', fill=False)
        #else:
    #for j in range(0,50):
    #circle = plt.Circle([x[kk, j] , y[kk, j]], rs, color='r', fill=False) 
    TT = np.linspace(0, 2 * np.pi, 100)  # Angle parameter
    #ax.add_artist(circle)
    for j in range(N):
        radius = math.sin(math.pi/N)*math.sqrt(800)
        #print("radius : ", radius)
    # Parametric equations for a circle in 3D space
        XX = x[kk, j] + radius * np.cos(TT)
        YY = y[kk, j] + radius * np.sin(TT)
        ZZ = np.zeros_like(TT)  # Constant z-coordinate (in this case, at zero)

    # Plot the circle

        ax.plot(XX, YY, ZZ, color='b')
    
    
    ax.plot(xb, yb, color='k')
    ax.plot(xX, yY, zZ, 'r', linewidth=1)
    for j in range(N):
            ax.scatter(x[kk, j] , y[kk, j], z[kk,j], c='b', marker='o')
    
    min_distance = minimum_distance_among_robots(x[kk, :] , y[kk, :], z[kk,:])
    print("Minimum distance among all pairs of robots:", min_distance)
    print("time:", dt*kk)

    #for i in range(len(Ad)):
    #    for j in range(i + 1, len(Ad[i])):
    #        if Ad[i, j] == 1:
    #            plt.plot([x[kk,i], x[kk,j]], [y[kk,i], y[kk,j]], 'o-')
    #print(theta[:,kk])
    plt.axis('equal')
    plt.show()
    plt.pause(0.01)
    ax.clear()
