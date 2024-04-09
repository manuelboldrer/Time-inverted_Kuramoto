import numpy as np


def circular_undirected_adjacency_matrix(N):
    adjacency_matrix = np.zeros((N, N))
    for i in range(N):
        adjacency_matrix[i, (i + 1) % N] = 1
        adjacency_matrix[(i + 1) % N, i] = 1
    return adjacency_matrix

def kuramoto(x, K, N, Omega, Ad):
    # Your implementation of the Kuramoto model here
    x = x.reshape(-1,1)
    f = Omega + (K / N) * np.sum(np.sin(x * Ad - Ad.T * x.T), axis=1)
    return f

def fun(X, alpha):
    X = alpha * X
    t = np.mod(X, 2 * np.pi)
    x = 35 * np.cos(7*t)
    y = 60 * np.sin(6*t)
    z = 5 * np.sin(5*t)
    return x,y,z

# Parameters
N = 13
p = 6
Ad = circular_undirected_adjacency_matrix(N) 
A = 35
B = 60
neigh = [np.where(Ad[i] == 1)[0] for i in range(N)]
iter = 2000
alpha = 1
rs = np.sin(np.pi/N)*np.sqrt(35**2+60**2)

# Initializations
theta = np.zeros((N, iter+1))
Kappa = 1 * N * np.ones((1, N))
Omega = .1 / alpha * np.ones((N,iter+1))
ff = np.zeros((N, iter+1))
x, y, z = np.zeros((iter+1, N)), np.zeros((iter+1, N)), np.zeros((iter+1, N))
xnext, ynext, znext = np.zeros((iter+1, N)), np.zeros((iter+1, N)), np.zeros((iter+1, N))

# Initialize theta
for j in range(N):
    theta[j, 0] = (j-1) * 2 * np.pi * p / N #+ 3 * np.random.rand()
# Time integration
dt = 0.1

co = 1
xX, yY, zZ = [], [], []
for c in np.arange(0, 2 * np.pi, 0.001):
    xX.append(fun(c, 1)[0])
    yY.append(fun(c, 1)[1])
    zZ.append(fun(c, 1)[2])
 


for kk in range(iter):
    for i in range(N):
        counte = 0
        for j in range(N):
            if Ad[i, j] == 1:
                counte += 1
                neigh[i][counte-1] = j

    for j in range(N):
        x[kk, j], y[kk, j], z[kk, j] = fun(theta[j, kk], alpha)

    k1 = kuramoto(theta[:, kk], Kappa, N, Omega[:,kk], Ad)
    k2 = kuramoto(theta[:, kk] + 0.5 * dt * k1, Kappa, N, Omega[:,kk], Ad)
    k3 = kuramoto(theta[:, kk] + 0.5 * dt * k2, Kappa, N, Omega[:,kk], Ad)
    k4 = kuramoto(theta[:, kk] + dt * k3, Kappa, N, Omega[:,kk], Ad)
    theta[:, kk+1] = theta[:, kk] + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    ff[:, kk+1] = (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    for j in range(N):

        xnext[kk, j], ynext[kk, j], znext[kk, j] = fun(theta[j, kk+1], alpha)

        vel = np.sqrt((xnext[kk,j]-x[kk,j])**2 +(ynext[kk,j]-y[kk,j])**2 +(ynext[kk,j]-y[kk,j])**2)
       
