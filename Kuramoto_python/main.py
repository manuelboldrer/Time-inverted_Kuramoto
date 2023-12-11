import numpy as np
import matplotlib.pyplot as plt



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
    x = np.cos(t)#*np.cos(5*t)
    y = np.sin(t)#*np.sin(5*t)
    z = np.cos(t)

    return x,y,z

# Parameters
N = 5
p = 0
Ad = circular_undirected_adjacency_matrix(N) 

neigh = [np.where(Ad[i] == 1)[0] for i in range(N)]
iter = 2000
alpha = 1
rs = 0.45

# Initializations
theta = np.zeros((N, iter+1))
K = 1 * N * np.ones((1, N))
Omega = 0.1 / alpha * np.ones((N,iter+1))
ff = np.zeros((N, iter+1))
x, y, z = np.zeros((iter+1, N)), np.zeros((iter+1, N)), np.zeros((iter+1, N))

# Initialize theta
for j in range(N):
    theta[j, 0] = (j-1) * 2 * np.pi * p / N + 3 * np.random.rand()

# Time integration
dt = 0.1


plt.ion()
fig, ax = plt.subplots()

for kk in range(iter):
    for i in range(N):
        counte = 0
        for j in range(N):
            if Ad[i, j] == 1:
                counte += 1
                neigh[i][counte-1] = j

    for j in range(N):
        x[kk, j], y[kk, j], z[kk, j] = fun(theta[j, kk], alpha)

    k1 = kuramoto(theta[:, kk], K, N, Omega[:,kk], Ad)
    k2 = kuramoto(theta[:, kk] + 0.5 * dt * k1, K, N, Omega[:,kk], Ad)
    k3 = kuramoto(theta[:, kk] + 0.5 * dt * k2, K, N, Omega[:,kk], Ad)
    k4 = kuramoto(theta[:, kk] + dt * k3, K, N, Omega[:,kk], Ad)
    theta[:, kk+1] = theta[:, kk] + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    ff[:, kk+1] = (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


    if kk == 0:  # Plots
        co = 1
        xX, yY, zZ = [], [], []
        for c in np.arange(0, 2 * np.pi, 0.001):
            xX.append(fun(c, 1)[0])
            yY.append(fun(c, 1)[1])
            zZ.append(fun(c, 1)[2])
    ax.plot(xX, yY, 'r', linewidth=1)
    a1 = plt.plot(x[kk, :], y[kk, :], 'o', markersize=8)
    for i in range(len(Ad)):
        for j in range(i + 1, len(Ad[i])):
            if Ad[i, j] == 1:
                plt.plot([x[kk,i], x[kk,j]], [y[kk,i], y[kk,j]], 'o-')
            
    plt.axis('equal')
    plt.draw()
    plt.pause(0.01)
    ax.clear()
