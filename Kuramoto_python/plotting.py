import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import pdist

def minimum_distance_among_robots(x, y, z):
    positions = np.column_stack((x, y, z))
    distances = pdist(positions)
    min_distance = np.min(distances)
    return min_distance

def circular_undirected_adjacency_matrix(N):
    adjacency_matrix = np.zeros((N, N))
    for i in range(N):
        adjacency_matrix[i, (i + 1) % N] = 1
        adjacency_matrix[(i + 1) % N, i] = 1
    return adjacency_matrix

def kuramoto(x, K, N, Omega, Ad):
    x = x.reshape(-1, 1)
    f = Omega + (K / N) * np.sum(np.sin(x * Ad - Ad.T * x.T), axis=1)
    return f

def fun(X, alpha):
    X = alpha * X
    t = np.mod(X, 2 * np.pi)
    x = 35 * np.cos(7 * t)
    y = 60 * np.sin(6 * t)
    z = 5 * np.sin(5 * t)
    return x, y, z

N = 13
iter = 2000
dt = 0.1
Ad = circular_undirected_adjacency_matrix(N) 
alpha = 1
Omega = .1 / alpha * np.ones((N, iter + 1))

theta = np.zeros((N, iter + 1))
K = 1 * N * np.ones((1, N))
x, y, z = np.zeros((iter + 1, N)), np.zeros((iter + 1, N)), np.zeros((iter + 1, N))
xnext, ynext, znext = np.zeros((iter + 1, N)), np.zeros((iter + 1, N)), np.zeros((iter + 1, N))

for j in range(N):
    theta[j, 0] = (j - 1) * 2 * np.pi * 6 / N

plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xX, yY, zZ = [], [], []
for c in np.arange(0, 2 * np.pi, 0.001):
    xX.append(fun(c, 1)[0])
    yY.append(fun(c, 1)[1])
    zZ.append(fun(c, 1)[2])
        
ax.plot(xX, yY, zZ, 'r', linewidth=1)

for kk in range(iter):
    for j in range(N):
        x[kk, j], y[kk, j], z[kk, j] = fun(theta[j, kk], alpha)

    k1 = kuramoto(theta[:, kk], K, N, Omega[:, kk], Ad)
    k2 = kuramoto(theta[:, kk] + 0.5 * dt * k1, K, N, Omega[:, kk], Ad)
    k3 = kuramoto(theta[:, kk] + 0.5 * dt * k2, K, N, Omega[:, kk], Ad)
    k4 = kuramoto(theta[:, kk] + dt * k3, K, N, Omega[:, kk], Ad)
    theta[:, kk + 1] = theta[:, kk] + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    for j in range(N):
        xnext[kk, j], ynext[kk, j], znext[kk, j] = fun(theta[j, kk + 1], alpha)
        vel = np.sqrt((xnext[kk, j] - x[kk, j])**2 + (ynext[kk, j] - y[kk, j])**2 + (znext[kk, j] - z[kk, j])**2)
        print(vel)

    ax.scatter(x[kk, :], y[kk, :], z[kk, :], c='b', marker='o')
    min_distance = minimum_distance_among_robots(x[kk, :], y[kk, :], z[kk, :])
    print("Minimum distance among all pairs of robots:", min_distance)
    plt.axis('equal')
    plt.show()
    plt.pause(0.01)
    ax.clear()
