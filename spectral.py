import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA
from pylab import *


def Lap(A):
  L = -copy(A)
  n = len(A)
  for i in range(n):
    d = sum(A[i])
    L[i, i] = d
  # print(L)
  return L


def eigvect(M):
  vals, vects = LA.eig(M)
  idx = vals.argsort()[::-1]
  vects = vects[:, idx]
  print("end lapeig")
  return vects[:, -2]


# Normalized laplacian eigenvector. Untested.
def normalized_lap_eigvect(A):
  print("start lapeig")
  L = Lap(A)
  n = len(A)
  # D inverse half
  Dih = zeros((n, n))
  for i in range(n):
    Dih[i, i] = L[i, i] ** (-1 / 2)
  return Dih @ eigvect(Dih @ L @ Dih)


def spectral(A):
  return thres(A, normalized_lap_eigvect(A))


# The epsilon ball graph, augmented with R and R2 for C and S.
# Neck, or the trough height, is needed to find R.
def eps_R_MCS(points, R, eps):
  n = len(points)
  M = zeros((n, n))
  C = zeros((n, n))
  S = zeros((n, n))
  for i in range(n):
    for j in range(n):
      if i != j and norm(points[i] - points[j]) <= eps:
        # iR = abs(points[i][0]) + neck - abs(points[i][1])
        # jR = abs(points[j][0]) + neck - abs(points[j][1])
        # Cheap hack, just how far from the neck are you
        # NOTE: THESE ARE SET TO 1 TO SEE IF THE REST OF THE CODE WORKS
        # iR = jR = 1
        C[i, j] = (R[i] + R[j]) / 2
        S[i, j] = ((R[i] + R[j]) / 2) ** 2
        # C[i, j] = (R[i] * R[j]) ** (1 / 2)
        # S[i, j] = R[i] * R[j]
        # Will get added to M[j, j] since i and j will change places.
        M[i, i] += 1
  print(Lap(C))
  print("Epsilon R MCS graph connected:", connected(C))
  return M, C, S


def plot_2d_points(points, file_name):
  x = []
  y = []
  for point in points:
    x.append(point[0])
    y.append(point[1])
  axes = plt.gca()
  # axes.set_aspect(0.35)
  axes.set_aspect(1)
  plt.scatter(x, y, color="black", s=4)
  plt.savefig(file_name)


def plot_3d_points(points, file_name):
  x = []
  y = []
  z = []
  for point in points:
    x.append(point[0])
    y.append(point[1])
    z.append(point[2])
  fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
  ax.scatter(x, y, z, color="black")

  ax.set(xticklabels=[], yticklabels=[], zticklabels=[])
  plt.savefig(file_name)


def plot_3d_point_clusters(points, cluster, file_name):
  n = len(points)
  x0 = []
  y0 = []
  z0 = []
  x1 = []
  y1 = []
  z1 = []
  for i in range(n):
    point = points[i]
    if i in cluster:
      x0.append(point[0])
      y0.append(point[1])
      z0.append(point[2])
    else:
      x1.append(point[0])
      y1.append(point[1])
      z1.append(point[2])
  fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
  ax.scatter(x0, y0, z0, color="red")
  ax.scatter(x1, y1, z1, color="black")

  ax.set(xticklabels=[], yticklabels=[], zticklabels=[])
  plt.savefig(file_name)


def plot_2d_point_clusters(points, cluster, file_name):
  n = len(points)
  x0 = []
  y0 = []
  x1 = []
  y1 = []
  for i in range(n):
    point = points[i]
    if i in cluster:
      x0.append(point[0])
      y0.append(point[1])
    else:
      x1.append(point[0])
      y1.append(point[1])
  axes = plt.gca()
  # axes.set_aspect(0.35)
  axes.set_aspect(1)
  plt.scatter(x0, y0, color="red", s=4)
  plt.scatter(x1, y1, color="black", s=4)
  plt.savefig(file_name)


def spectral_points_sqrttrough_MCS(points, valley_height, eps):
  R = [abs(point[0]) + valley_height for point in points]
  M, C, S = eps_R_MCS(points, R, eps)
  return spectralMCS(M, C, S)


def spectral_points_bowtie_MCS(points, neck, eps):
  R = [abs(point[0]) + neck - abs(point[1]) for point in points]
  # R = [abs(point[0]) + neck for point in points]
  print(R[(len(R) - 1) // 2])
  print("R", R)
  M, C, S = eps_R_MCS(points, R, eps)
  return spectralMCS(M, C, S)


# Mass matrix is always set to sum of degreees in this case.
# M, C and S are full matrices here. C and S are adjacency! M is diagonal.
def spectralMCS(M, C, S):
  return thresMC(M, C, lap_eigvectMS(M, S))


def IncMC(M, C, s, v, cut, vol):
  n = len(M)
  for i in range(n):
    if i in s:
      cut -= C[v][i]
    else:
      cut += C[v][i]
  vol += M[v, v]
  s.add(v)
  return (cut, vol)


def thresMC(M, C, v):
  order = argsort(v)
  sparsest = 10**9
  vol = 0
  cut = 0
  s = set()
  total_vol = sum(M)
  sparsest_i = -100
  for i, v in enumerate(order[:-1]):
    cut, vol = IncMC(M, C, s, v, cut, vol)
    current_sparsity = cut / min(vol, total_vol - vol)
    if current_sparsity < sparsest:
      sparsest = current_sparsity
      sparsest_threshold, sparsest_i = v, i
  return set(order[: sparsest_i + 1]), sparsest


# Solves Sv = M lambda_2 v in v.
def lap_eigvectMS(M, S):
  # M inverse halfpower
  n = len(M)
  LS = Lap(S)
  Mih = zeros((n, n))
  for i in range(n):
    Mih[i, i] = M[i, i] ** (-1 / 2)
  return Mih @ eigvect(Mih @ LS @ Mih)


def thres(G, v):
  # Volume measured in sum of edges on each side. Vol of a single-vertex side is
  # 0.
  def Inc(G, s, v, cut, vol):
    n = len(G)
    for i in range(n):
      if i in s:
        vol += G[v][i]
        cut -= G[v][i]
      else:
        cut += G[v][i]
    s.add(v)
    return (cut, vol)

  order = argsort(v)
  sparsest = 10**9
  vol = 0
  cut = 0
  s = set()
  total_vol = sum(G) / 2
  # Don't need sparsest_threshold, but sparsest_i by itself is confusing.
  sparsest_threshold = 0
  sparsest_i = -100
  # Don't actually need
  for i, v in enumerate(order[:-1]):
    # print(v, s)
    cut, vol = Inc(G, s, v, cut, vol)
    current_sparsity = cut / min(2 * vol + cut, 2 * (total_vol - vol) - cut)
    # print(cut, vol, total_vol - vol, current_sparsity)
    if current_sparsity < sparsest:
      # print(current_sparsity, cut, vol, total_vol - vol - cut)
      sparsest = current_sparsity
      sparsest_threshold, sparsest_i = v, i
  return set(order[: sparsest_i + 1]), sparsest


def accuracy(set, truth, n):
  nset = list(range(n))
  setC = nset.difference(set)
  truthC = nset.difference(truth)
  return max(
      set.intersection(truth) + setC.intersection(truthC),
      set.intersection(truthC) + setC.intersection(truth),
  )


def sample():
  return


def get_eps():
  return


def connected(G):
  visited = set()
  dfs(G, visited, 0)
  return len(visited) == len(G)


def dfs(G, visited, current):
  n = len(G)
  if current not in visited:
    visited.add(current)
    neighbors = [i for i in range(n) if G[current][i] > 0]
    for i in neighbors:
      dfs(G, visited, i)
  return
