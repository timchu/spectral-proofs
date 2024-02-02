from pylab import *
from spectral import *


def dumbbell(a, b):
  n = a + b
  A = zeros((n, n))
  for i in range(a):
    for j in range(a):
      if i != j:
        A[i, j] = 1
  for i in range(a, n):
    for j in range(a, n):
      if i != j:
        A[i, j] = 1
  A[a - 1, a] = A[a, a - 1] = 1
  return A


# a = 4
# b = 10
# D = dumbbell(a, b)
# print(D)
# print(spectral(D))


def eps_graph(A, eps):
  n = len(A)
  B = zeros((n, n))
  for i in range(n):
    for j in range(n):
      if i != j and A[i, j] <= eps:
        # Unit length graph.
        B[i, j] = 1
  return B


def knn(A):
  n = len(A)
  B = zeros((n, n))

  for i in range(n):
    knns = argsort(A[i])[: k + 1]
    for j in knns:
      B[i, j] += A[i, j]
  return B


def norm(point):
  val = 0
  for i in point:
    val += i**2
  return sqrt(val)


# Bidirectional knn, so in the case of two points there's a double weighted
# edge.
def unweighted_knn(points, k):
  n = len(points)
  A = zeros((n, n))
  for i in range(n):
    lengths = []
    for j in range(n):
      if i != j:
        lengths.append([norm(points[i] - points[j]), j])
    lengths = sorted(lengths)
    for q in range(k):
      length = lengths[q][0]
      j = lengths[q][1]
      A[i, j] += 1  # Unweighted knn.
      A[j, i] += 1
  return A


def affinity_graph(points):
  n = len(points)
  A = zeros((n, n))
  for i in range(n):
    for j in range(i + 1, n):
      A[i, j] = A[j, i] = norm(points[i] - points[j])
  return A
