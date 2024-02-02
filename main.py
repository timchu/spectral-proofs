import sys
from graph import *
from pylab import *
from sample import *
from spectral import *


def cluster_gaussian():
  n = 600
  mean1 = [0, 0]
  cov1 = [[3, 0], [0, 1]]
  mean2 = [6, 0]
  cov2 = [[3, 0], [0, 1]]
  samples = sample_mixture_two_gaussians_2d(n, mean1, cov1, mean2, cov2)
  eps = 3
  cluster, _ = spectral(eps_graph(affinity_graph(samples), eps))
  plot_2d_point_clusters(samples, cluster, "two_gaussians.png")


# cluster_gaussian()

print("Spectral MCS")


def cluster_sqrttrough_MCS_and_normal():
  # Sample sqrttrough
  n = 3000
  eps = 0.2
  # Required for DFS to not break Python's recurison limit.
  sys.setrecursionlimit(int(1.5 * n))
  valley_height = 1 / 90000
  depth = 10
  samples = sample_sqrttrough(n, valley_height, depth)
  plot_2d_points(samples, "sqrttrough.png")

  cluster, sparse = spectral_points_sqrttrough_MCS(samples, valley_height, eps)
  print(cluster, len(cluster), sparse)
  plot_2d_point_clusters(samples, cluster, "spectral_R_sqrttrough.png")

  cluster, sparse = spectral(eps_graph(affinity_graph(samples), eps))
  print(cluster, len(cluster), sparse)
  plot_2d_point_clusters(samples, cluster, "spectral_normal_sqrttrough.png")


def plot_ellipse():
  n = 500
  center_x = 1
  center_y = 2
  radius_x = 3
  radius_y = 4
  samples = sample_ellipse(n, center_x, center_y, radius_x, radius_y)
  plot_2d_points(samples, "ellipse.png")


def plot_two_ellipses():
  n = m = 1500
  (center_x1, center_y1) = (-4, 0)
  (center_x2, center_y2) = (4, 0)
  radius_x1 = radius_x2 = 4 - 0.1  # + 0.06
  radius_y1 = radius_y2 = 18
  samples1 = sample_ellipse(n, center_x1, center_y2, radius_x1, radius_y1)
  samples2 = sample_ellipse(n, center_x2, center_y2, radius_x2, radius_y2)
  samples = []
  samples += [[samples1[i][0], samples1[i][1]] for i in range(n)]
  samples += [[samples2[i][0], samples2[i][1]] for i in range(m)]
  samples = array(sorted(samples))
  plot_2d_points(samples, "two_ellipses.png")

  eps = 1
  cluster, sparse = spectral(eps_graph(affinity_graph(samples), eps))
  print(cluster, len(cluster), sparse)
  plot_2d_point_clusters(samples, cluster, "spectral_two_ellipses.png")


def plot_bowtie_rect():
  samples = bowtie_rect(2, 100, 4, 1, 0.5)
  print(samples)
  print(len(samples))
  plot_2d_points(samples, "bowtie_rect.png")

  cluster, sparse = spectral(unweighted_knn(samples, 4))
  print(cluster, len(cluster), sparse)
  plot_2d_point_clusters(samples, cluster, "spectral_bowtie_rect.png")


def plot_densitybridges():
  total_x = 3
  total_y = 7
  bridge_x = 0.8
  bridge_y = 0.2
  bridge_spacing = 2
  n = 1000
  samples = densitybridges(
      total_x, total_y, bridge_x, bridge_y, bridge_spacing, n
  )
  print(samples)
  print(len(samples))
  plot_2d_points(samples, "density_bridges.png")
  G = unweighted_knn(samples, 8)
  cluster, sparse = spectral(G)
  print(cluster, len(cluster), sparse)
  plot_2d_point_clusters(samples, cluster, "spectral_density_bridges.png")

  # Create R density graph
  RG = copy(G)
  N = len(samples)
  for i in range(N):
    if samples[i][0] < -bridge_x:
      for j in range(N):
        RG[i, j] *= sqrt(-samples[i][0] - bridge_x)
        RG[j, i] *= sqrt(-samples[i][0] - bridge_x)
    elif -bridge_x <= samples[i][0] <= bridge_x:
      for j in range(N):
        RG[i, j] *= sqrt(bridge_y)
        RG[j, i] *= sqrt(bridge_y)
    else:  # samples[i][0] > bridge_x
      for j in range(N):
        RG[i, j] *= sqrt(samples[i][0] - bridge_x)
        RG[j, i] *= sqrt(samples[i][0] - bridge_x)
  cluster, sparse = spectral(RG)
  print(cluster, len(cluster), sparse)
  plot_2d_point_clusters(samples, cluster, "R_spectral_density_bridges.png")


def plot_rectbridges():
  total_x = 3
  total_y = 7
  bridge_x = 0.8
  bridge_y = 0.2
  bridge_spacing = 2
  point_spacing = 0.2
  samples = rectbridges(
      total_x, total_y, bridge_x, bridge_y, bridge_spacing, point_spacing
  )
  print(samples)
  print(len(samples))
  plot_2d_points(samples, "bridges_rect.png")
  G = unweighted_knn(samples, 4)
  cluster, sparse = spectral(G)
  print(cluster, len(cluster), sparse)
  plot_2d_point_clusters(samples, cluster, "spectral_bridges_rect.png")

  # Create an R_weighted graph.
  RG = copy(G)
  N = len(samples)
  for i in range(N):
    if samples[i][0] < -bridge_x:
      for j in range(N):
        RG[i, j] *= sqrt(-samples[i][0] - bridge_x)
        RG[j, i] *= sqrt(-samples[i][0] - bridge_x)
    elif -bridge_x <= samples[i][0] <= bridge_x:
      for j in range(N):
        RG[i, j] *= sqrt(bridge_y)
        RG[j, i] *= sqrt(bridge_y)
    else:  # samples[i][0] > bridge_x
      for j in range(N):
        RG[i, j] *= sqrt(samples[i][0] - bridge_x)
        RG[j, i] *= sqrt(samples[i][0] - bridge_x)
  cluster, sparse = spectral(RG)
  print(cluster, len(cluster), sparse)
  plot_2d_point_clusters(samples, cluster, "R_spectral_bridges_rect.png")


def plot_rectbridge():
  samples = rectbridge(5, 40, 0.7, 3, 0.5)
  print(samples)
  print(len(samples))
  plot_2d_points(samples, "rectbridge.png")
  cluster, sparse = spectral(unweighted_knn(samples, 4))
  print(cluster, len(cluster), sparse)
  plot_2d_point_clusters(samples, cluster, "spectral_rectbridge.png")


# Plots two rectangles with (2, 20) heights, centeres (-1.5, 0), (1.5, 0), # bridged by a 1 by 1 rectanlge in the center.
def plot_rectangles():
  length_x = 2
  length_y = 60
  small_length_x = 1
  small_length_y = 1
  center_x = (length_x + small_length_x) / 2
  n_bigrect = 900
  n_smallrect = int(
      n_bigrect * small_length_x * small_length_y / (length_x * length_y) + 1
  )
  samples1 = sample_rectangle(n_bigrect, center_x, 0, length_x, length_y)
  samples2 = sample_rectangle(n_bigrect, -center_x, 0, length_x, length_y)
  samples3 = sample_rectangle(n_smallrect, 0, 0, small_length_x, small_length_y)
  samples = []
  samples += [[samples1[i][0], samples1[i][1]] for i in range(n_bigrect)]
  samples += [[samples2[i][0], samples2[i][1]] for i in range(n_bigrect)]
  samples += [[samples3[i][0], samples3[i][1]] for i in range(n_smallrect)]
  samples = array(sorted(samples))
  plot_2d_points(samples, "rectangles.png")

  # cluster, sparse = spectral(eps_graph(affinity_graph(samples), eps))
  cluster, sparse = spectral(unweighted_knn(samples, 10))
  print(cluster, len(cluster), sparse)
  plot_2d_point_clusters(samples, cluster, "spectral_rectangles")


# cluster_sqrttrough_MCS_and_normal()


# Sample bowtie
def cluster_bowtie_MCS_and_normal():
  length = 1
  neck = 1 / 1000
  # Actually if depth is too large, since epsilon is big, the "augmented" cut
  # method shouldn't cut the trough.. because lots of high length edges cross
  # it.
  depth = 5
  n = 3000
  # Required for DFS to not break Python's recurison limit.
  # Making this not have +1000 made this break on small data.
  sys.setrecursionlimit(int(1.5 * n + 1000))

  eps = 2.7 * (2 * length**2 * depth / n) ** (1 / 3)
  samples = sample_bowtie(n, neck, length, depth)
  print("Samples", samples)
  print("eps", eps)
  cluster, sparse = spectral_points_bowtie_MCS(samples, neck, eps)
  truncated_samples = np.array([[sample[0], sample[1]] for sample in samples])
  print(cluster, len(cluster), sparse)
  # plot_points(truncated_samples, "bowtie.png")
  # plot_2d_point_clusters(truncated_samples, cluster, "spectral_R_bowtie.png")
  plot_3d_points(samples, "bowtie_samples.png")
  plot_3d_point_clusters(samples, cluster, "spectral_R_bowtie.png")
  cluster, sparse = spectral(eps_graph(affinity_graph(samples), eps))
  print(cluster, len(cluster), sparse)
  plot_3d_point_clusters(samples, cluster, "spectral_normal_bowtie.png")
  # plot_2d_point_clusters(
  #    truncated_samples, cluster, "spectral_normal_bowtie.png"
  # )


plot_densitybridges()
# plot_two_ellipses()
# cluster_bowtie_MCS_and_normal()
