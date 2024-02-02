from pylab import *


# Bowtie rect seems to cut properly with spectral!
def bowtie_rect(rect_x, rect_y, bowtie_x, bowtie_neck, eps):
  x = -bowtie_x / 2 - rect_x
  coords = []
  while x <= bowtie_x / 2 + rect_x:
    y = -rect_y / 2
    while y <= rect_y / 2:
      if (
          x <= -bowtie_x / 2
          or x >= bowtie_x / 2
          or (y <= x + bowtie_neck / 2 and -y <= x + bowtie_neck / 2)
          or (y <= -x + bowtie_neck / 2 and -y <= -x + bowtie_neck / 2)
      ):
        coords.append([x, y])
      y += eps
    x += eps
  return array(sorted(coords))


def densitybridges(max_x, max_y, bridge_x, bridge_y, bridge_spacing, n):
  samples = 0
  coords = []
  while samples < n:
    x = (rand() - 1 / 2) * max_x
    y = (rand() - 1 / 2) * max_y
    q = y % (bridge_spacing + bridge_y)
    if (
        -bridge_x / 2 >= x
        or x >= bridge_x / 2
        or (q < bridge_y / 2 or q > bridge_spacing + bridge_y / 2)
    ):
      coords.append([x, y])
      samples += 1
  return array(sorted(coords))


# Use plot_rectbridges to see what this does. It's basically a bunch of bridges
# across two rectangles.
def rectbridges(max_x, max_y, bridge_x, bridge_y, bridge_spacing, eps):
  x = -max_x / 2
  coords = []
  while x <= max_x / 2:
    y = -max_y / 2
    while y <= max_y / 2:
      q = y % (bridge_spacing + bridge_y)
      if (
          -bridge_x / 2 >= x
          or x >= bridge_x / 2
          or (q < bridge_y / 2 or q > bridge_spacing + bridge_y / 2)
      ):
        coords.append([x, y])
      y += eps
    x += eps
  return array(sorted(coords))


# Gaps of eps
# Run plot_rectbridge to see what this does.
def rectbridge(max_x, max_y, bridge_x, bridge_y, eps):
  x = -max_x / 2
  coords = []
  while x <= max_x / 2:
    y = -max_y / 2
    while y <= max_y / 2:
      if not (
          -bridge_x / 2 <= x <= bridge_x / 2
          and (y > bridge_y / 2 or y < -bridge_y / 2)
      ):
        coords.append([x, y])
      y += eps
    x += eps
  # We Sort this because when we output clusters, this lets us figure out if the
  # clusters are correct.
  return array(sorted(coords))


def single_sample_sqrttrough(valley_height, depth):
  y = rand() * depth
  x = 2 * rand() - 1
  dummy = rand() * (1 + valley_height)
  while dummy > abs(x + valley_height**2) ** (1 / 2):
    x = 2 * rand() - 1
    dummy = rand() * (1 + valley_height)
  return [x, y]


def sample_rectangle(n, center_x, center_y, length_x, length_y):
  coords = []
  for i in range(n):
    x = length_x * (rand() - 0.5)
    y = length_y * (rand() - 0.5)
    coords.append([x + center_x, y + center_y])
  return array(sorted(coords))


def sample_ellipse(n, center_x, center_y, radius_x, radius_y):
  coords = []
  for i in range(n):
    x = radius_x
    y = radius_y
    while x**2 / (radius_x) ** 2 + y**2 / (radius_y) ** 2 >= 1:
      x = radius_x * (2 * rand() - 1)
      y = radius_y * (2 * rand() - 1)
    coords.append([x + center_x, y + center_y])
  return array(sorted(coords))


def sample_sqrttrough(n, valley_height, depth):
  coords = sorted(
      [single_sample_sqrttrough(valley_height, depth) for _ in range(n)]
  )
  return array(coords)


def sample_mixture_two_gaussians_2d(n, mean1, cov1, mean2, cov2):
  coords = []
  print(multivariate_normal(mean1, cov1))
  for i in range(n):
    coin = rand() < 0.5
    if coin:
      print("Mean 1 achieved")
      coords.append(list(multivariate_normal(mean1, cov1)))
    else:
      print("Mean 2 achieved")
      coords.append(list(multivariate_normal(mean2, cov2)))
  return array(sorted(coords))


# Bowtie: |y| - neck <= |x| from [0] to [length], [0] to [depth]
# Returns a sample or nothing.
def single_sample_bowtie(neck, length, depth):
  x = 2 * rand() * length - length
  y = 2 * rand() * (length + neck) - (length + neck)
  while abs(y) - neck > abs(x):
    x = 2 * rand() * length - length
    y = 2 * rand() * (length + neck) - (length + neck)
  z = rand() * depth
  return [x, y, z]


def sample_bowtie(n, neck, length, depth):
  coords = sorted([single_sample_bowtie(neck, length, depth) for _ in range(n)])
  return array(coords)


def norm(point):
  val = 0
  for i in point:
    val += i**2
  return sqrt(val)
