#!/usr/bin/env python

from math import asin, acos, atan, cos, exp, log, pi, sin, sqrt
import random

def r(point):
  return sqrt(sum([ c**2 for c in point ]))

#  *******************************  CALCULATORS  *******************************

class TFluorCalcer:
  def __init__(self, r, step):
    self.r = r
    self.step = step
    self.n = int(r / step) - 1
    self.__data = [0.0] * self.n

  def add(self, points):
    assert(len(points) == 2)
    rads = [ sqrt(sum([ c**2 for c in p ])) for p in points ]
    rads.sort()
    rads = [ int(r / self.step + (1 if r != r / self.step else 0)) for r in rads ]
    c = rads[0]
    for i in range(max((1, rads[0])), min((rads[1], self.n + 1))):
      self.__data[i - 1] += 1

  def get_fl(self):
    return [ x / (2 * pi * self.step * (i + 1)) for i, x in enumerate(self.__data) ]

class TIntersectionCalcer:
  def __init__(self, cell_radius, tube_radius, step):
    self.tube_radius = tube_radius
    self.cell_radius = cell_radius
    self.step = step
    self.intersections = [0] * 3
    self.n = int(cell_radius / step)
    self.__intersections = [ [0] * int(2 * cell_radius / 10.0) for i in range(int(2 * cell_radius / 10.0)) ]
    self.__data = [ [] for i in range(self.n) ]

  def intersects(self, line, i, j):
    EPS = 1e-4

    line2 = self.__data[i][j]
    if max((line[0][0], line[1][0])) < min((line2[0][0], line2[1][0])):
      return
    a1 = 0
    a2 = 0
    try:
      a1 = (line[1][1] - line[0][1]) / (line[1][0] - line[0][0])
      a2 = (line2[1][1] - line2[0][1]) / (line2[1][0] - line2[0][0])
    except:
      return
    b1 = line[0][1] - a1 * line[0][0]
    b2 = line2[0][1] - a2 * line2[0][0]
    if abs(a1 - a2) < EPS:
      if abs(b1 - b2) < EPS:
        return
      else:
        return
    x = (b2 - b1) / (a1 - a2)
    y = a1 * x + b1
    if r((x, y)) < i * self.step or r((x, y)) >= (i + 1) * self.step:
      return
    if x < max(( min((line[0][0], line[1][0])), min((line2[0][0], line2[1][0])))) or x > min((max((line[0][0], line[1][0])), max((line2[0][0], line2[1][0])))):
      return
    else:
      k = int((x + self.cell_radius) / 10)
      l = int((y + self.cell_radius) / 10)
      if self.__intersections[k][l] == 0:
        self.__intersections[k][l] = 1
        self.intersections[0] += 1
        self.intersections[1 if r((x, y)) < self.cell_radius / 2 else 2] += 1

  def add(self, line):
    assert(len(line) == 2)
    rads = [ sqrt(sum([ c**2 for c in point ])) for point in line ]
    rads.sort()
    rads = [ int(r / self.step) for r in rads ]
    rads[1] += 1
    for i in range(*rads):
      for j in range(len(self.__data[i])):
        self.intersects(line, i, j)
      self.__data[i].append(line)

#  *******************************  SECTIONS GENERATORS  *******************************

class TTubeSectionsGenerator:
  def __init__(self, cell, pos = (0, 0, 0)):
    self.cell = cell
    self.validate()
    self.direction = (0, 0)
    self.position = pos
    self.sections_built = 0
    self.__t = (0, 0, 0)
    self.__old_position = self.position

  def validate(self):
    try:
      assert(self.cell.psi0 >= -0.5*pi and self.cell.psi0 <= 0.5*pi)
    except AttributeError:
      pass

  """
    def __setattr__(self, name, value):
      self.__dict__[name] = value
      good_len = 0
      if name == 'direction':
        good_len = 2
      elif name == 'position':
        good_len = 3
      else:
        return
      assert(type(self.__getattr__(name)) == tuple and len(self.direction) == good_len)
  """

  def get_next_part(self):
    phi = 2 * pi * random.random()
    psi = 0
    if tuple(self.direction) == (0, 0):
      max_sin = sin(self.cell.cell_height / sqrt(self.cell.cell_height**2 + self.cell.cell_radius**2))
      psi = asin(2 * max_sin * random.random() - max_sin)
    else:
      psi = self.approx_inverse_value(random.random(), 0.0001)
    r = cos(psi) * self.cell.tube_length
    x2 = [ cos(phi) * r, sin(phi) * r, sin(psi) * self.cell.tube_length ]
    x1 = [0]*3
    x1[0] = x2[2] * sin(0.5 * pi - self.direction[1]) + x2[0] * cos(0.5 * pi - self.direction[1])
    x1[1] = x2[1]
    x1[2] = x2[2] * cos(0.5 * pi - self.direction[1]) - x2[0] * sin(0.5 * pi - self.direction[1])
    x0 = [0]*3
    x0[0] = x1[0] * cos(self.direction[0]) - x1[1] * sin(self.direction[0])
    x0[1] = x1[0] * sin(self.direction[0]) + x1[1] * cos(self.direction[0])
    x0[2] = x1[2]
    return x0

  def __iter__(self):
    return self

  def next(self):
    assert(self.sections_built <= self.cell.m)
    if self.sections_built == self.cell.m:
      raise StopIteration
    self.sections_built += 1
    addition = self.get_next_part()
    self.__old_position = self.position
    self.position = tuple([ self.__old_position[i] + addition[i] for i in range(0, 3) ])
    if abs(self.position[0]) >= self.cell.cell_height / 2:
      projected_length = sqrt(self.__t[1]**2 + self.__t[2]**2)
      if projected_length < 1:
        raise StopIteration
      multiplier = self.cell.tube_length / projected_length
      addition[0] = -50 if self.position[0] > 0 else 50
      addition[1] = self.__t[1] * multiplier
      addition[2] = self.__t[2] * multiplier
      self.position = tuple([ self.__old_position[i] + addition[i] for i in range(0, 3) ])
    if sqrt(self.position[1]**2 + self.position[2]**2) >= self.cell.cell_radius:
      raise StopIteration

    # Updating all auxiliary fields
    self.__t = tuple([ self.position[i] - self.__old_position[i] for i in range(3) ])
    new_direction = [0, 0]
    if self.__t[0] != 0:
      new_direction[0] = atan(self.__t[1] / self.__t[0])
      if self.__t[0] < 0:
        new_direction[0] += pi
    if self.__t[0] * self.__t[1] != 0:
      new_direction[1] = atan(self.__t[2] / sqrt(self.__t[0]**2 + self.__t[1]**2))
    self.direction = tuple(new_direction)
    return (self.__old_position, self.position)

  def approx_inverse_value(self, p, eps):
    coords = [ self.cell.psi0, 0.5 * pi ]
    while True:
      mean = 0.5 * sum(coords)
      if self.cumulative_distribution(coords[1]) - \
        self.cumulative_distribution(coords[0]) < eps:
        return mean
      if self.cumulative_distribution(mean) > p:
        coords[1] = mean
      else:
        coords[0] = mean

  def cumulative_distribution(self, psi):
    return (psi - self.cell.psi0) / (0.5*pi - self.cell.psi0)

class TPowerTubeSectionsGenerator(TTubeSectionsGenerator):
  def validate(self):
    TTubeSectionsGenerator.validate(self)
    try:
      assert(self.cell.c > 0)
    except AttributeError:
      pass

  def cumulative_distribution(self, psi):
    x = (psi - self.cell.psi0) / (0.5*pi - self.cell.psi0)
    x = x**(self.cell.c)
    return x

class TExpTubeSectionsGenerator(TTubeSectionsGenerator):
  def validate(self):
    TTubeSectionsGenerator.validate(self)
    try:
      assert(self.cell.c > 1)
    except AttributeError:
      pass

  def cumulative_distribution(self, psi):
    x = (self.cell.c**((psi - self.cell.psi0) / (0.5*pi - self.cell.psi0)) - 1) / (self.cell.c - 1)
    return x

def construct_generator(cell, generator):
  if generator == 'exponential':
    cls = TExpTubeSectionsGenerator
  elif generator == 'power':
    cls = TPowerTubeSectionsGenerator
  elif generator == uniform:
    cls = TTubeSectionsGenerator
  return cls(cell)

#  *******************************  CELLS  *******************************

class TCell:
  def __init__(self, n, m, c, psi0, kernel_radius, tube_radius, cell_radius, \
               cell_height, tube_length, generator_type):
    self.n = n
    self.m = m
    self.c = c
    self.psi0 = psi0
    self.kernel_radius = kernel_radius
    self.tube_radius = tube_radius
    self.cell_radius = cell_radius
    self.cell_height = cell_height
    self.tube_length = tube_length
    self.generator_type = generator_type
    self.validate()
    self.start()

  def validate(self):
    try:
      assert(self.kernel_radius <= self.cell_radius)
    except AttributeError:
      pass
    try:
      self.generator.validate()
    except AttributeError:
      pass

  def start(self):
    self.tubes_built = 0

  def __setattr__(self, name, value):
    self.__dict__[name] = value
    self.validate()

  def __iter__(self):
    return self

  def next(self):
    assert(self.tubes_built <= self.n)
    if self.tubes_built == self.n:
      raise StopIteration
    self.tubes_built += 1
    p = (0, 0, 0)
    if self.kernel_radius != 0:
      radius = self.kernel_radius * sqrt(random.random())
      phi = 2 * pi * random.random()
      h = 0.9 * self.cell_height * random.random() - 0.45 * self.cell_height
      p = (h, radius * sin(phi), radius * cos(phi))
    self.generator = construct_generator(self, self.generator_type)
    return self.generator

  def info(self):
    out  = str(self.cell_height) + '\t' + str(self.cell_radius) + '\n'
    out += str(self.tube_length) + '\t' + str(self.tube_radius) + '\n'
    return out

"""
def generate_tube_part_uniform(start, length):
  phi = 2 * pi * random.random()
  psi = asin(2 * random.random() - 1)
  r = cos(psi) * length
  delta = [ cos(phi) * r, sin(phi) * r, sin(psi) * length ]
  return [ start[i] + delta[i] for i in range(0, 3) ]

def cumulative_distribution(c, psi0, psi):
  #x = (c**((psi - psi0) / (0.5*pi - psi0)) - 1) / (c - 1)
  #return x
  x = (psi - psi0) / (0.5*pi - psi0)
  x = x**c
  return x
  x = (0.25 * log(2))**(1/c)
  x /= 0.5 * pi - psi0
  x *= psi
  x -= psi0
  return exp(4 * (x**c)) - 1

def approx_inverse_value(p, c, psi0, eps):
  coords = [ psi0, 0.5 * pi ]
  while True:
    mean = 0.5 * sum(coords)
    if cumulative_distribution(c, psi0, coords[1]) - \
      cumulative_distribution(c, psi0, coords[0]) < eps:
      return mean
    if cumulative_distribution(c, psi0, mean) > p:
      coords[1] = mean
    else:
      coords[0] = mean

def generate_tube_part_linear(length, c, psi0, direction, max_sin=0):
  phi = 2 * pi * random.random()
  psi = 0
  if max_sin != 0:
    psi = asin(2 * max_sin * random.random() - max_sin)
  else:
    psi = approx_inverse_value(random.random(), c, psi0, 0.0001)
  r = cos(psi) * length
  x2 = [ cos(phi) * r, sin(phi) * r, sin(psi) * length ]
  x1 = [0]*3
  x1[0] = x2[2] * sin(0.5 * pi - direction[1]) + x2[0] * cos(0.5 * pi - direction[1])
  x1[1] = x2[1]
  x1[2] = x2[2] * cos(0.5 * pi - direction[1]) - x2[0] * sin(0.5 * pi - direction[1])
  x0 = [0]*3
  x0[0] = x1[0] * cos(direction[0]) - x1[1] * sin(direction[0])
  x0[1] = x1[0] * sin(direction[0]) + x1[1] * cos(direction[0])
  x0[2] = x1[2]
  return x0
"""
