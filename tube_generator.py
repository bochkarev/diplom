#!/usr/bin/env python

from math import asin, acos, atan, cos, exp, log, pi, sin, sqrt
import random
import sys
import time

def r(point):
  return sqrt(sum([ c**2 for c in point ]))

class TInterpolatedFunction:
  def __init__(self):
    self.__data = []
    self.sorted = True

  def __setitem__(self, key, value):
    self.__data.append((float(key), float(value)))
    self.sorted = False

  def sort(self):
    self.__data.sort()
    self.sorted = True

  def __getitem__(self, key):
    if len(self.__data) < 2:
      raise IndexError
    if not self.sorted:
      self.sort()
    beg = 0
    end = len(self.__data) - 1
    if key < self.__data[beg][0] or key > self.__data[end][0]:
      raise IndexError
    while end - beg > 1:
      mid = (end + beg) / 2
      if key < self.__data[mid][0]:
        end = mid
      elif key > self.__data[mid][0]:
        beg = mid
      else:
        return self.__data[mid][1]
    assert(key >= self.__data[beg][0] and key <= self.__data[end][0])
    alpha = (key - self.__data[beg][0]) / (self.__data[end][0] - self.__data[beg][0])
    assert(alpha >= 0 and alpha <= 1)
    return alpha * self.__data[beg][1] + (1 - alpha) * self.__data[end][1]

#  *******************************  CALCULATORS  *******************************

class TFluorCalcer:
  def __init__(self, r, step):
    self.r = r
    self.step = step
    self.__data = [ list() for i in range(int(r / step) - 1) ]

  def add_cell(self, cell):
    for measurements in self.__data:
      measurements.append(0.0)
    for tube in cell:
      for i in range(len(tube) - 1):
        self.add_section([ x[1:3] for x in tube[i:i+2] ])
    for idx, measurements in enumerate(self.__data):
      measurements[-1] /= 2 * pi * self.step * (idx + 1)

  def add_section(self, points):
    assert(len(points) == 2)
    rads = [ sqrt(sum([ c**2 for c in p ])) for p in points ]
    rads.sort()
    rads = [ int(r / self.step + (1 if r != r / self.step else 0)) for r in rads ]
    c = rads[0]
    for i in range(max((1, rads[0])), min((rads[1], len(self.__data) + 1))):
      self.__data[i - 1][-1] += 1000

  def get_fl(self):
    res = [ [0.0] * 2 for i in range(len(self.__data)) ]
    for idx, measurements in enumerate(self.__data):
      res[idx][0] = sum(measurements) / len(measurements)
      res[idx][1] = sum([ (m - res[idx][0])**2 for m in measurements ]) / len(measurements)
    return res

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
    if x < max((min((line[0][0], line[1][0])), min((line2[0][0], line2[1][0])))) or \
       x > min((max((line[0][0], line[1][0])), max((line2[0][0], line2[1][0])))):
      return
    else:
      k = int((x + self.cell_radius) / 10)
      l = int((y + self.cell_radius) / 10)
      if self.__intersections[k][l] == 0:
        self.__intersections[k][l] = 1
        self.intersections[0] += 1
        self.intersections[1 if r((x, y)) < self.cell_radius / 2 else 2] += 1

  def add_cell(self, cell):
    for tube in cell:
      for i in range(len(tube) - 1):
        self.add_section([ x[1:3] for x in tube[i:i+2] ])

  def add_section(self, line):
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

class TTube:
  def __init__(self, cell, pos=(0, 0, 0)):
    self.cell = cell
    self.direction = (0, 0)
    self.position = pos
    self.__sections = [ pos ]
    self.__t = (0, 0, 0)
    self.__old_position = self.position
    self.__iteration_stopped = False
    self.validate()
    self.build()

  def validate(self):
    try:
      assert(self.cell.psi0 >= -0.5*pi and self.cell.psi0 <= 0.5*pi)
    except AttributeError:
      pass

  def build(self):
    for section in self:
      pass

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
    if len(self.__sections) > 1:
      self.wait()
      return self.__sections.__iter__()
    else:
      return self

  def __len__(self):
    self.wait()
    return self.__sections.__len__()

  def __getitem__(self, key):
    return self.__sections.__getitem__(key)

  def __str__(self):
    self.wait()
    out = '\n'.join([ '\t'.join([ str(x) for x in [p[2], p[1], p[0]]]) for p in self ])
    return out

  def wait(self):
    while not self.__iteration_stopped:
      time.sleep(0.001)

  def next(self):
    assert(len(self.__sections) <= self.cell.m + 1)
    if len(self.__sections) == self.cell.m + 1:
      self.__iteration_stopped = True
      raise StopIteration
    addition = self.get_next_part()
    self.__old_position = self.position
    self.position = tuple([ self.__old_position[i] + addition[i] for i in range(0, 3) ])
    if abs(self.position[0]) >= self.cell.cell_height / 2:
      projected_length = sqrt(self.__t[1]**2 + self.__t[2]**2)
      if projected_length < 1:
        self.__iteration_stopped = True
        raise StopIteration
      multiplier = self.cell.tube_length / projected_length
      addition[0] = -50 if self.position[0] > 0 else 50
      addition[1] = self.__t[1] * multiplier
      addition[2] = self.__t[2] * multiplier
      self.position = tuple([ self.__old_position[i] + addition[i] for i in range(0, 3) ])
    if sqrt(self.position[1]**2 + self.position[2]**2) >= self.cell.cell_radius:
      self.__iteration_stopped = True
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
    self.__sections.append(self.position)
    return self.__sections[-1]

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

class TPowerTube(TTube):
  def validate(self):
    TTube.validate(self)
    try:
      assert(self.cell.c > 0)
    except AttributeError:
      pass

  def cumulative_distribution(self, psi):
    x = (psi - self.cell.psi0) / (0.5*pi - self.cell.psi0)
    x = x**(self.cell.c)
    return x

class TExpTube(TTube):
  def validate(self):
    TTube.validate(self)
    try:
      assert(self.cell.c > 1)
    except AttributeError:
      pass

  def cumulative_distribution(self, psi):
    x = (self.cell.c**((psi - self.cell.psi0) / (0.5*pi - self.cell.psi0)) - 1) / (self.cell.c - 1)
    return x

#  *******************************  CELLS  *******************************

class TCell:
  def __init__(self, args, generator_cls, **more_args):
    for s in ('n', 'm', 'c', 'psi0', 'kernel_radius', 'tube_radius', \
              'cell_radius', 'cell_height', 'tube_length'):
      if s in more_args:
        self.__setattr__(s, more_args[s])
      else:
        self.__setattr__(s, args[s])
    self.set_generator_cls(generator_cls)
    self.__tubes = list()
    self.__iteration_stopped = False
    self.validate()
    self.build()

  def validate(self):
    try:
      assert(self.kernel_radius <= self.cell_radius)
    except AttributeError:
      pass
    try:
      self.generator.validate()
    except AttributeError:
      pass

  def set_generator_cls(self, generator_cls):
    if generator_cls == 'exponential':
      cls = TExpTube
    elif generator_cls == 'power':
      cls = TPowerTube
    elif generator_cls == uniform:
      cls = TTube
    self.generator_cls = cls

  def __setattr__(self, name, value):
    self.__dict__[name] = value
    self.validate()

  def __iter__(self):
    if self.__tubes:
      self.wait()
      return self.__tubes.__iter__()
    else:
      return self

  def __len__(self):
    return self.__tubes.__len__()

  def __getitem__(self, key):
    return self.__tubes.__getitem__(key)

  def __str__(self):
    self.wait()
    out  = str(self.cell_height) + '\t' + str(self.cell_radius) + '\n'
    out += str(self.tube_length) + '\t' + str(self.tube_radius) + '\n'
    out += '\n'
    out += '\n\n'.join([ tube.__str__() for tube in self.__tubes ])
    return out

  def wait(self):
    while not self.__iteration_stopped:
      time.sleep(0.001)

  def next(self):
    assert(len(self.__tubes) <= self.n)
    if len(self.__tubes) == self.n:
      self.__iteration_stopped = True
      sys.stderr.write('DONE GENERATING\n')
      raise StopIteration
    p = (0, 0, 0)
    if self.kernel_radius != 0:
      radius = self.kernel_radius * sqrt(random.random())
      phi = 2 * pi * random.random()
      h = 0.9 * self.cell_height * random.random() - 0.45 * self.cell_height
      p = (h, radius * sin(phi), radius * cos(phi))
    self.__tubes.append(self.generator_cls(self, p))
    return self.__tubes[-1]

  def info(self):
    out  = str(self.cell_height) + '\t' + str(self.cell_radius) + '\n'
    out += str(self.tube_length) + '\t' + str(self.tube_radius) + '\n'
    return out

  def build(self):
    if not self.__iteration_stopped and not self.__tubes:
      for tube in self:
        pass

  def draw(self):
    import visual

    visual.scene.background = (0.6, 0.6, 0.6)
    """visual.cylinder(pos=(0, 0, 0),
             axis=(1000, 0, 0),
             radius=2,
             color=visual.color.red)
    visual.cylinder(pos=(0, 0, 0),
             axis=(0, 1000, 0),
             radius=2,
             color=visual.color.green)
    visual.cylinder(pos=(0, 0, 0),
             axis=(0, 0, 1000),
             radius=2,
             color=visual.color.blue)"""
    visual.cylinder(pos=(0, 0, -0.5 * self.cell_height),
                    axis=(0, 0, self.cell_height),
                    radius=self.cell_radius,
                    color=visual.color.blue,
                    opacity=0.25)
    #visual.sphere(pos=(0, 0, 0), radius=50, color=visual.color.green)
    for tube in self.__tubes:
      for i in range(len(tube) - 1):
        p1, p2 = [ (p[2], p[1], p[0]) for p in tube[i:i+2] ]
        direction = [p2[i] - p1[i] for i in range(3)]
        visual.cylinder(pos=p1,
                        axis=direction,
                        radius=self.tube_radius,
                        color=visual.color.red)
