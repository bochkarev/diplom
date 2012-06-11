#!/usr/bin/env python

from math import asin, acos, atan, cos, exp, log, pi, sin, sqrt
import random
import argparse
import sys
from tube_generator import *

if __name__ == '__main__':
  start_pos = (0, 0, 0)
  parser = argparse.ArgumentParser(description='Generate/draw cell microtubes.')
  parser.add_argument('-d', dest='draw', action='store_true', default=False,
                      help='draw one generated cell')
  parser.add_argument('-f', dest='calc_fl', action='store_true',
                      default=False, help='calculate fluorescence decline mode')
  parser.add_argument('-i', dest='calc_int', action='store_true',
                      default=False, help='calculate intersections mode')
  parser.add_argument('-k', dest='kernel_radius', action='store', default=0,
                      type=float, help='cell kernel radius')
  parser.add_argument('-H', '--cell-height', dest='cell_height',
                      action='store', default=2000, type=float,
                      help='cell height')
  parser.add_argument('-R', '--cell-radius', dest='cell_radius',
                      action='store', default=10000, type=float,
                      help='cell radius')
  parser.add_argument('-l', '--tube-length', dest='tube_length',
                      action='store', default=500, type=float,
                      help='tube length')
  parser.add_argument('-r', '--tube-radius', dest='tube_radius',
                      action='store', default=50, type=float,
                      help='tube radius')
  parser.add_argument('-n', dest='n', action='store', default=300, type=int,
                      help='number of microtubes')
  parser.add_argument('-m', dest='m', action='store', default=50, type=int,
                      help='number of parts in one microtube')
  parser.add_argument('-c', dest='c', action='store', default=1, type=float,
                      help='c')
  parser.add_argument('-z', '--psi-zero', dest='psi0', action='store',
                      default=0, type=float,
                      help='starting value for microtube psi angle')
  parser.add_argument('-x', dest='x', action='store', default=1, type=int,
                      help='number of experiments')
  args = parser.parse_args()
  if args.draw:
    from visual import *
    cell_height, cell_radius = [int(x) for x in
                                sys.stdin.readline().split('\t')]
    tube_length, tube_radius = [int(x) for x in
                                sys.stdin.readline().split('\t')]
    scene.background = (0.6, 0.6, 0.6)
    """cylinder(pos=(0, 0, 0),
             axis=(1000, 0, 0),
             radius=2,
             color=color.red)
    cylinder(pos=(0, 0, 0),
             axis=(0, 1000, 0),
             radius=2,
             color=color.green)
    cylinder(pos=(0, 0, 0),
             axis=(0, 0, 1000),
             radius=2,
             color=color.blue)"""
    cell = cylinder(pos=(0, 0, -0.5 * cell_height),
                    axis=(0, 0, cell_height),
                    radius=cell_radius,
                    color=color.blue,
                    opacity=0.25)
    sphere(pos=(0, 0, 0), radius=50, color=color.green)
    sys.stdin.readline()
    line = sys.stdin.readline()
    p1 = [ float(x) for x in line.split('\t') ]
    while line:
      if line == '\n':
        line = sys.stdin.readline()
        p1 = [ float(x) for x in line.split('\t') ]
        line = sys.stdin.readline()
        continue
      p2 = [float(x) for x in line.split('\t')]
      if len(p2) == 2:
        line = sys.stdin.readline()
        continue
      direction = [p2[i] - p1[i] for i in range(3)]
      cylinder(pos=p1,
               axis=direction,
               radius=tube_radius,
               color=color.red)
      p1 = p2
      line = sys.stdin.readline()
    sys.stderr.write('DONE DRAWING\n')
    exit()
  calc_fl = args.calc_fl
  calc_int = args.calc_int
  print_tubes = not (calc_fl or calc_int)
  x = args.x
  if calc_int:
    calc_fl = False
    x = 20
  fl_calcer = None
  if calc_fl:
    fl_calcer = TFluorCalcer(args.cell_radius, args.cell_radius / 200.0)
  int_out = None
  if calc_int:
    int_out = [ open('data/int0.data', 'w'), \
                open('data/int1.data', 'w'), \
                open('data/int2.data', 'w')\
              ]
  for iteration in range(x):
    int_calcer = None
    cell_c = args.c
    if calc_int:
      #cell.c = iteration + 1.0) / 2
      cell_c = 0.1 + 2**iteration
      int_calcer = TIntersectionCalcer(cell.cell_radius, cell.tube_radius, cell.tube_length / 10)
    cell = TCell(n = args.n, m = args.m, c = args.c, psi0 = args.psi0, \
                 kernel_radius = args.kernel_radius, tube_radius = args.tube_radius, \
                 cell_radius = args.cell_radius, cell_height = args.cell_height, \
                 tube_length = args.tube_length, generator_cls = 'power')
    if print_tubes:
      sys.stdout.write(cell.info())
    for tube in cell:
      if print_tubes:
        print
      first_pass = True
      for old_position, p in tube:
        if print_tubes:
          if first_pass:
            print '\t'.join([ str(x) for x in [old_position[2], old_position[1], old_position[0]]])
          print '\t'.join([ str(x) for x in [p[2], p[1], p[0]]])
        if calc_fl:
          fl_calcer.add((old_position[1:3], p[1:3]))
        if calc_int:
          int_calcer.add((old_position[1:3], p[1:3]))
    if calc_int:
      int_out[0].write('%f %f\n' % (log(cell.c, 2), int_calcer.intersections[0]))
      int_out[1].write('%f %f\n' % (log(cell.c, 2), int_calcer.intersections[1]))
      int_out[2].write('%f %f\n' % (log(cell.c, 2), int_calcer.intersections[2]))
  if calc_fl:
    for i, x in enumerate(fl_calcer.get_fl()):
      print (i + 1) * fl_calcer.step / 10, x * 1000
