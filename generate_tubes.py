#!/usr/bin/env python

from math import asin, acos, atan, cos, exp, log, pi, sin, sqrt
import random
import argparse
import sys
from tube_generator import *

def main():
  start_pos = (0, 0, 0)
  parser = argparse.ArgumentParser(description='Generate/draw cell microtubes.')
  parser.add_argument('-p', dest='print_cell', action='store_true', default=False,
                      help='print out one generated cell')
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
  draw_cell = not (args.calc_fl or args.calc_int or args.print_cell)
  if args.calc_int:
    args.calc_fl = False
    args.x = 20
  fl_calcer = None
  if args.calc_fl:
    fl_calcer = TFluorCalcer(args.cell_radius, args.cell_radius / 200.0)
  int_out = None
  if args.calc_int:
    int_out = [ open('data/int0.data', 'w'), \
                open('data/int1.data', 'w'), \
                open('data/int2.data', 'w')\
              ]
  for iteration in range(args.x):
    int_calcer = None
    cell_c = args.c
    if args.calc_int:
      #cell.c = iteration + 1.0) / 2
      cell_c = 0.1 + 2**iteration
      int_calcer = TIntersectionCalcer(args.cell_radius, args.tube_radius, args.tube_length / 10)
    cell = TCell(n=args.n, m=args.m, c=args.c, psi0=args.psi0, \
                 kernel_radius=args.kernel_radius, tube_radius=args.tube_radius, \
                 cell_radius=args.cell_radius, cell_height=args.cell_height, \
                 tube_length=args.tube_length, generator_cls='power')
    if args.print_cell:
      print cell
    if args.calc_fl:
      fl_calcer.add_cell(cell)
    if args.calc_int:
      int_calcer.add_cell(cell)
      int_out[0].write('%f %f\n' % (log(cell.c, 2), int_calcer.intersections[0]))
      int_out[1].write('%f %f\n' % (log(cell.c, 2), int_calcer.intersections[1]))
      int_out[2].write('%f %f\n' % (log(cell.c, 2), int_calcer.intersections[2]))
  if args.calc_fl:
    for i, x in enumerate(fl_calcer.get_fl()):
      print (i + 1) * fl_calcer.step / 10, x[0]#, x[1]
  if draw_cell:
    cell.draw()
    sys.stderr.write('DONE DRAWING\n')

if __name__ == '__main__':
  main()
