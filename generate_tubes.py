#!/usr/bin/env python

from math import asin, acos, atan, cos, exp, log, pi, sin, sqrt
import random
import argparse
import sys
import os.path
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
  parser.add_argument('-M', '--match', dest='match', action='store', type=str, default='',
                      help='match parameters with experimental data from file')
  args = parser.parse_args()
  draw_cell = not (args.calc_fl or args.calc_int or args.print_cell or args.match)
  if args.calc_int:
    int_out = [ open('data/int0.data', 'w'), \
                open('data/int1.data', 'w'), \
                open('data/int2.data', 'w')\
              ]
    for iteration in range(20):   # Number of experiments
      int_calcer = TIntersectionCalcer(args.cell_radius, args.tube_radius, args.tube_length / 10)
      cell_c = (iteration + 1.0) / 2
      cell = TCell(vars(args), 'power', c=cell_c)
      int_calcer.add_cell(cell)
      int_out[0].write('%f %f\n' % (log(cell.c, 2), int_calcer.intersections[0]))
      int_out[1].write('%f %f\n' % (log(cell.c, 2), int_calcer.intersections[1]))
      int_out[2].write('%f %f\n' % (log(cell.c, 2), int_calcer.intersections[2]))
  elif args.calc_fl:
    fl_calcer = TFluorCalcer(args.cell_radius, args.cell_radius / 200.0)
    for iteration in range(args.x):
      cell = TCell(vars(args), 'power')
      fl_calcer.add_cell(cell)
    for i, x in enumerate(fl_calcer.get_fl()):
      print (i + 1) * fl_calcer.step / 10, x[0]#, x[1]
  else:
    cell = TCell(vars(args), 'power')
    if args.print_cell:
      print cell
    else:
      cell.draw()
      sys.stderr.write('DONE DRAWING\n')

if __name__ == '__main__':
  main()
