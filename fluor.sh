#!/bin/bash -v

args='-f -x 2 -k 1000 -n 480 -m 150 -R5000'

./generate_tubes.py -c 1   $args > data/1.data
./generate_tubes.py -c 1.1 $args > data/1.5.data
./generate_tubes.py -c 1.2   $args > data/2.data
./generate_tubes.py -c 4   $args > data/4.data
./generate_tubes.py -c 10  $args > data/10.data
./plot.gplot
