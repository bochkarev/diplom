#!/bin/bash -v

args='-f -x 2 -k 2000 -n 100 -m 1000'

./generate_tubes.py -c 0.5 $args > data/0.5.data
./generate_tubes.py -c 1   $args > data/1.data
./generate_tubes.py -c 2.5 $args > data/2.5.data
./generate_tubes.py -c 5   $args > data/5.data
./generate_tubes.py -c 10  $args > data/10.data
./plot.gplot
