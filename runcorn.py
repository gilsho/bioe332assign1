#!/usr/bin/python

import os
import datetime 

import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-Np', type = int, default = 2048)
parser.add_argument('-nsim', type = int, default = 100)
parser.add_argument('-rleak', type = int, default=0)
parser.add_argument('-nstart', type = int, default=1)
parser.add_argument('-decay', type = float, default=0.99)
pa = parser.parse_args()

for sim in xrange(pa.nstart, pa.nsim+1): 
	job = 'python ~/Desktop/Classes/BIOE332/bioe332assign1/assign1.py -seed %d -Np %d -rleak %d -decay %d' % (sim, pa.Np, pa.rleak, pa.decay)
	os.system(job)

