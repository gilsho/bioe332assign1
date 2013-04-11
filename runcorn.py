#!/usr/bin/python

import os
import datetime 

import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-Np', type = int, default = 2048)
parser.add_argument('-nsim', type = int, default = 100)
parser.add_argument('-rleak', type = int, default=0)
pa = parser.parse_args()

for sim in xrange(1, pa.nsim+1): 
	if (pa.rleak):
		job = 'python ~/Desktop/Classes/BIOE332/assign1/assign1.py -seed %d -Np %d -rleak 1' % (sim, pa.Np)
	else:
		job = 'python ~/Desktop/Classes/BIOE332/assign1/assign1.py -seed %d -Np %d' % (sim, pa.Np)
	os.system(job)

