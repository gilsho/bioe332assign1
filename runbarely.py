#!/usr/bin/python

import os
import datetime 

families = [1024,2048,4096]
simulations = 100


startjob = datetime.datetime.now() 

for fam in families:
	for sim in xrange(1, simulations+1): 
		job = 'qsub -N python ~/Desktop/Classes/BIOE332/assign1/assign1.py -seed %d -Np %d' % (sim, fam)
		print job
		os.system(job)

endjob = datetime.datetime.now()

print 'time to submit all jobs: ', endjob-startjob
