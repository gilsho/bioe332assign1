import os
from matplotlib.pyplot import *
from numpy import *
from collections import namedtuple

SIMTIME = 6
BASE_FILE_NAME = 'sim'
FAMILY = 2048
GROUPS = 3
NUM_SAMPLES		 = 30000-1
NUM_FILES = 30

d = dict()
d[0] = 'sim' + str(FAMILY) + '/'
d[1] = 'sim' + str(FAMILY) + 'Leak1/'
d[2] = 'sim' + str(FAMILY) + 'Leak2/'
d[3] = 'sim' + str(FAMILY) + 'Leak3/'

t = arange(0,SIMTIME,1.*SIMTIME/NUM_SAMPLES)	


def process_leak_group(group):
	print 'processing leak group' + str(group) + '...'
	variance = zeros(NUM_SAMPLES)
	for nsim in range(1,NUM_FILES):
		drift = process_sim_file(group,nsim)
		variance += array(drift)**2								#add current simulation drifts to avg
	variance /= NUM_FILES
	return variance
	

def process_sim_file(group,nsim):
	drift = []
	print 'processing file: ' + d[group] + str(FAMILY) + BASE_FILE_NAME + str(nsim) + '.dat ...'
	f = open(d[group] + str(FAMILY) + BASE_FILE_NAME + str(nsim) + '.dat','r')
	ang_cue = float(f.readline())
	for line in f:
		drift.append(ang_cue - float(line))
	f.close()
	figure(group)
	plot(t,drift,'k')
	return drift


for i in range(3):
	variance = process_leak_group(i)
	figure(0)
	plot(t,variance,label='Leak variance ' + str(i) + 'mV')
ylabel('Population Vector Variance')
xlabel('Time (sec)')
legend(loc=2)
savefig('hetro')
show()



