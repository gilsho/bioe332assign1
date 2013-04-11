import os
from matplotlib.pyplot import *
from numpy import *

SIMTIME = 6
FAMILIES = [4096]
NUM_FILES = 1
BASE_FILE_NAME = 'sim'
NUM_SAMPLES		 = 30000-1

t = arange(0,SIMTIME,1.*SIMTIME/NUM_SAMPLES)	

def process_family(fam):
	variance = zeros(NUM_SAMPLES)
	figure(fam)
	title('NE='+str(fam))
	ylabel('Population Vector')
	xlabel('Time (sec)')
	axis([0,SIMTIME,-40,40])
	for nsim in range(1,NUM_FILES+1):
		drift = process_sim_file(fam,nsim)
		variance += array(drift)**2								#add current simulation drifts to avg
		plot(t,drift,'k')
	variance /= NUM_FILES
	return variance
	

def process_sim_file(ndir,nsim):
	drift = []
	f = open(str(ndir) + BASE_FILE_NAME + str(nsim) + '.dat','r')
	ang_cue = float(f.readline())
	for line in f:
		drift.append(ang_cue - float(line))
	f.close()
	return drift


for fam in FAMILIES:
	variance = process_family(fam)
	figure(0)
	plot(t,variance,label=str(fam))
ylabel('Population Vector Variance')
xlabel('Time (sec)')
legend(loc=2)
show()



