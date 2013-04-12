import os
from matplotlib.pyplot import *
from numpy import *

SIMTIME = 6
FAMILIES = [1024,2048,4096]
NUM_FILES = 26
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
		print 'processing simulation number ' + str(nsim) + '...'
		drift = process_sim_file(fam,nsim)
		variance += array(drift)**2								#add current simulation drifts to avg
		plot(t,drift,'k')
	variance /= NUM_FILES
	savefig('trace' + str(fam))
	return variance
	

def process_sim_file(ndir,nsim):
	drift = []
	f = open('sim' + str(ndir) + '/' + str(ndir) + BASE_FILE_NAME + str(nsim) + '.dat','r')
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
savefig('variance')
show()



