import brian_no_units
from brian import *
from scipy.special import erf
import numpy
from numpy.fft import rfft, irfft
from matplotlib.pyplot import *
import argparse


POPULATION_OUTPUT_FILE_BASE = 'sim'
NEURON_OUTPUT_FILE_BASE = 'neuron'

# parse arguments according to format:
# -N <excitatory population size>
# -s <random seed>
parser = argparse.ArgumentParser()
parser.add_argument('-Np', type = int, default = '2048')
parser.add_argument('-seed', type = int, default = 0)
parser.add_argument('-rleak', type = int, default = False)
parser.add_argument('-icue', type = double, default=0.000000001)
pa = parser.parse_args()
numpy.random.seed(pa.seed)

# parameters (pg 2 of Compte et al. 2000)

Np = pa.Np               # number of excitatory (pyramidal) neurons
Ni = Np/4                # number of inhibitory interneurons

dt_sim = 0.2*ms
sim_duration = 6.35 * second    # in seconds

#simulation_clock = Clock(dt = 0.02*ms)
current_clock = Clock(dt = dt_sim)
#dt_current = 0.1 * msecond
#dt_simulation = 0.02 * msecond
#current_clock = Clock(dt=dt_current)
#simulation_clock = Clock(dt=dt_simulation)

stim_start = 100 * msecond
stim_stop = 350 * msecond


# each type of cell is characterized by 6 intrinsic parameters:

# pyramidal cells
Cm_p = 0.5 * nF           # total capacitance Cm
gL_p = 25 * nS            # total leak conductance gL
El_p = -70 * mvolt        # leak reversal potential El
Vt_p = -50 * mvolt        # threshold potential Vt
Vr_p = -60 * mvolt        # reset potential Vr
tauRef_p = 2 * msecond    # refractory time tau

tau_p = gL_p/Cm_p         # membrane time constant


# interneurons
Cm_i = 0.2 * nF           # total capacitance Cm
gL_i = 20 * nS            # total leak conductance gL

# leak reversal potential El
if (pa.rleak):
    El_i = numpy.random.normal(-70,1) * mvolt
else:
    El_i = -70 * mvolt        

Vt_i = -50 * mvolt        # threshold potential Vt
Vr_i = -60 * mvolt        # reset potential Vr
tauRef_i = 1 * msecond    # refractory time tau

tau_i = gL_i/Cm_i         # membrane time constant


# external input is modeled as uncorrelated Poisson spike trains
# to each neuron at a rate of v_ext = 1800 Hz per cell
fext = 1800 * hertz

# equation 3 constants
a = 0.062 * 1/mvolt
b = 1/3.57

# external input mediated exclusively by AMPA receptors
g_ext_p = 3.1 * nS          # max conductance on pyramidal (excitatory) cells
g_ext_i = 2.38 * nS         # max conductance on interneurons (inhibitory)


# synaptic responses 
tau_syn_ampa = 2 * msecond     # time constants for AMPA
tau_syn_gaba = 10 * msecond    # time constants for GABA
tau_syn_nmda = 100 * msecond   # time constant for NMDA decay time
tau_syn_x = 2 * msecond        # time constant for NMDAR rise time
alfa = 500 * hertz             # controls saturation properties of NMDAR

# synaptic reversal potential
E_ampa = 0 * mvolt          # excitatory syn reversal potential
E_gaba = -70 * mvolt        # inhibitory syn reversal potential
E_nmda = 0 * mvolt          # excitatory syn reversal potential


# recurrent connection strengths (conductances)  (page 3 of Compte)
Gee = 0.381 * nS            # excitatory (pyramid-to-pyramid)
Gei = 0.292 * nS            # excitatory (pyramid-to-interneuron)
Gie = 1.336 * nS            # inhibitory (interneuron-to-pyramid)
Gii = 1.024 * nS            # inhibitory (interneuron-to-interneuron)

# connectivity profile and normalization condition
sigma_ee = 14.4     # degrees
Jp_ee = 1.62

# distribution of input current depending on orientation
i_cue_amp = pa.icue * amp  #  0.0000000001  is threshold
i_cue_ang = 180
i_cue_width = 10

sample_size = 5

sample_active_start = i_cue_ang-i_cue_width
sample_active_end = i_cue_ang+i_cue_width
sample_step = int((sample_active_end-sample_active_start)/sample_size)
sample_exc_neurons = range(sample_active_start,sample_active_end,sample_step)

sample_inactive_start = (i_cue_ang+180 % 360)
sample_inactive_end = sample_inactive_start+i_cue_width
sample_step = int((sample_inactive_end-sample_inactive_start)/sample_size)
sample_exc_neurons.append([n % 360 for n in range(sample_inactive_start, sample_inactive_end,sample_step)])

sample_inh_neurons = range(1,Ni,int(Ni/sample_size))



# specify the interneuron model
eqs_i = '''
dV/dt = (-gL_i*(V - El_i) - g_ext_i*s_ext*(V - E_ampa) \
- Gii*s_gaba*(V - E_gaba) \
- Gei*s_tot*(V - E_nmda)/(1 + b*exp(-a*V)))/Cm_i : volt
ds_ext/dt = -s_ext/tau_syn_ampa : 1
ds_gaba/dt = -s_gaba/tau_syn_gaba : 1
s_tot : 1
'''

# specify the pyramidal neuron model
eqs_e = '''
dV/dt = (-gL_p*(V - El_p) - g_ext_p*s_ext*(V - E_ampa) \
- Gie*s_gaba*(V - E_gaba) \
- Gee*s_tot*(V - E_nmda)/(1 + b*exp(-a*V)) + I_e)/Cm_p : volt
ds_ext/dt = -s_ext/tau_syn_ampa : 1
ds_gaba/dt = -s_gaba/tau_syn_gaba : 1
ds_nmda/dt = -s_nmda/tau_syn_nmda + alfa*x*(1 - s_nmda) : 1
dx/dt = -x/tau_syn_x : 1
s_tot : 1
I_e : 1 * amp
'''

# modify reset to increment auxilary variable x at the time of each spike
reset_nmda = '''
V = Vr_p
x+=1*1
'''

# return the distance between neurons on the ring
def circ_distance(deltaTheta):
    if (deltaTheta > 0):
        return min(deltaTheta,360-deltaTheta)
    else:
        return max(deltaTheta,deltaTheta-360)

currents = lambda i,j: i_cue_amp* \
exp(-0.5*circ_distance((i-j)*360./Np)**2/i_cue_width**2)

# precompute input current for the cue direction i_cue_ang
current_e=zeros(Np)
j = i_cue_ang*Np/360.
for i in xrange(Np):
    current_e[i]=currents(i,j)




# make the interneurons
Pi = NeuronGroup(N=Ni, model=eqs_i, threshold=Vt_i, reset=Vr_i, \
refractory=tauRef_i, order=2)
# deliver Poisson input to each neuron
PGi = PoissonGroup(Ni, fext) 
Cpi = IdentityConnection(PGi, Pi, 's_ext', weight=1.0)
# set up recurrent inhibitory connections
Cii = Connection(Pi, Pi, 's_gaba', weight=1.0)
# ???? why do we use the synaptic variable s_gaba to determine sparseness??
# 's_gaba' is not sparseness, it's the state variable



# make the pyramidal neurons
Pe = NeuronGroup(N=Np, model=eqs_e, threshold=Vt_p, reset=reset_nmda, \
refractory=tauRef_p, order=2)   # reset=reset_nmda causes overflow of spiking
# add the external Poisson drive to the E population
PGe = PoissonGroup(Np, fext)
Cpe = IdentityConnection(PGe, Pe, 's_ext', weight=1.0)
# set up recurrent inhibition from I to E cells with uniform strength = 1
Cie = Connection(Pi, Pe, 's_gaba', weight=1.0)


# Jm_ee is calculated from the normalization condition
tmp = sqrt(2*pi)*sigma_ee*erf(360.*0.5/sqrt(2.)/sigma_ee)/360.
Jm_ee = (1.-Jp_ee*tmp)/(1.-tmp)
weight = lambda i:(Jm_ee+(Jp_ee-Jm_ee)* \
exp(-0.5*(360.*min(i,Np-i)/Np)**2/sigma_ee**2))

weight_e = zeros(Np)
for i in xrange(Np):
    weight_e[i]=weight(i)

# weight_e*=Np/sum(weight_e)  # Normalization
fweight = rfft(weight_e)    # Fourier transform


# a vector representing the tuned angle of all neurons
orientation = lambda i : i*360./Np
orient_weight = numpy.zeros(Np)
for i in range(Np):
    orient_weight[i]=orientation(i)


#keep an exponentially moving average of number of spikes per neuron
decay = 0.99
spikebin = numpy.zeros(Np)


# keeps track of current angle encoded in the population
population_ang = 0



# implement the update of NMDA conductance s_tot on each time step
@network_operation(current_clock,when='start')
def update_nmda(current_clock):
#    print Pe.x
    fs_NMDA = rfft(Pe.s_nmda)   #s_NMDA = Pe.s_nmda.sum()
#    fs_tot = fs_NMDA*fweight    # fweight is a vector
#    s_tot = irfft(fs_tot)        
#    s_tot = s_tot.sum()
    Pe.s_tot = irfft(fweight * fs_NMDA)
#    print weight_e
    Pi.s_tot = Pe.s_nmda.sum()  # get same results if replace with Pe.s_tot.sum()

@network_operation(current_clock, when='start')
def update_current(current_clock):
    c_time = current_clock.t
    if stim_start < c_time < stim_stop:
        Pe.I_e = current_e
    else:
        Pe.I_e = 0

@network_operation(current_clock,when='end')
def update_population_stats(current_clock):
    global spikebin
    spikebin = (decay)*spikebin + (1-decay)*numpy.array(Pe.x)
    norm_spikebin = spikebin/(sum(spikebin))
    population_ang = sum(norm_spikebin * orient_weight)
    if current_clock.t > stim_stop:
        fpopv.write(str(population_ang) + '\n')

@network_operation(current_clock,when='end')
def update_neuron_stats(current_clock):
    if current_clock.t > stim_stop:
        for n in sample_exc_neurons:
            V = Pe.V[n] * mvolt
            s_tot = Pe.s_tot[n]
            s_ext = Pe.s_ext[n]
            
            #not quite sure if these are the right statistics
            current_exc = Gee*s_tot*(V - E_nmda)/(1 + b*exp(-a*V))  # recurrent excitation
            current_bck = g_ext_p*s_ext*(V - E_ampa) # background excitatory current
            current_inh = Gie*s_gaba*(V - E_gaba)   # recurrent inhibition

            fneuron.write(str(n)+',')
            fneuron.write(str(current_exc) + ',')
            fneuron.write(str(current_bck) + ',')
            fneuron.write(str(current_inh))
            fneuron.write(' ')

        fneuron.write('\n')
        for n in sample_inh_neurons:
            V = Pi.V[n]  * mvolt
            s_gaba = Pi.s_gaba[n]
            s_ext = Pi.s_ext[n]

            #not quite sure if these are the right statistics
            current_inh =  Gii*s_gaba*(V - E_gaba) 
            current_bck = g_ext_p*s_ext*(V - E_ampa)

            fneuron.write(str(n)+',')
            fneuron.write(str(current_inh) + ',')
            fneuron.write(str(current_bck))
            fneuron.write(' ') 
        fneuron.write('\n')


M = SpikeMonitor(Pe)
Q = SpikeMonitor(Pi)

Pi.V = Vr_i + rand(Ni) * (Vt_i - Vr_i)
Pe.V = Vr_p + rand(Np) * (Vt_p - Vr_p)


# open file for storing population vector statistics 
fpopv = open(str(Np) + POPULATION_OUTPUT_FILE_BASE + str(pa.seed) + '.dat','w+')     
fpopv.write(str(i_cue_ang)+'\n')            #write cue angle on first line

#open file for storing information on individual neuron activity
fneuron = open(str(Np) + NEURON_OUTPUT_FILE_BASE + str(pa.seed) + '.dat','w+')
fneuron.write(str(len(sample_exc_neurons)))
fneuron.write(str(len(sample_inh_neurons)))

run(sim_duration)

#close files
fpopv.close()                        


#subplot(211)
#raster_plot(M)
#subplot(212)
#raster_plot(Q)
#show()

