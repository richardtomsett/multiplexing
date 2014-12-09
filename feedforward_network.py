from brian2 import *
from brianmodels import * # models stored as multi-line strings in brianmodels.py
from spiketrains import * # functions to generate spike trains
import time
import sys
import os
try:
  import cPickle as pickle
except:
  import pickle

#prefs.codegen.target = 'weave'         # Codegen doesn't seem to boost performance significantly

# Scale the inhibition as well as excitation?
if len(sys.argv) > 5:
  scale_inhibition = int(sys.argv[5])
else:
  scale_inhibition = False

# Set the random seed for the simulation
if len(sys.argv) > 4:
  rseed = int32(sys.argv[4])
else:
  rseed = 1234
np.random.seed(4321) # Initially the seed should always be the same, to create the same network every time

if len(sys.argv) > 3:
  connections_per_layer = int32(sys.argv[3])
else:
  connections_per_layer = 50

if len(sys.argv) > 2:
  model = sys.argv[2]
else:
  model = 'LIF'

if len(sys.argv) > 1:
  simtype = sys.argv[1]
else:
  simtype = 'rate'

if simtype == 'sync':
  pulse_connections = 50
  pulse_a = 2
  pulse_std = 5
  pulse_time = 300.
else:
  input_rate = 1500. * Hz

defaultclock.dt = 0.05 * ms             # timestep
num_layers = 9                          # number of layers
N_layer_e = 100                         # number of excitatory neurons per layer
N_layer_i = 25                          # number of inhibitory neurons per layer
I_rand = 0. * amp                       # set I_rand (in the model equations) to 0
FF_delay = 2. * ms                      # axonal delay time
weight_scale = np.linspace(.5, 2., 16)  

# The excitatory neuron group
print('Creating neurons... ')
N_e = N_layer_e * num_layers
N_i = N_layer_i * num_layers

if model.startswith('ML'):
  eqs   = eqs_ML  + eqs_g_exp
  P_E = NeuronGroup(N_e, eqs, threshold=eqs_threshold, refractory=eqs_refractory, method="rk4")
  area = 1000. * umetre**2
  P_E.C    = ( 2. * uF * cm**-2) * area
  P_E.g_Na = (20. * msiemens * cm**-2) * area
  P_E.g_K  = (20. * msiemens * cm**-2) * area
  P_E.E_Na =  50. * mV
  P_E.E_K  = -100. * mV
  P_E.E_l  = -70. * mV
  P_E.phi = 0.15
  P_E.V_1 = -1.2 * mV
  P_E.V_2 = 23. * mV
  P_E.V_4 = 21. * mV
  P_E.alpha = 0.005 / ms
  P_E.gamma = 5. * mV
  P_E.V_3 = 10. * mV
  P_E.g_l = 2. * msiemens * cm**-2 * area

  if model.endswith('M'):
    P_E.beta = -35. * mV
    P_E.g_adapt = 2. * msiemens * cm**-2 * area
    FF_weight_E = .575 * nS
  elif model.endswith('AHP'):
    P_E.beta = 0. * mV
    P_E.g_adapt = 15. * msiemens * cm**-2 * area
    FF_weight_E = .52 * nS

elif model == 'LIF':
  eqs   = eqs_LIF + eqs_g_exp
  P_E = NeuronGroup(N_e, eqs, threshold='v >= v_t', reset='v = v_r', refractory=2*ms, method="rk4")
  P_E.E_l = -70. * mV
  P_E.g_l =  20 * nS
  P_E.C   = 200. * pF
  P_E.v_t = -57. * mV
  P_E.v_r = -70. * mV
  FF_weight_E = .76 * nS
  
P_E.I_inject = 0. * nA
P_E.E_syn_e = 0. * mV
P_E.E_syn_i = -70. * mV
P_E.tau_syn_e = 1.5 * ms
P_E.tau_syn_i = 10. * ms

eqs_I = eqs_LIF + eqs_g_exp
P_I = NeuronGroup(N_i, eqs_I, threshold='v >= v_t', reset='v = v_r', refractory=2*ms, method="rk4")
P_I.E_l = -70. * mV
P_I.g_l =  20 * nS
P_I.C   = 200. * pF
P_I.v_t = -57. * mV
P_I.v_r = -70. * mV
P_I.I_inject = 0. * nA
P_I.E_syn_e = 0. * mV
P_I.E_syn_i = -70. * mV
P_I.tau_syn_e = 1.5 * ms
P_I.tau_syn_i = 10. * ms
FF_weight_I = .76 * nS

# Initialise membrane potentials
P_E.v = P_E.E_l
P_I.v = P_I.E_l

# Create spike monitors to record spikes from the neuron groups
s_mon_E = SpikeMonitor(P_E)
s_mon_I = SpikeMonitor(P_I)
vmon = StateMonitor(P_E, ['v', 'I_syn'], 0)

# Create feedforward connections
print('Creating feedforward connections... ')
C_FF_EE = Synapses(P_E,P_E,model='w : siemens', pre='g_e += w') # excitatory to excitatory
C_FF_EI = Synapses(P_E,P_I,model='w : siemens', pre='g_e += w') # excitatory to inhibitory
C_FF_IE = Synapses(P_I,P_E,model='w : siemens', pre='g_i += w') # inhibitory to excitatory
for i in range(num_layers - 1):
  potential_sources = np.arange(N_layer_e)+(i*N_layer_e)
  for j_e in range(N_layer_e):
    target = j_e + N_layer_e * (i+1)
    chosen_sources = np.random.choice(potential_sources, connections_per_layer, replace=False)
    C_FF_EE.connect(chosen_sources, target)
  for j_i in range(N_layer_i):
    target = j_i + N_layer_i * (i+1)
    chosen_sources = np.random.choice(potential_sources, connections_per_layer, replace=False)
    C_FF_EI.connect(chosen_sources, target)
    C_FF_IE.connect(target - N_layer_i, potential_sources)
C_FF_EE.delay = FF_delay
C_FF_EI.delay = FF_delay
C_FF_IE.delay = FF_delay
# Set initial weights to zero (modified prior to each simulation run)
C_FF_EE.w = 0. * nS
C_FF_EI.w = 0. * nS
C_FF_IE.w = 0. * nS

# Create the external random Poisson input
print('Creating external input... ')
external_rate_e = 3000 * Hz
P_external = PoissonGroup(N_e+N_i, external_rate_e)
C_extE = Synapses(P_external, P_E, model='w : siemens', pre='g_e += w')
C_extI = Synapses(P_external, P_I, model='w : siemens', pre='g_e += w')
C_extE.connect(np.arange(N_e), np.arange(N_e))
C_extI.connect(np.arange(N_i)+N_e, np.arange(N_i))
C_extE.w = FF_weight_E
C_extI.w = FF_weight_I

# Create the signal input
if simtype == 'sync':
  print('Creating pulse packet inputs... ')
  potential_sources = np.arange(N_layer_e)
  chosen_sources_e = []
  chosen_sources_i = []
  for j_e in range(N_layer_e):
    chosen_sources_e.append(np.random.choice(potential_sources, pulse_connections, replace=False))
  for j_i in range(N_layer_i):
    chosen_sources_i.append(np.random.choice(potential_sources, pulse_connections, replace=False))

  start_time = time.time()
  (ids, times) = pulse_packet(N_layer_e, pulse_a, pulse_time, pulse_std)
  P_pulse = SpikeGeneratorGroup(N=N_layer_e, indices=ids, times=times * ms)
  C_pulseE = Synapses(P_pulse, P_E, model='w : siemens', pre='g_e += w')
  C_pulseI = Synapses(P_pulse, P_I, model='w : siemens', pre='g_e += w')
  for j_e in range(N_layer_e):
    C_pulseE.connect(chosen_sources_e[j_e], j_e)
  for j_i in range(N_layer_i):
    C_pulseI.connect(chosen_sources_i[j_i], j_i)
  # Set initial weights to zero (modified prior to each simulation run)
  C_pulseE.w = 0. * nS
  C_pulseI.w = 0. * nS
else:
  print('Creating rate inputs... ')
  P_rate = PoissonGroup(N_layer_e+N_layer_i, input_rate)
  C_rateE = Synapses(P_rate, P_E, model='w : siemens', pre='g_e += w')
  C_rateI = Synapses(P_rate, P_I, model='w : siemens', pre='g_e += w')
  C_rateE.connect(np.arange(N_layer_e), np.arange(N_layer_e))
  C_rateI.connect(np.arange(N_layer_i)+N_layer_e, np.arange(N_layer_i))
  # Set initial weights to zero (modified prior to each simulation run)
  C_rateE.w = 0. * nS
  C_rateI.w = 0. * nS 

# Simulation control: only turn on excitatory synapses after 0.1 seconds to reduce network startup
# effects. Switch on rate input at 0.6 seconds if this is a rate increase simulation.
@network_operation()
def syn_on():
  if defaultclock.t == 0.1 * second:
    C_FF_EE.w = scaled_weight_EE
  if defaultclock.t == 0.6 * second and simtype != 'sync':
    C_rateE.w = scaled_weight_EE
    C_rateI.w = scaled_weight_EI

# Store the initialised network so we can reset it for repeated simulation runs
net = Network(collect())
net.store()

print('Beginning simulations... ')
# Create arrays to store the recorded spike times
spike_times_E = []
spike_times_I = []

# Run a simulation per synaptic weight scaling
for iW in xrange(len(weight_scale)):
  start_time = time.time() # time the simulation runs to see if it's worth going to get a coffee
  
  # Calculate the scaled synaptic weights
  scaled_weight_EE = FF_weight_E * weight_scale[iW]
  if scale_inhibition == 1:
    scaled_weight_EI = FF_weight_I * weight_scale[iW] * 3.5
    scaled_weight_IE = FF_weight_E * weight_scale[iW]
  elif scale_inhibition == -1:
    scaled_weight_EI = 0. * nS
    scaled_weight_IE = 0. * nS
  else:
    scaled_weight_EI = FF_weight_I * 3.5
    scaled_weight_IE = FF_weight_E * weight_scale[iW]
  
  # Restore the network and set the synaptic weights
  net.restore()
  np.random.seed(rseed) # set the seed to the user-specified seed for random background inputs
  C_FF_EE.w = scaled_weight_EE
  C_FF_EI.w = scaled_weight_EI
  C_FF_IE.w = scaled_weight_IE * 2.
  if simtype == 'sync':
    C_pulseE.w = scaled_weight_EE
    C_pulseI.w = scaled_weight_EI
    net.run(.5 * second)
  else:
    net.run(1.1 * second)

  spike_times_E.append([  s_mon_E.i.get_item(slice(None), level=1), \
                          s_mon_E.t_.get_item(slice(None), level=1)])
  spike_times_I.append([  s_mon_I.i.get_item(slice(None), level=1), \
                          s_mon_I.t_.get_item(slice(None), level=1)])
  print('Done simrun ' + str(iW+1) + ', time: ' + str(time.time()-start_time))

# Save the spike times in a pickled file
print('Saving spike times... ')
if not os.path.exists(model):
    os.makedirs(model)
fname = model + '/synfire_noinhib_sync_' + str(rseed) + '_' + str(connections_per_layer) + '.pic'
with open(fname,'wb') as fp:
  pickle.dump(spike_times_E, fp)
  pickle.dump(spike_times_I, fp)
