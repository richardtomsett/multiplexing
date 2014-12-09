# Contains functions for generating different kinds of input spike train
# for the Brian2 SpikeGeneratorGroup class. Currently only pulse_packet,
# but will add others as needed

import numpy as np

# Generator for pulse packet, adapted from brian/examples/misc/pulsepacket.py
def pulse_packet(n, a, t, sigma):
  total_spikes = n * a
  times = np.zeros(total_spikes)
  ids   = np.zeros(total_spikes)
  for i in xrange(n):
    times[(i*a):(i*a)+a] = np.random.normal(t, sigma, a)
    ids[(i*a):(i*a)+a]   = np.ones(a) * i
  return (ids, times)