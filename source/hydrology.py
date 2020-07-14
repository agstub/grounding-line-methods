#-------------------------------------------------------------------------------
# This function defines:
# (1) volume change rate of the subglacial lake problem.
#    *Default = smoothed sawtooth wave
#
# (2) sea level change timeseries for marine ice sheet problem.
#    *Default = sinusoidal tidal cycle if 'tides' turned 'on', OR...
#             = zero if 'tides' turned 'off'
#
# trg(), sqr(), and swt() functions are from ybeltukov's answer on:
# https://mathematica.stackexchange.com/questions/38293/make-a-differentiable-smooth-sawtooth-waveform
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import scipy.misc as scm
from params import t_final,nt_per_year,tides,model_setup

#------------------Define rate of subglacial lake volume change-----------------
d0 = 0.1            # Smoothing parameter

# Smoothed triangle wave
def trg(t):
    return 1 - 2*np.arccos((1 - d0)*np.sin(2*np.pi*t))/np.pi

# Smooth square wave
def sqr(t):
    return 2*np.arctan(np.sin(2*np.pi*t)/d0)/np.pi

# Smoothed sawtooth wave
def swt(t):
    return (1 + trg((2*t - 1)/4)*sqr(t/2))/2

# Sawtooth volume change
def Vol(t,lake_vol_0):
    if model_setup == 'wedge_test':
        V = lake_vol_0*(1.0+1.0*t/t_final)**2
    else:
        V = 2.0*lake_vol_0*swt(t/(3.154e7))

    return V

# Compute rate of subglacial lake volume change.
def Vdot(lake_vol_0,t):
    # Corresponds to sawtooth volume change
    dt_fine = 3.154e7/5000.0
    Vd = scm.derivative(Vol,t,dx=dt_fine,args=(lake_vol_0,))
    return Vd


#---------------------Define sea level change timeseries------------------------
def sl_change(t):
    if tides == 'on':
        SLC = np.sin(4*np.pi*t/(3.154e7/12.0/30.0))
    else:
        SLC = 0.0
    return SLC
