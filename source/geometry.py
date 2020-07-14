#-------------------------------------------------------------------------------
# Define bed topography and initial ice-water interface functions.
# Note: Bed and ice-water interface should be equal on margins of the domain
# for the lake problem! They should be equal on one margin of the domain for the
# marine ice sheet problem (i.e., the grounded portion).
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from params import Lngth,Hght,model,model_setup

#-------------------- Generate Bed Topography-----------------------------------
def bed(x):

    #-------------------Default subglacial lake bed geometry--------------------
    if model == 'lake':
        # Generate Gaussian bed topography

        if model_setup == 'wedge_test':
            Bed = 10*np.abs(x-0.5*Lngth)/(0.5*Lngth)-10

        else:
            Bed = -10.0*(np.exp((-(x-Lngth/2.0)**2)/((0.25*Lngth)**2)))

    #-------------------Default marine ice sheet bed geometry-------------------
    elif model == 'marine':
        # Generate linear bed topography
        a0 = 5.0
        Bed = 2.5 - a0*x/Lngth
    return Bed

#-------------------------------------------------------------------------------

#------------------Generate initial ice-water/ice-bed interface-----------------
def interface(x):
    if model == 'lake':
        Int = 0.5*(bed(x) - 5 + np.abs(bed(x) - (-5)))

        if model_setup == 'wedge_test':
            Int = 0.5*(bed(x) - 7.5 + np.abs(bed(x) - (-7.5)))


    elif model == 'marine':
        Int = 0.5*(bed(x)  + np.abs(bed(x) ))
    return Int
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#
import matplotlib.pyplot as plt
from params import X_fine

# plt.plot((X_fine-0.5*Lngth)/1000.0,bed(X_fine))
# plt.plot((X_fine-0.5*Lngth)/1000.0,interface(X_fine))
# plt.show()
