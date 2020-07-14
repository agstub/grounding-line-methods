# This script computes the grounding line timeseries in the limit of
# very slow flow (assuming a flat ice-water interface), given the bed geometry
# and volume change timeseries.

import sys
sys.path.insert(0, './source')
from geometry import bed,interface
from params import X_fine, Lngth,t_final,nt_per_year,tol
import scipy.integrate as scpint
from hydrology import Vol
import scipy.optimize as optimize
import numpy as np


lake_vol_0 = scpint.quad(lambda x: interface(x)-bed(x),0,Lngth,full_output=1)[0]

def get_glines(z):
    # Computes minimum and maximum grounding line positions given the
    # lower surface elevation function s and the bed geometry.

    s = np.maximum(z*np.ones(np.size(X_fine)),bed(X_fine))
    s_new = np.copy(s)

    x = np.copy(X_fine)
    x_new = np.copy(X_fine)

    key = 1e10

    # Mark points on ice-bed boundary (by assigning a ridiculously large value)
    s[s-bed(x)<tol] = key
    s_new[s_new-bed(x)<tol] = key

    # Loop over points on lower boundary
    for j in range(np.shape(s)[0]-1):
        if s[j+1] < 0.9*key and  s[j-1] < 0.9*key:
        # If *both* neighboring points are on the ice-water boundary, then
        # mark these points! (these cannot be grounding lines)
            s_new[j] = key

    # Mark last point
    if s[-2] < 0.9*key:
        s_new[-1] = key

    # All points except grounding lines have now been marked.
    glines = x_new[s_new<0.9*key] - 0.5*Lngth
    glines_pos = glines[glines>0]
    glines_neg = glines[glines<0]
    if np.size(glines)>0:
        XR = np.sort(glines_pos)[0]       # Right grounding line
        XL = -np.sort(-glines_neg)[0]     # Left grounding line
    else:
        XL = 0              # Left grounding line
        XR = Lngth          # Right grounding line

    return XL/1000.0,XR/1000.0


def R(z,t):

    func = Vol(t,lake_vol_0)-scpint.quad(lambda x: np.maximum(z,bed(x))-bed(x),0,Lngth,full_output=1)[0]

    return func


t_arr = np.linspace(0,t_final,num=int(nt_per_year*t_final/3.154e7))

Xg = np.zeros(np.size(t_arr))
res = np.zeros(np.size(t_arr))

for i in range(np.size(t_arr)):
    #print(i)
    if i == 0:
        x0 = -5
    else:
        x0 = res[i-1]

    res[i] = optimize.newton(R, x0, args=(t_arr[i], ),maxiter=100,disp=False)

    XL,XR = get_glines(res[i])

    Xg[i] = XR

np.savetxt('Xg_slow',Xg)
