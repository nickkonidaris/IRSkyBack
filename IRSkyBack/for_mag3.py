

from astropy import units as ur
from astropy import constants as cc


import importlib


import IRSkyBack.reader as reader
import observe as OBS

importlib.reload(OBS)


ll, ss = reader.load()

AA = ur.angstrom
Fl_units = ur.photon/ur.micron/ur.m**2/ur.second
ss = ss.to(Fl_units)
ss += 300*Fl_units

ss = ss.to(ur.photon / ur.s / ur.cm**2 / ur.angstrom)


# FRom Mirmos-1000-v3
Y = [851, 10500]
J = [685, 12530]
H = [523, 16370]
K = [396, 22000]

import numpy as np

# Mag3 Spectrum
Res = 4000 # with 3.0 pixel slit
beam_d = 118
frac = 0.9
for bandname in ['Y', 'J', 'H', 'K']:
    band = eval(bandname)
    lp, sp, spl, prof_l, prof = OBS.to_observed(ll, ss, band[1]/Res*AA,
                                                N=band[0]*beam_d * frac,
                                                nsamp=2.8)
    print(bandname, np.median(ep))



