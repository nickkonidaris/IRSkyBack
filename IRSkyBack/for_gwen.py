

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

# Chuck's units: erg / s / cm2 / angstrom / spatial pixel for a 0.7" slit in ETC
ss = ss.to(ur.photon / ur.s / ur.cm**2 / ur.angstrom)

hc = cc.h * cc.c

Res = 4000 # with 2.5 pixel slit

# FRom Mirmos-1000-v3
Y = [851, 10500]
J = [685, 12530]
H = [523, 16370]
K = [396, 22000]

import numpy as np

MOS_PIX_SCALE = 0.18 # 0.18 arcsecond per pixel in spaital
MOS_SLT_WIDTH = 0.7 # Slit width
aerial = MOS_PIX_SCALE * MOS_SLT_WIDTH

# IRMOS spectrum<
beam_d = 118
frac = 0.9
for bandname in ['Y', 'J', 'H', 'K']:
    band = eval(bandname)
    lp, sp, spl, prof_l, prof = OBS.to_observed(ll, ss, band[1]/Res*AA,
                                                N=band[0]*beam_d * frac,
                                                nsamp=2.8)
    ep = (hc/lp).to(ur.erg) * spl * aerial
    np.savetxt("/Volumes/GoogleDrive/Team Drives/F-IRMOS/for_gwen/band_%s_v4.txt" % bandname, 
        np.array([lp, ep]).T)
    print(bandname, np.median(ep))



frac  = .3
lp, sp, spl, prof_l, prof = OBS.to_observed(ll, ss, 16000/3200*AA,
                                        N=125*110*frac*4,
                                        nsamp=2.8)

ep = (hc/lp).to(ur.erg) * spl * aerial
np.savetxt("/Volumes/GoogleDrive/Team Drives/F-IRMOS/for_gwen/mosfire_H_v4.txt", np.array([lp, ep]).T)

# 2: Added energy
# 3: Fixed typo with no scattering
# 4: Fixed aerial coverage






