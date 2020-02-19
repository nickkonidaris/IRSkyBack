
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


# From Mirmos-1000-v3
Y = [851, 10500]
J = [685, 12530]
H = [523, 16370]
K = [396, 22000]


# What Gwen Wants:
# 1 pixel slit, 2 pixel, and 2.5 pixel slit.

import numpy as np

MOS_PIX_SCALE = 0.18 # 0.18 arcsecond per pixel in spaital
MOS_SLT_WIDTH = 0.7 # Slit width
aerial = MOS_PIX_SCALE * MOS_SLT_WIDTH

# IRMOS spectrum<
Res = 4000 # with 2.5 pixel slit
slitwidth = 2.5
beam_d = 118
frac = 0.9

import datetime

now = datetime.datetime.now()


for bandname in ['Y', 'J', 'H', 'K']:
    band = eval(bandname)
    for Res in np.arange(2000, 5500, 500):
        for nsamp in [1.0, 2.0, 2.5]:
            lp, sp, spl, prof_l, prof = OBS.to_observed(ll, ss, band[1]/Res*AA,
                                                        N=band[0]*beam_d * frac,
                                                        nsamp=nsamp)
            epp = (hc/lp).to(ur.erg)
            ep = epp * spl * aerial
            header = f"""
Created on {now}

Spectral Resolution: {Res} for a slit of {slitwidth} pixel.

Slit is {nsamp} wide.
Assume {frac}x lines are illuminated.
Beam diameter {beam_d} mm


Delivered resolution is thus {Res*slitwidth/nsamp}
Units per column:
    {lp.unit} 
    {(epp.unit/ur.photon*ss).unit}
            """
            np.savetxt("/Volumes/GoogleDrive/Team Drives/F-IRMOS/for_gwen/band_%s_%s_%s_v5.txt" % \
                    (bandname, Res, nsamp), np.array([lp, ep]).T,
                    header=header)
            print(bandname, np.median(ep), Res, nsamp)



frac  = .3
lp, sp, spl, prof_l, prof = OBS.to_observed(ll, ss, 16000/3200*AA,
                                        N=125*110*frac*4,
                                        nsamp=2.8)

ep = (hc/lp).to(ur.erg) * spl * aerial
np.savetxt("/Volumes/GoogleDrive/Team Drives/F-IRMOS/for_gwen/mosfire_H_v4.txt", np.array([lp, ep]).T)

# 2: Added energy
# 3: Fixed typo with no scattering
# 4: Fixed aerial coverage






