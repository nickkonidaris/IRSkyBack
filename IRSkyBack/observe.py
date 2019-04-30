import numpy as np

from astropy import units as ur
import pylab as pl

import astropy.convolution as CC
from scipy.signal import convolve as scipy_convolve


import IRSkyBack.lsf_grating as lsf_grating

import importlib


import IRSkyBack.reader as reader


ll, ss = reader.load()
# ss = ss.to(ur.photon / ur.micron / ur.meter**2 / ur.second)


AA = ur.angstrom
Hband = (14600*AA, 18070*AA)

ROI = (ll > Hband[0]) & (ll < Hband[1])

ll = ll[ROI]
ss = ss[ROI]
dl = np.median(np.diff(ll))

Fl_units = ur.photon/ur.AA/ur.m**2/ur.second


def to_new_fwhm_box1d(ll, ss, width, nsamp=3):
    """ Use box kernel to convolve down and resample spectrum """

    nres = np.round(width/dl)
    npix = np.int(np.round(nres/nsamp))
    # print("Width of {0:1.4f} is {1:1.2f} pix".format(width, nres))
    kernel = CC.Box1DKernel(nres)
    sp = scipy_convolve(ss, kernel, mode='same', method='direct')

    return ll, sp, npix


def create_GRF(l0, N=125*1200*0.15, A_B=3e-9):
    """ Create a GRF with a line density of N (recall 15% is lower limit from Woods)
    and A_B Lorentizan floor """

    dl = np.arange(-50, +50, 0.1)
    l0d = l0.to(AA).value

    GRF = lsf_grating.fitted_profile_approx(dl+l0d, l0d, N, A_B)
    GRF /= np.sum(GRF)
    kernel = CC.Kernel1D(array=GRF, width=1)

    return dl, GRF, kernel


def to_observed(ll, ss, width, N=125*1200, A_B=3e-9, nsamp=3):

    l0 = np.median(ll).to(AA)
    dl, GRF, kernel = create_GRF(l0, N, A_B)

    _, sp, npix = to_new_fwhm_box1d(ll, ss, width, nsamp=nsamp)

    spl = scipy_convolve(sp, kernel, mode='same', method='direct')

    return ll[::npix], sp[::npix], spl[::npix], dl, GRF


def observe_maihara_region():

    pl.clf()
    pl.figure(2, figsize=(10,5))

    pl.clf()
    pl.figure(1, figsize=(10,5))

    for frac in [0.3]:
        # MOSFIRE: lpm=110, order=3
        # lpm=54 is FIRE, order=16, 3.3 pixel slit

        lp, sp, spl, prof_l, prof = to_observed(ll, ss, 16000/3200*AA,
                                                N=125*110*frac*4,
                                                nsamp=2.8)
        pl.figure(1)
        pl.step(lp, -2.5*np.log10(spl), color='red')
        pl.figure(2)
        pl.step(lp, spl, color='red')

        lp, sp, spl, prof_l, prof = to_observed(ll, ss, 16000/3200*AA,
                                                N=125*512*frac*1,
                                                nsamp=2.8)
        pl.figure(1)
        pl.step(lp, -2.5*np.log10(spl), color='orange')
        pl.figure(2)
        pl.step(lp, spl, color='orange')

        lp, sp, spl, prof_l, prof = to_observed(ll, ss, 16000/9000*AA,
                                                N=125*1440*frac*1,
                                                nsamp=2.8)
        # floor = 10**(-3.75/2.5)
        # pl.step(lp, -2.5*np.log10(spl+floor), color='orange')

        lp, sp, spl, prof_l, prof = to_observed(ll, ss, 16000/6000*AA,
                                                N=50*54*frac*16,
                                                nsamp=3.3)
        pl.figure(1)
        pl.step(lp, -2.5*np.log10(spl), color='blue')
        pl.figure(2)
        pl.step(lp, spl, color='blue')

    pl.figure(1)
    pl.step(lp, -2.5*np.log10(sp), lw=2, color='grey', alpha=0.5)
    pl.ylim(5.5, -1.5)
    pl.xlim(16400, 16900)
    pl.axvline(16600, color='red')
    pl.axvline(16700, color='red')
    pl.axhline(3+.75, color='black')
    pl.grid(True)

    pl.legend(["MOSFIRE", "IRMOS", "FIRE", "FIRE no scatter"])
    pl.xlabel(r"$Wavelength [\AA]")
    pl.ylabel(r"Magnitude")

    pl.figure(2)
    pl.step(lp, sp, lw=2, color='grey', alpha=0.5)
    pl.xlim(16400, 16900)
    pl.ylim(-0.1, 0.7)
    pl.axvline(16600, color='red')
    pl.axvline(16700, color='red')
    pl.axhline(3+.75, color='black')
    pl.grid(True)

    pl.legend(["MOSFIRE", "IRMOS", "FIRE", "FIRE no scatter"])
    pl.xlabel(r"$Wavelength [\AA]")
    pl.ylabel(r"Signal (linear)")



def fraction_of_pixels_above(frac=0.3, thresh=3, THPT=0.4, scale=0.3,
                             nsamp=2.8, Lam0=16000, et=5*ur.minute):
    """ frac[0-1] is grating quality, thresh is threshold above background level
    Returns resolution, fraction of pixels above thresh(res), and signal
    level(res)"""

    Rs = 10**np.linspace(np.log10(100), np.log10(15000), 35)
    Fracs = np.zeros_like(Rs)
    Sigs = np.zeros_like(Rs)

    bgd = 0.03 * Fl_units
    atel = 30 * ur.m**2
    skyarea = scale**2

    for ix, R in enumerate(Rs):
        lp, sp, spl, prof_l, prof = to_observed(ll, ss, 16000/R*AA,
                                                N=125*512*frac*3, nsamp=nsamp)
        spl *= Fl_units

        BW = Lam0 * AA / R / nsamp
        sig = (THPT * bgd * BW * atel * skyarea).to(ur.photon/ur.second)

        Sigs[ix] = sig.value

        G = BW * atel * skyarea * et
        npp_sky = spl * G
        npp_bgd = bgd * G
        ratio = npp_sky/npp_bgd

        Fracs[ix] = len(np.where(ratio > thresh)[0])/len(spl)

    return Rs, Fracs, Sigs*ur.photon/ur.second


def create_figure_for_irmos_v1():
    """ """

    scales = [1.4/2.8, 1.0/2.8, 0.9/2.8]
    pl.figure(3, figsize=(8, 4))
    pl.clf()
    pl.figure(2, figsize=(8, 4))
    pl.clf()
    pl.figure(1, figsize=(8, 4))
    pl.clf()

    THPT = .35
    GOAL = 2.0
    THRESH = 3
    GRATING = 0.9

    pl.figure(2)
    pl.xlabel("Spectral Resolution")
    pl.ylabel(r"$\frac{\left(%s x RN\right)^2}{sky\ flux} [min]$" % GOAL)
    pl.title("Time to sky shot noise %sx RN for Olivia bgd $300\ \gamma/m^2/s/as^2/\mu m$" % GOAL)
    labels = []
    for RN in [3.0, 5.5]:
        for scale in scales:
            labels.append("RN: %s, Scale: %1.2f\"/px" % (RN, scale))
            R, f, s = fraction_of_pixels_above(thresh=THRESH, frac=GRATING,
                                               scale=scale, et=0.1*ur.min,
                                               THPT=THPT)

            pl.figure(2)

            # BGD Limit is sqrt(signal x exptime) == RN*goal
            # exptime == (RN*goal)**2/signal
            goal = (RN*GOAL)
            exptime = goal**2/s.value*ur.second

            pl.plot(R, exptime.to(ur.minute))
            pl.grid()

    pl.legend(labels)
    pl.grid(True)
    pl.ylim(0,50)

    pl.figure(1)
    threshs = [3, 10]
    gratings  = {"vphg": .9, "surface relief": 0.3}
    legends = []
    for thresh in threshs:
        for gname, frac in gratings.items():
            legends.append("Grating: %s thresh: %sx" % (gname, thresh))
            R, f, s = fraction_of_pixels_above(thresh=thresh, frac=frac,
                                           scale=scale, et=1*ur.min,
                                               THPT=THPT)


            pl.plot(R, 1-f, '+-')
            pl.ylim(0, 1)
            pl.grid()

    pl.xlabel("Spectral resolution")
    pl.ylabel("Fraction of pixels some threshold above Olivia bgd")
    pl.title("Fraction of pixels above the background level for different R")
    pl.legend(legends)
    pl.title("Fraction of pixels above background")
    pl.grid(True)


    pl.figure(3)
    pl.title("Noise")
    thresh = 2
    scale = 0.43
    et = 120*ur.second
    R, f, s = fraction_of_pixels_above(thresh=thresh, frac=0.3,
                                        scale=scale, et=1*ur.min,
                                        THPT=THPT)
    s = s/ur.second
    shotnoise = np.sqrt(s * et / thresh)

    pl.loglog(R, shotnoise)
    pl.axhline(3)

    fmt = "Noise bgd thresh: %sx | thp: %s | et: %s | scale: %sas/pix | grating: %s"
    legends = [fmt % (thresh, THPT, et, scale, "surface relief"),
               "Read noise"]
    pl.legend(legends)
    pl.ylabel("Noise")
    pl.grid(True)



def observe_maihara_region_custom(combos=[ (3200, .3, "k")], labels=[""]):

    pl.figure(1)

    for combo in combos:
        R, frac, fmt = combo
        lp, sp, spl, prof_l, prof = to_observed(ll, ss, 16000/R*AA,
                                                N=125*512*frac*3, nsamp=2.8)
        pl.step(lp, -2.5*np.log10(np.abs(spl)), fmt)
        #pl.step(lp, spl)

    pl.ylim(5, -2)
    #pl.ylim(-0.1, 0.5)
    pl.legend(labels)
    pl.xlim(16400, 16900)
    pl.axvline(16600, color='red')
    pl.axvline(16700, color='red')
    pl.axhline(3+.75, color='black')
    pl.grid(True)

    pl.xlabel(r"Wavelength [$\rm \AA$]")
    pl.ylabel(r"$F_\lambda$ [magnitude]")


def create_figure_for_irmos_vSAC():
    """ """

    scales = [1.4/2.8, 1.0/2.8, 0.9/2.8]
    pl.figure(3, figsize=(8, 4))
    pl.clf()
    pl.figure(2, figsize=(8, 4))
    pl.clf()
    pl.figure(1, figsize=(8, 4))
    pl.clf()

    THPT = .35
    GOAL = 2.0
    THRESH = 3
    GRATING = 0.9

    pl.figure(2)
    pl.xlabel("Spectral Resolution")
    pl.ylabel(r"$\frac{\left(%s x RN\right)^2}{sky\ flux} [min]$" % GOAL)
    pl.title("Time to sky shot noise %sx RN for Olivia bgd $300\ \gamma/m^2/s/as^2/\mu m$" % GOAL)
    labels = []
    for RN in [3.0, 5.5]:
        for scale in scales:
            labels.append("RN: %s, Scale: %1.2f\"/px" % (RN, scale))
            R, f, s = fraction_of_pixels_above(thresh=THRESH, frac=GRATING,
                                               scale=scale, et=0.1*ur.min,
                                               THPT=THPT)

            pl.figure(2)

            # BGD Limit is sqrt(signal x exptime) == RN*goal
            # exptime == (RN*goal)**2/signal
            goal = (RN*GOAL)
            exptime = goal**2/s.value*ur.second

            pl.plot(R, exptime.to(ur.minute))
            pl.grid()

    pl.legend(labels)
    pl.grid(True)
    pl.ylim(0,50)

    pl.figure(1)
    threshs = [3]
    gratings  = {"vphg": .9, "surface relief": 0.3}
    legends = []
    for thresh in threshs:
        for gname, frac in gratings.items():
            legends.append("Grating: %s thresh: %sx" % (gname, thresh))
            R, f, s = fraction_of_pixels_above(thresh=thresh, frac=frac,
                                           scale=scale, et=1*ur.min,
                                               THPT=THPT)


            pl.plot(R, 1-f, '+-')
            pl.ylim(0, 1)
            pl.grid()

    pl.xlabel("Spectral resolution")
    pl.ylabel("Fraction of pixels some threshold above Olivia bgd")
    pl.title("Fraction of pixels above the background level for different R")
    pl.legend(legends)
    pl.title("Fraction of pixels above background")
    pl.grid(True)


    pl.figure(3)
    pl.title("Noise")
    thresh = 2
    scale = 0.43
    et = 120*ur.second
    R, f, s = fraction_of_pixels_above(thresh=thresh, frac=0.3,
                                        scale=scale, et=1*ur.min,
                                        THPT=THPT)
    s = s/ur.second
    shotnoise = np.sqrt(s * et / thresh)

    pl.loglog(R, shotnoise)
    pl.axhline(3)

    fmt = "Noise bgd thresh: %sx | thp: %s | et: %s | scale: %sas/pix | grating: %s"
    legends = [fmt % (thresh, THPT, et, scale, "surface relief"),
               "Read noise"]
    pl.legend(legends)
    pl.ylabel("Noise")
    pl.grid(True)
