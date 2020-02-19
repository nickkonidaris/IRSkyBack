import numpy as np

""" Reference used :
    Woods, Thomas N., Raymond T. Wrigley, Gary J. Rottman,
        and Robert E. Haring.
    “Scattered-Light Properties of Diffraction Gratings.” Applied Optics
33, no. 19 (July 1, 1994): 4273. https://doi.org/10.1364/AO.33.004273. """


def fitted_profile(lam, lam_blaze, lam0, f, N, A_B=1e-9):
    """ lam is the wavelength, _blaze is blaze wavelength, lam0 is the on-axis
    wavelength, and f is a dimensionless number from 0 to 1 representing the
    ratio of the groove width to groove spacing, N is number of elements
    illuminated. Note typical surface relief gratings have f=0.25 and VPHG
    are f=0.9"""

    a = np.pi * lam/lam0
    b = np.pi*(lam - lam_blaze)/lam0 * f
    b0 = np.pi*(lam0 - lam_blaze)/lam0 * f

    Yfit = (np.sinc(b)/np.sinc(b0))**2 * 0.5/(N**2 * np.sin(a)**2) + A_B

    return Yfit


def fitted_profile_approx(lam, lam0, N, A_B=3e-9):
    """ lam, lam0, N, A_B=3e-9.
    N ~ 0.25 x N_actual
    Good when dlam is smaller than 10% lam0 """
    X = lam-lam0

    w = lam0/(1.414*np.pi*N)

    Yfit = w**2/(X**2 + w**2) + A_B

    return Yfit


def test_fitted_approx():
    # From Woods' paper testing PE grating #8 upper half
    # 2880 lines / mm, blazed at 250 nm, 58 m size, alpha = 58.8 deg
    # N = 25,000, f=0.36, A_B=3.1e-9

    lams = np.arange(-100, 100, 0.01) + 250

    GDF = fitted_profile_approx(lams, 250, 25000, A_B=3.1e-9)

    return (lams-250, GDF)
