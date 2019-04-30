import numpy as np
import os

import scipy
from astropy import units as ur


def load():
    " Returns Table 4 from Olivia+ 2015 "
    dir_path = os.path.dirname(os.path.realpath(__file__))

    dat = np.loadtxt("%s/../etc/Oliva2015/table4.dat" % dir_path)
    ll, ss = dat.T

    ll = ll * ur.angstrom
    ss = ss * ur.photon / ur.meter**2 / ur.second / ur.Angstrom

    if False:
        MD = np.loadtxt("%s/../etc/Mosfire-sky-background/mosfire-130429.dat"
                        % dir_path)
        llm, ssm = MD.T
        llm = (llm*1.62695+14500-1.62695)

    return ll, ss
