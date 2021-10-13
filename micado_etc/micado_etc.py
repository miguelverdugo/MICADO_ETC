import numpy as np

import astropy.units as u
from astropy.table import Table

from spextra import Spextrum, Passband

from .utils import *

"""
Expected usage:

initialize
target = SetTarget(magnitude=15, filter_curve="R", redshift=1)

target.template(template_name="kc96/s0")  # defines the spectra, default a vega spectrum
terget.sersic(parameters) # defines a spatial distribution  # default point source 

target.get_source # get the source

"""


SED_dict = dict(template=template,
                MARCS=MARCS,
                emission_line=emission_line,
                black_body=black_body,
                powerlaw=powerlaw,
                uniform=flat_spec)

FLUX_DISTRO_dict = dict(sersic=sersic,
                        point_source=point_source,
                        uniform=uniform_flux)



class ETC:

    def __init__(self, mode="imaging", pixel_size=0.004, ao_mode="SCAO"):

        self.mode = mode
        self.pixel_size = pixel_size
        self.ao_mode = ao_mode

        self.spectrum = None
        self.source = None
        self.sky = None

    def set_instrument(self):
        """
        Reset instrument and perform checks

        """
        pass

    def set_sed(self, spectrum_type, **kwargs):
        if spectrum_type not in SED_dict.keys():
            raise ValueError("SED not available")

        self.spectrum = SED_dict[spectrum_type](**kwargs)

    def set_target_flux_distribution(self, distribution, **kwargs):
        if distribution not in FLUX_DISTRO_dict.keys():
            raise ValueError("FLUX distribution not availalbe")

        self.source = FLUX_DISTRO_dict[distribution](self.spectrum, **kwargs)

    def set_sky_conditions(self):
        """
        TODO: Check how ScopeSim can read this
        :return:
        """

        self.sky = None


