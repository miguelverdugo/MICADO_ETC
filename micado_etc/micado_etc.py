import numpy as np
from spextra import Spextrum, Passband
import astropy.units as u

#set_target.emission_line()
#set_target.template()
#set_target.black_body()


class Set_Target:

    def __init__(self, redshift=0, magnitude=20, filter_curve="V", filter_system="Vega"):
        self.redshift = redshift
        self.filter_system = filter_system
        if self.filter_system == "Vega"
            mag_unit = u.mag
        else:
            mag_unit = u.ABmag

        self.magnitude = magnitude * mag_unit
        self.filter_curve = filter_curve

        self.spectrum = None
        self.spatial = None

    def template(self, template_name='pickles/a0v'):
        sp = Spextrum(template_name=template_name).redshift(self.redshift)
        self.spectrum = sp.scale_to_magnitude(amplitude=self.magnitude, filter_curve=self.filter_curve)

    def MARCS(self):
#        self.spectrum = MARCS spectrum

    def emission_line(self, center, fwhm, flux):

        sp = Spextrum.emission_line_spectra(center=center, fwhm=fwhm, flux=flux)
        self.spectrum = sp

    def black_body(self, temperature):
        sp = Spextrum.black_body_spectrum(temperature=temperature, amplitude=self.magnitude, filter_curve=self.filter_curve)

        self.spectrum = sp

    def powerlaw(self, alpha):
        sp = Spextrum.powerlaw(alpha=alpha)

    def sersic(self):
        self.spatial = "sersic"

    def set_template(self):
        source = "source"
        return source





