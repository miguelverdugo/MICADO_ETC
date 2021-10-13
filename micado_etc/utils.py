from spextra import Spextrum
from scopesim_templates.basic.galaxy import galaxy
from scopesim import Source





def template(template_name='pickles/a0v', magnitude=20, filter_curve="V", redshift=0):
    sp = Spextrum(template_name=template_name).redshift(redshift)
    return sp.scale_to_magnitude(amplitude=magnitude, filter_curve=filter_curve)


def MARCS(self):
    #  self.spectrum = MARCS spectrum
    return NotImplementedError


def emission_line(center, fwhm, flux):
    sp = Spextrum.emission_line_spectra(center=center, fwhm=fwhm, flux=flux)
    return sp


def black_body(temperature, magnitude, filter_curve):
    sp = Spextrum.black_body_spectrum(temperature=temperature, amplitude=magnitude, filter_curve=filter_curve)
    return sp


def powerlaw(alpha, magnitude, filter_curve):
    sp = Spextrum.powerlaw(alpha=alpha, amplitude=magnitude, filter_curve=filter_curve)
    return sp


def flat_spec(magnitude):
        sp = Spextrum.flat_spectrum(amplitude=magnitude)
        return sp


def point_source():
    pass


def sersic(sed, redshift, magnitude, filter_curve, r_eff, n, ellip=0.1, theta=0, extend=3):
    src = galaxy(sed=sed,
                 z=redshift,
                 amplitude=magnitude,
                 filter_curve=filter_curve,
                 pixel_scale=0.0015,  # !!!!
                 r_eff=r_eff,
                 n=n,
                 ellip=ellip,
                 theta=theta,
                 extend=extend)
    return src


def uniform_flux():
    pass
