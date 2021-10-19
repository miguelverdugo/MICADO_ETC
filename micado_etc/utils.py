from spextra import Spextrum
from scopesim_templates.basic.galaxy import galaxy
from scopesim_templates.basic.stars import star
from scopesim import Source
from astropy.table import Table


def template(template_name='pickles/a0v', magnitude=20, filter_curve="V", redshift=0):
    sp = Spextrum(template_name=template_name).redshift(redshift)
    return sp.scale_to_magnitude(amplitude=magnitude, filter_curve=filter_curve)


def from_file(filename, magnitude, filter_name, redshift):

    sp = Spextrum.from_file(filename=filename)
    sp = sp.redshift(z=redshift).scale_to_magnitude(amplitude=magnitude, filter_curve=filter_name)
    return sp


def MARCS():
    """
    TODO: Implement in spextra
    Parameters
    ----------
    self

    Returns
    -------

    """

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


def point_source(magnitude, sed, filter_name, redshift=0):
    sp = Spextrum(sed).redshift(z=redshift).scale_to_magnitude(amplitude=magnitude, filter_curve=filter_name)
    tbl = Table(names=["x", "y", "ref", "weight", "spec_types"],
                data=[[0], [0], [0], [1], [sed]])

    src = Source(spectra=[sp], table=tbl)
    # src.meta.update()

    return src


def sersic(magnitude, sed, filter_name,  redshift, r_eff, n, ellip=0.1, theta=0, extend=3):
    src = galaxy(sed=sed,
                 z=redshift,
                 amplitude=magnitude,
                 filter_curve=filter_name,
                 pixel_scale=0.0015,  # !!!!
                 r_eff=r_eff,
                 n=n,
                 ellip=ellip,
                 theta=theta,
                 extend=extend)
    return src


def uniform_flux():
    pass
