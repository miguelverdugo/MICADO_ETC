import numpy as np

import astropy.units as u


from spextra import Spextrum, Passband

import scopesim as sim

from .utils import *


"""
Expected usage:

initialize
target = SetTarget(magnitude=15, filter_curve="R", redshift=1)

target.template(template_name="kc96/s0")  # defines the spectra, default a vega spectrum
terget.sersic(parameters) # defines a spatial distribution  # default point source 

target.get_source # get the source

"""
sim.server.database.download_package(["locations/Armazones.zip",
                                      "telescopes/ELT.zip",
                                      "instruments/MAORY.zip",
                                      "instruments/MICADO.zip"])

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

    def __init__(self, mode="imaging", pixel_size=0.004, ao_mode="SCAO", target_separation=0):

        self.mode = mode
        self.pixel_size = pixel_size
        self.ao_mode = ao_mode
        self.target_separation = target_separation
        self.set_instrument(mode=self.mode, pixel_size=self.pixel_size,
                            ao_mode=self.ao_mode, target_separation=self.target_separation)

        # This could be used to read dictionaries directly from a yml file or a interactively
        self.spectrum = None
        self.source = None
        self.sky = None
        self.setup = None
        self.ao_params = None

    def set_instrument(self, mode, pixel_size, ao_mode, target_separation):
        """
        Reset instrument and perform checks

        """
        if mode.lower() not in ["spectroscopy", "imaging"]:
            raise ValueError
        if pixel_size not in [0.0015, 0.004]:
            raise ValueError
        if ao_mode.lower() not in ["scao", "mcao"]:
            raise ValueError
        if pixel_size == 0.015 and mode.lower() == "spectroscopy":
            raise ValueError
        if mode.lower() == "spectroscopy" and ao_mode.lower() == "mcao":
            raise ValueError
        if target_separation > np.sqrt(2*25**2):   # Only important for SCAO
            raise ValueError

        self.mode = mode
        self.pixel_size = pixel_size
        self.ao_mode = ao_mode
        self.target_separation = target_separation


    def set_sed(self, spectrum_type, **kwargs):
        """

        Parameters
        ----------
        spectrum_type:
        kwargs

        Returns
        -------

        """
        if spectrum_type not in SED_dict.keys():
            raise ValueError("SED not available")

        self.spectrum = dict(spec_name=spectrum_type, params=kwargs)

    def set_target_flux_distribution(self, distribution, **kwargs):
        """
        Same here
        Parameters
        ----------
        distribution
        kwargs

        Returns
        -------

        """
        if distribution not in FLUX_DISTRO_dict.keys():
            raise ValueError("FLUX distribution not availalbe")
        if self.spectrum is None:
            raise ValueError("Please define the SED first")

        self.source = dict(distro_name=distribution, params=kwargs)

    def _get_spectrum(self):
        sed_func = SED_dict[self.spectrum["spec_name"]]
        params = self.spectrum["params"]
        sp = sed_func(**params)

        return sp

    def _get_source(self):
        src_func = FLUX_DISTRO_dict[self.source["distro_name"]]
        params = self.source["params"]

        src = src_func(self._get_spectrum(),
                       **params)

        src.shift(dx=dist, dy=dist)

        return src

    def set_ao(self):
        pass

    def _get_scao_psf(self):
        dist = np.sqrt(0.5 * self.target_separation ** 2)
        psf = scao_psf(x=dist, y=dist)

        return psf


    def set_sky_conditions(self, airmass=1.5, moon_phase=0.5, pwv=10):
        """
        TODO: Check how ScopeSim can read this
        TODO: turbulence and iq affect the PSF!
        """

        self.sky = dict(airmass=airmass, moon_phase=moon_phase, pwv=pwv)

    def set_ao(self, distance, turbulence, iq):

        self.ao_params = dict(disntance=distance, turbulence=turbulence, iq=iq)

    def _get_psf(self):
        kernel = scao_psf()
        psf_effect = sim.effects.psfs.PSF()
        psf_effect.kernel = kernel

        return psf_effect

    def set_setup_obs(self, dit=60, ndit=1, filter_name="Ks"):

        self.setup = dict(dit=dit, ndit=ndit, filter_name=filter_name)

    def run(self, filename):
        """
        Runs the simulation and creates a report

        Returns
        -------
        fits file with the simulation
        report with the results

        """

        src = self._get_source()
        psf_effect = self._get_psf()

        micado = sim.OpticalTrain("MICADO")
        micado.cmds["!OBS.filter_name"] = self.filter_name  # observing filter
        micado.cmds["!INST.pixel_scale"] = self.pixel_size
        micado["armazones_atmo_skycalc_ter_curve"].include = True
        micado["armazones_atmo_default_ter_curve"].include = False
        micado['detector_linearity'].include = False
        micado.cmds["!OBS.dit"] = self.setup_obs["dit"]  # dit & ndit
        micado.cmds["!OBS.ndit"] = self.setup_obs["ndit"]
        micado["relay_psf"].include = False
        micado.optics_manager["default_ro"].add_effect(psf_effect)

        micado.observe(src)
        micado.cmds["!OBS.modes"] = ['SCAO', 'IMG_1.5mas']

        hdus = micado.readout(filename=filename)



