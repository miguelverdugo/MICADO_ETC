import numpy as np

import astropy.units as u

from spextra import Spextrum, Passband

import scopesim as sim

from .utils import *

#from scopesim_templates.basic.basic import empty_sky
#from scopesim_templates.b
"""
Expected usage:

initialize

etc = ETC()



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



#class ETC:
#    common class for setting parameters
#class ETC_HAWKI(ETC)
#class ETC_MICADO(ETC)
#class ETC_METIS(ETC)



class ETC_base:
    def __init__(self, mode, pixel_size, ao_mode=None, target_separation=0):

        self.mode = mode
        self.pixel_size = pixel_size
        self.ao_mode = ao_mode
        self.target_separation = target_separation

        self.sed_params = None
        self.source_params = None
        self.ao_params = None
        self.sky_params = None
        self.obs_params = None

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

        self.sed_params = dict(spec_name=spectrum_type, params=kwargs)

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

        self.source_params = dict(distro_name=distribution, params=kwargs)

    def set_sky_conditions(self, airmass=1.5, moon_phase=0.5, pwv=10):
        """
        TODO: Check how ScopeSim can read this
        TODO: turbulence and iq affect the PSF!
        """

        self.sky_params = dict(airmass=airmass,
                        moon_phase=moon_phase,
                        pwv=pwv)

    def set_ao(self, profile_name="EsoQ4", zenDist=0, seeing=1):

        self.ao_params = dict(profile_name=profile_name,
                              zenDist=zenDist,
                              #turbulence=turbulence,
                              #iq=iq,
                              seeing=seeing)

    def set_setup_obs(self, dit=60, ndit=1, filter_name="Ks"):

        self.obs_params = dict(dit=dit, ndit=ndit, filter_name=filter_name)

    def _get_spectrum(self):
        sed_func = SED_dict[self.spectrum["spec_name"]]
        params = self.spectrum["params"]
        sp = sed_func(**params)

        return sp

    def _get_source(self):
        src_func = FLUX_DISTRO_dict[self.source["distro_name"]]
        params = self.source["params"]

        params.update({"sed": self._get_spectrum(),
                       "filter_name": self.spectrum["params"]["filter_curve"],
                       "magnitude": self.spectrum["params"]["magnitude"]
                       })

        src = src_func(**params)
        src.shift(dx=self.dx,
                  dy=self.dy)

        return src


class HAWKI_ETC(ETC_base):

    def __init__(self, mode="imaging", pixel_size=0.106, ao_mode="no_ao" ):

        super(HAWKI_ETC, self).__init__(mode=mode, pixel_size=pixel_size, ao_mode=ao_mode)

        self.set_instrument()

    def set_instrument(self, ao_mode):   # This is very instrument specific
        if ao_mode.lower() not in ["no_ao", "ao"]:
            raise ValueError







class MICADO_ETC(ETC_base):

    def __init__(self, mode="imaging", pixel_size=0.004, ao_mode="SCAO", target_separation=0, dx=0, dy=0):



        self.mode = mode
        self.pixel_size = pixel_size
        self.image_mode = "IMG_4mas"
        self.ao_mode = ao_mode
        self.target_separation = target_separation
        self.dx = dx
        self.dy = dy
        self.set_instrument(mode=self.mode, pixel_size=self.pixel_size,
                            ao_mode=self.ao_mode, target_separation=self.target_separation,
                            dx=self.dx, dy=self.dy)

        # This could be used to read dictionaries directly from a yml file or a interactively
        self.spectrum = None
        self.source = None
        self.sky = None
        self.setup = None
        self.set_ao()

    def set_instrument(self, mode, pixel_size, ao_mode, target_separation, dx, dy):
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
        if target_separation > 0 and dx == 0 and dy == 0:
            dist = np.sqrt(0.5 * self.target_separation ** 2)
            self.dx = dist
            self.dy = dist
        if pixel_size == 0.015:
            self.image_mode = 'IMG_1.5mas'


        self.mode = mode
        self.pixel_size = pixel_size
        self.ao_mode = ao_mode
        self.target_separation = target_separation






    def _get_psf(self):

        if self.ao_mode.lower() == "scao":
            kwargs = self.ao_params
            kwargs.update({"x":self.dx, "y":self.dy})
            kernel = scao_psf(**kwargs)
            psf_effect = sim.effects.psfs.PSF()
            psf_effect.kernel = kernel
        else:
            filename = sim.utils.find_file("PSF_MCAO_ConstPSF_40_18_6.fits")
            psf_effect = sim.effects.psfs.FieldConstantPSF(filename=filename)

        return psf_effect



    def run(self, filename=None):
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
        micado.cmds["!OBS.filter_name"] = self.setup["filter_name"]  # observing filter
        micado.cmds["!INST.pixel_scale"] = self.pixel_size
        micado.cmds["!OBS.modes"] = [self.ao_mode.upper(), self.image_mode]
        micado.cmds["!OBS.dit"] = self.setup["dit"]  # dit & ndit
        micado.cmds["!OBS.ndit"] = self.setup["ndit"]
        micado["armazones_atmo_skycalc_ter_curve"].include = True
        micado["armazones_atmo_default_ter_curve"].include = False
        micado['detector_linearity'].include = False
        micado.cmds["!DET.y"] = self.dy / self.pixel_size
        micado.cmds["!DET.x"] = self.dx / self.pixel_size
        micado["relay_psf"].include = False
        micado.optics_manager["default_ro"].add_effect(psf_effect)

        micado.observe(src)
        noiseless_obj = micado.image_planes[0].data
    #    signal_to_noise = np.sqrt(noiseless_image)
        hdu_obj = micado.readout(filename=filename)
        observed_image = hdu_obj[0][1]

        sky_src = empty_sky()
        micado.observe(sky_src)
        noiseless_sky = micado.image_planes[0].data
        hdu_sky = micado.readout(filename=filename)
       # signal_to_noise = np.sqrt(hdu_sky[0][1])


        signal_to_noise = np.sqrt(self.setup["ndit"]) * noiseless_obj / (noiseless_obj + noiseless_sky)

        return observed_image, signal_to_noise
  #      return noiseless_obj, noisless_sky



