import pytest
import inspect
import os

import numpy as np
import astropy.units as u
from spextra import Spextrum

from ..micado_etc import ETC


class TestETCBasic:

    def test_initialization(self):
        etc = ETC()
        assert isinstance(etc, ETC)

    def test_initialization_ii(self):
        etc = ETC(mode="imaging", pixel_size=0.0015, ao_mode="MCAO")
        assert etc.pixel_size == 0.0015

    def test_set_instrument(self):
        etc = ETC(mode="spectroscopy", pixel_size=0.004, ao_mode="SCAO")
        assert etc.mode =="spectroscopy"


    def test_set_sed(self):
        etc = ETC()
        params = dict(template_name="kc96/s0", amplitude=18, filter_curve="R", redshift=1)
        etc.set_sed("template", **params)
        assert etc.spectrum["params"].keys() == params.keys()

    def test_get_sed(self):
        etc = ETC()
        etc.set_sed("template", template_name="kc96/s0")
        print(etc.spectrum)
        assert isinstance(etc._get_spectrum(), Spextrum)


