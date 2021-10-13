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

    def test_select_sed(self):
        etc = ETC()
        etc.set_sed("template")
        assert isinstance(etc.spectrum, Spextrum)


