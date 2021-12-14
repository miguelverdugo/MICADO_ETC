import os
import inspect

from .etc import HAWKI_ETC, MICADO_ETC, METIS_ETC

__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))
__data_dir__ = os.path.join(__pkg_dir__, "data")
