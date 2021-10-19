# -*- coding: utf-8 -*-
"""
    Setup file for MICADO ETC.
"""
import sys
import setuptools
from datetime import datetime


# Version number
MAJOR = 0
MINOR = 1
ATTR = '0'

VERSION = '%d.%d%s' % (MAJOR, MINOR, ATTR)


def write_version(filename='micado_etc/version.py'):
    """Write a file version.py"""
    cnt = """# THIS FILE GENERATED BY speXtra setup.py 
version = '{}' 
date    = '{}' 
     """
    timestamp = datetime.utcnow().strftime('%Y-%m-%d %T GMT')
    with open(filename, 'w', encoding='utf-8') as fd:
        if sys.version_info.major == 2:
            fd.write(cnt.format(VERSION, timestamp).decode('utf-8'))
        else:
            fd.write(cnt.format(VERSION, timestamp))


with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()


def setup_sp():

    write_version()
    setuptools.setup(
        version=VERSION,
        name="micado_etc",
        descriptio="MICADO ETC",
        long_description=long_description,
        author="Miguel Verdugo",
        author_email="mverduol@gmail.com",
        license="MIT",
        url="https://github.com/miguelverdugo/MICADO_ETC",
        package_dir={'micado_etc': 'micado_etc'},
        packages=['micado_etc'],
        package_data={'micado_etc': ['micado_etc/data/*']},
        include_package_data=True,
        install_requires=["numpy",
                          "astropy>=3.1",
                          "synphot>=0.2.0",
                          "pyyaml", ],   # Also ScopeSim environment
        classifiers=["Programming Language :: Python :: 3",
                     "License :: OSI Approved :: MIT License",
                     "Operating System :: OS Independent",
                     "Intended Audience :: Science/Research",
                     "Topic :: Scientific/Engineering :: Astronomy", ],
    )


if __name__ == "__main__":
    setup_sp()
