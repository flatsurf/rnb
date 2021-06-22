import os
from distutils.core import setup

setup(
    name = 'rnb',
    version = '0.1.0',
    packages = ['rnb',],
    license = 'GPL 3.0+',
    long_description = open('README.md').read(),
    include_package_data=True,
)
