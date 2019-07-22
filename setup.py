from setuptools import setup

setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name='peakpressure',
    url='',
    author='Christopher Howlett',
    author_email='chowlet5@uwo.ca',
    # Needed to actually package something
    packages=['peakpressure'],
    # Needed for dependencies
    install_requires=['numpy','scipy'],
    # *strongly* suggested for sharing
    version='1.0',
    # The license can be anything you like
    license='MIT',
    description='Python package provided functions which are used to estimate peak values from time series of pressures and wind induced internal forces. The functions are based on MATLAB functions written by Joseph A. Main and Fahim Sadek.',
    # We will also need a readme eventually (there will be a warning)
    # long_description=open('README.txt').read(),