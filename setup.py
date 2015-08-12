# -*- coding: utf-8 -*-

# Copyright (C) 2015 Michael Hogg

# This file is part of pyvXRAY - See LICENSE.txt for information on usage and redistribution

import pyvXRAY

from distutils.core import setup
from distutils.extension import Extension
import numpy
try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True
    
cmdclass    = {}
ext_modules = []
if use_cython:  
    ext_modules += [ Extension("pyvXRAY.cythonMods", sources=["cython/cythonMods.pyx"],include_dirs=[numpy.get_include()],language="c++")]
    cmdclass.update({ 'build_ext':build_ext })
else:
    ext_modules += [ Extension("pyvXRAY.cythonMods", sources=["cython/cythonMods.cpp"],include_dirs=[numpy.get_include()],language="c++")]
    
setup(
    name = 'pyvXRAY',
    version = pyvXRAY.__version__,
    description = 'ABAQUS plug-in to create virtual x-rays from 3D finite element bone/implant models',
    license = 'MIT license',
    keywords = ["ABAQUS","plug-in","virtual","x-rays","finite","element","bone","python","cython"],
    author = 'Michael Hogg',
    author_email = 'michael.christopher.hogg@gmail.com',
    url = "https://github.com/mhogg/pyvxray",
    download_url = "https://github.com/mhogg/pyvxray/releases", 
    packages = ['','pyvXRAY'],
    package_data = {'':['LICENSE.txt','README.md'],'pyvXRAY': ['cythonMods.pyd',]},
    classifiers = [
        "Programming Language :: Python",                                  
        "Programming Language :: Cython",         
        "Programming Language :: Python :: 2",             
        "Programming Language :: Python :: 2.6",                                                    
        "Development Status :: 4 - Beta",                                  
        "Environment :: Other Environment", 
        "Environment :: Plugins", 
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",   
        "License :: OSI Approved :: MIT License", 
        "Operating System :: OS Independent",     
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Visualization",
        ],
    ext_modules = ext_modules,
    cmdclass = cmdclass,
    long_description = """

An ABAQUS plug-in to generate virtual x-rays from 3D finite element bone/implant models.

Creates a time series of virtual x-rays of a bone/implant model by projecting a scalar representing bone density onto three orthogonal planes. This is intended for post-processing of numerical bone remodelling analyses, the main objective of which typically is to detect temporal changes in the bone density around an implant and predict the long term stability of the implant in the absence of clinical data.   

The resulting series of virtual x-ray images are in a format similar to clinical results making the result (1) more easily understood and accepted by clinicians and (2) able to be directly compared to  clinical data, which enables validation of numerical bone remodelling algorithms.
 
The series of images (which are outputted in common image formats such as bmp, jpeg and png) can be analysed using `BMDanalyse <http://pypi.python.org/pypi/BMDanalyse>`_, which was developed for this purpose. This tool enables the quick creation of regions of interest (ROIs) over which changes in bone density over time can be determined.
""",
)
