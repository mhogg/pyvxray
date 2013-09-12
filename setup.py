# -*- coding: utf-8 -*-

# Copyright (C) 2013 Michael Hogg

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
    ext_modules += [ Extension("pyvXRAY.cythonMods", sources=["cython/cythonMods.pyx"],include_dirs=[numpy.get_include()]), ]
    cmdclass.update({ 'build_ext':build_ext })
else:
    ext_modules += [ Extension("pyvXRAY.cythonMods", sources=["cython/cythonMods.c"],include_dirs=[numpy.get_include()]), ]
    
setup(
    name = 'pyvXRAY',
    version = pyvXRAY.__version__,
    description = 'ABAQUS plug-in to create virtual x-rays from bone models',
    license = 'MIT license',
    keywords = ["python", "virtual", "x-rays", "ABAQUS", "finite", "element", "bone"],    
    author = 'Michael Hogg',
    author_email = 'michael.christopher.hogg@gmail.com',
    url = "https://code.google.com/p/pyvxray",
    download_url = "https://code.google.com/p/pyvxray/downloads", 
    packages = ['pyvXRAY'],
    package_data = {'pyvXRAY': ['LICENSE.txt','README.txt','cythonMods.pyd']},
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
About
-----

An ABAQUS plug-in to generate virtual x-rays from 3D finite element bone/implant models.

Creates a time series of virtual x-rays of a bone/implant model by projecting a scalar representing bone density onto three orthogonal planes. This is intended for post-processing of numerical bone remodelling analyses, the main objective of which typically is to detect temporal changes in the bone density around an implant and predict the long term stability of the implant in the absence of clinical data.   

The resulting series of virtual x-ray images are in a format similar to clinical results making the result (1) more easily understood and accepted by clinicians and (2) able to be directly compared to  clinical data, which enables validation of numerical bone remodelling algorithms.
 
The series of images (which are outputted in common image formats such as bmp, jpeg and png) can be analysed using `BMDanalyse <http://pypi.python.org/pypi/BMDanalyse>`_, which was developed for this purpose. This tool enables the quick creation of regions of interest (ROIs) over which changes in bone density over time can be determined.  

Installation
------------

(1) Installation of pyvXRAY with prebuilt extension module

    - Download the binary distribution (zip file) version of pyvXRAY appropriate for your system. This is currently available on the `project website <https://code.google.com/p/pyvxray/>`_ for Windows 64-bit only 

    - Unzip the zip file to a convenient location

    - Copy the abaqus_plugins directory to one of the locations required by ABAQUS. See README.txt or the ABAQUS documentation for further details

(2) Installation of pyvXRAY from source

    - Download source distribution (zip file) version of pyvXRAY

    - Unzip the zip file to a convenient location

    - Open a command window and browse to pyvXRAY folder containing setup.py

    - Build the extension modules using:
      
      ``$ python setup.py build_ext --inplace``
      
      See README.txt for more details on how to do this. Note that this should create file pyvXRAY/cythonMods.pyd (on Windows).
      
    - Copy the entire pyvXRAY to one of the abaqus_plugins directories that are searched by ABAQUS for plug-ins. See README.txt or the ABAQUS documentation for further details  

How to use
----------

To use pyvXRAY, follow these steps:

   1. Open ABAQUS/CAE or ABAQUS/Viewer

   2. Open the ABAQUS odb file in the current viewport

   3. Run pyvXRAY by clicking the following from the ABAQUS toolbar:
      
      ``Plug-ins -> Bone modelling tools -> Create virtual x-rays``

   4. Fill in the required information in the pyvXRAY GUI. The required information includes:

      - The part set name of the bone and implant regions

      - The name of the fieldoutput representing bone density

      - A list of steps for which the virtual x-rays will be created (only the last frame in the step is used)

      - Details corresponding to sampling resolution and image size

   5. Click 'OK' as the bottom of the pyvXRAY GUI to run

   6. Look at the Python Scripting Window at the bottom of the ABAQUS GUI for progress and error messages

   7. Analyse the images to investigate regional changes in bone density over time. `BMDanalyse`_ can be used for this purpose

Requirements
------------

Note that there are a several requirements for using the current version of pyvRAY:

   a. Requires ABAQUS, a commercially available, general purpose finite element analysis software. ABAQUS requires a paid license available from `SIMULIA <http://www.3ds.com/products/simulia>`_ 

   b. The finite element model must consist of tetrahedral elements (ABAQUS element base types C3D4 or C3D10)  

   c. The ABAQUS results file (odb) must contain a scalar variable that represents bone density

""",
)
