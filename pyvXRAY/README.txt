=======
pyvXRAY
=======

*An ABAQUS plug-in for the creation of virtual x-rays from 3D finite element bone models*

Copyright 2013, Michael Hogg

https://pypi.python.org/pypi/pyvXRAY/
https://code.google.com/p/pyvxray/

Author
------

Michael Hogg (michael.christopher.hogg@gmail.com)

Requirements
------------

ABAQUS >= 6.11  (a commerically available finite element software program) 
python  = 2.6   (either built-in to ABAQUS or separate installation)
numpy   = 1.4.0 (need to use same version used in your ABAQUS version)   
PIL    >= 1.1.6

For building cython modules from source (e.g. if not using versions with pre-built modules):
   - Cython >= 0.17
   - A c++ compiler such as gcc or msvc

NOTE: Either gcc or msvc can be used to build the modules if a separate python installation
      is used. However, only msvc (specifically Microsoft Visual Studio 2008) can be used if
      ABAQUS python is used. Typically it is recommended to use MSVC2008 because the version
      of python 2.6 from python.org were built using this compiler.

Support:
--------
Contact author by email (michael.christopher.hogg@gmail.com)

Installation Method
-------------------

The following steps should be followed in order to install the pyvXRAY plug-in:

   1. Download the zip file and extract the contents to a convenient location. 

      There are several versions of pyvXRAY available. This is because pyvXRAY contains Cython 
      modules that must be built before it can be used. The user has the option to:
      
      (i)  download versions with pre-built modules. These are available only for 32-bit
           and 64-bit Windows; and 
   
      (ii) download and build from source. See build instruction below.
      
   2. Setup pyvXRAY plug-in in ABAQUS 

      You can install a plug-in with ABAQUS in several ways. When you start ABAQUS/CAE, it searches
      for plug-in files in the following directories and all their subdirectories:

      - "abaqus_dir\abaqus_plugins", where "abaqus_dir" is the ABAQUS parent directory

      - "home_dir\abaqus_plugins", where "home_dir" is your home directory

      - "current_dir\abaqus_plugins", where "current_dir" is the current directory

      - "plugin_dir", where "plugin_dir" is a directory specified in the abaqus_v6.env file
        by the environment variable "plugin_central_dir". This if often a shared directory on
        a file system that all users can access. For example:

        plugin_central_dir = r'\\fileServer\sharedDirectory'

      I would recommend that, to start with, you use the first option. That is, copy the pyvXRAY 
      directory and files into the abaqus_dir\abaqus_plugins directory. The location of this directory
      depends on the version of ABAQUS being used. For example, the default location on Windows is:

      - ABAQUS 6.11: C:\SIMULIA\Abaqus\6.11-1\abaqus_plugins
      - ABAQUS 6.12: C:\SIMULIA\Abaqus\6.12-1\code\python\lib\abaqus_plugins

      Based on these default directories, the final path to the pyvXRAY directory would be

      - ABAQUS 6.11: C:\SIMULIA\Abaqus\6.11-1\abaqus_plugins\pyvXRAY
      - ABAQUS 6.12: C:\SIMULIA\Abaqus\6.12-1\code\python\lib\abaqus_plugins\pyvXRAY

   3. Build instructions (if installing from source, and not versions with pre-built Cython modules)

      - At command prompt, browse to the pyvXRAY folder.

      - Build using the following command:

        (i)  If using a separate python installation:

                 $ python setup.py build_ext --inplace

             Optionally, the compiler may be specified by using:

                 $ python setup.py build_ext --inplace --compiler=msvc

             for the msvc compiler, or

                 $ python setup.py build_ext --inplace --compiler=mingw32

             for the gcc compiler.

        (ii) If using ABAQUS python:

                 $ abaqus python setup.py build_ext --inplace

             Note that the --compiler option is not available when using ABAQUS python.
 
How to run
----------
To use pyvXRAY, follow these steps:

   1. Open ABAQUS/CAE or ABAQUS/Viewer

   2. Open the ABAQUS odb file in the current viewport

   3. Run pyvXRAY by clicking the following from the ABAQUS toolbar:     

      Plug-ins -> Bone modelling tools -> Create virtual x-rays

   4. Fill in the required information in the pyvXRAY GUI. The required information includes:
      - The part set name of the bone and implant regions
      - The name of the fieldoutput representing bone density
      - A list of steps for which the virtual x-rays will be created (only the last frame in the step is used)
      - Details corresponding to sampling resolution and image size

   5. Click 'OK' as the bottom of the pyvXRAY GUI to run

   6. Look at the Python Scripting Window at the bottom of the ABAQUS GUI for progress and error messages

   7. Analyse the images to investigate regional changes in bone density over time. BMDanalyse 
      (http://pypi.python.org/pypi/BMDanalyse) can be used for this purpose

Documentation
-------------
See website for information on pyvXRAY: http://pypi.python.org/pypi/pyvXRAY
