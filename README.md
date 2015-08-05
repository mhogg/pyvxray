# pyvXRAY

**An ABAQUS plug-in for the creation of virtual x-rays from 3D finite element bone/implant models.**

**Developed together with [bonemapy](https://github.com/mhogg/bonemapy) and [BMDanalyse](https://github.com/mhogg/BMDanalyse) to provide tools for preparation and post-processing of bone/implant computer models.**

Copyright 2013, Michael Hogg (michael.christopher.hogg@gmail.com)

MIT license - See LICENSE.txt for details on usage and redistribution

## Requirements

### Software requirements

* ABAQUS >= 6.11
* Python Image Library (PIL) >= 1.1.6 OR Pillow >= 2.2.0

For building cython modules from source (e.g. if not using releases with pre-built modules):
* A C compiler. Using ABAQUS Python on Windows requires Microsoft C++. Can use other compilers (i.e. mingw32) if you have a separate Python installation.
* Cython >= 0.17. This is optional, as .c files generated by Cython are provided

**NOTES:**

1.  ABAQUS is a commerical software package and requires a license from [Simulia](http://www.3ds.com/products-services/simulia/overview/)
2.  The authors of pyvXRAY are not associated with ABAQUS/Simulia 
3.  Python and numpy are heavily used by pyvXRAY. These are built in to ABAQUS. All of the last few releases (v6.11 - v6.13) use Python 2.6.x and numpy 1.4.x

### Model setup requirements

* The model must contain only linear or quadrilateral tetrahedral elements (ABAQUS element types C3D4, C3D4H, C3D10, C3D10H, C3D10I, C3D10M, and C3D10MH are all supported).
* The model must have a scalar fieldoutput variable representing bone density. This is typically a state variable such as SDV1. 
  This scalar variable must be available in the last frame of each step to be analysed, as only the last frame is used.

## Installation

####1. Installation of pyvXRAY plug-in

pyvXRAY is an ABAQUS plug-in. ABAQUS plug-ins may be installed in several ways. Only one of the ways is discussed here. For other options the user is referred to the ABAQUS user manuals.

* _Releases with pre-built modules_

  + Download the latest pyvXRAY release with pre-built modules. This is only made available for 32-bit and 64-bit Windows
  + Unpack this to a convenient location
  + Move the `abaqus_plugins\pyvXRAY` folder to the correct location of the `abaqus_plugins` directory within your ABAQUS installation. The location of this directory depends on your ABAQUS version. Some possible locations are:

      v6.11-x: `C:\SIMULIA\Abaqus\6.11-x\abaqus_plugins`

      v6.12-x: `C:\SIMULIA\Abaqus\6.12-x\code\python\lib\abaqus_plugins`

      v6.13-x: `C:\SIMULIA\Abaqus\6.13-x\code\python\lib\abaqus_plugins`

* _Installation from source_

  + Download the latest pyvXRAY source, typically called `pyvXRAY-x.x.x.zip` or `pyvXRAY-x.x.x.tar.gz`, where `x.x.x` is the version number
  + Unpack this to a convenient location
  + Open a command prompt and browse to directory `pyvXRAY-x.x.x` (containing file `setup.py`)
  + Run the following command:

            abaqus python setup.py build_ext --inplace

      which will build the Cython modules. If Cython is available, it will be used. Otherwise the .c files previously generated using Cython will be compiled directly.
  + Copy the pyvXRAY sub-folder to the `abaqus_plugins` directory within your ABAQUS installation, following the instructions above for pre-built distribution 

####2. Installation of pyvXRAY dependencies

The ABAQUS GUI is built on Python, and has its own Python installation. This Python installation is not the typically Python setup, so some guidance is provided here on how to install pyvXRAY's dependencies.

Currently pyvXRAY has only one dependency that is not part of the ABAQUS Python, which is PIL / Pillow (NOTE: Pillow is a PIL fork that appears to have largely superseded PIL). 
On Windows it is easiest to install PIL / Pillow using a binary installer, particularly because PIL / Pillow have many dependencies. There are currently several issues with this:

* ABAQUS Python is typically not registered in the Windows registry, and therefore installation with binary installers will not work by default because the ABAQUS Python
  installation does not appear in the list of available Python installations

* PIL binaries are readily available for 32-bit Windows (at [pythonware.com](http://www.pythonware.com/products/pil/)) 

* Pillow binaries for Python 2.6 on both 32-bit and 64-bit Windows are available on [PyPi](https://pypi.python.org/pypi/Pillow)

Given these limitations, there are two obvious choices for installating PIL / Pillow with ABAQUS Python.

* _Use a separate Python installation and binary installer_

  The SIMULIA support site suggests that a separate Python installation be setup. The dependencies can be installed easily into this separate Python installation, after which they can 
  then be used by ABAQUS Python. This separate Python installation version must match the ABAQUS Python version, which has been 2.6.x for the following few ABAQUS versions. 
  This method is the best solution if you do not feel comfortable modifying the Windows registry (as recommended next).

  + Download and install Python 2.6 (i.e. the plain vanilla version from www.python.org/download is recommended)
  + Download and run the corresponding binary installer for PIL / Pillow. Pillow is in active development and is recommended over PIL.
  + Create an environment variable `PYTHONPATH=C:\Python26\Lib\site-packages` that tells ABAQUS Python where these packages are installed, assuming that `C:\Python26` is the installation directory.

* _Edit the Windows registry and use a binary installer (Windows only)_

  By editing the Windows registry the binary installers will be able to detect the ABAQUS Python version and install as usual. Use caution when editing the Windows registry or backup your 
  registry before hand.

  You need to create registry entry `HKEY_LOCAL_MACHINE\Software\Python\Pythoncore\2.6\InstallPath` and make its (Default) value that of your ABAQUS Python directory location. Registry key 
  `HKEY_CURRENT_USER` also works. 
  This location depends on the ABAQUS version. For the default ABAQUS installation location, possible locations are:

      v6.11-x: `C:\\SIMULIA\\Abaqus\\6.11-x\\External\\Python`
     
      v6.12-x: `C:\\SIMULIA\\Abaqus\\6.12-x\\tools\\SMApy`
     
      v6.13-x: `C:\\SIMULIA\\Abaqus\\6.13-x\\tools\\SMApy\\python2.6`

  Editing the Windows registry can be done using the regedit utility.

## Usage

* Open ABAQUS/CAE
* Open an odb file
* To launch the pyvXRAY GUI, go to the menubar at the top of the screen and select:

        Plug-ins --> pyvXRAY --> Create virtual x-rays

* Complete the required inputs in the GUI to suit the current model. More information is given below about the inputs
* Click OK to run pyvXRAY
* Look at the message area at the bottom of the screen for messages. On completion 'Finished' will be shown.

## Required inputs

A basic description of each of the inputs required by pyvXRAY is listed here.

<table>
<th align="left">GUI tab</th><th>Input name </th><th>Input description</th>

<tr>
<td width="100">Select regions</td>
<td>Result file: Odb</td>
<td>The ABAQUS result file</td>
</tr>

<tr>
<td></td>
<td width="150">Bone region: Bone set</td>
<td>The name of the element set representing the bone</td>
</tr>

<tr>
<td></td>
<td>Bone region: Density variable</td>
<td>A scalar fieldoutput variable representing bone density.<br>This is most often a state variable i.e. SDV1</td>
</tr>

<tr>
<td></td>
<td>Implant region: Show implant on x-rays</td>
<td>Option to include implant on the virtual x-rays </td>
</tr>

<tr>
<td></td>
<td>Implant region: Implant set</td>
<td>The name of the element set representing the implant</td>
</tr>

<tr>
<td></td>
<td>Implant region: Density (kg/m^3)</td>
<td>The density of the implant material in kg/m^3 i.e. 4500 for Titanium Alloy</td>
</tr>

<tr>
<td>Inputs</td>
<td>Required inputs: Step list</td>
<td>A list of steps to be analysed i.e. 1, 2, 3. A virtual x-ray is created for the last frame of each step in this list.</td>
</tr>

<tr>
<td></td>
<td>Required inputs: Coordinate system</td>
<td>The name of the coordinate system used to create the projections. By default this is the global coordinate system. However, the views can be changed by creating a new coordinate 
system in ABAQUS and using it instead.</td>
</tr>

<tr>
<td></td>
<td>Required inputs: Mapping resolution (mm)</td>
<td>pyvXRAY works by mapping the results of the bone density variable onto a regular grid. The mapping resolution is the cell spacing of this regular grid. Decreasing this number 
increases the accuracy of the mapping, but also increases the calculation time. As a first pass, a value of around 2mm is recommended to ensure that output is as expected.</td>
</tr>

<tr>
<td>X-ray settings</td>
<td>Settings: Base name of xray file(s)</td>
<td>This is the base or root name of the virtual x-ray image files. That is, image files are labelled <code>basename_projection_stepnumber</code> i.e. <code>basename_XY_1</code> for 
the X-Y projection from Step 1.</td>
</tr>

<tr>
<td></td>
<td>Settings: Approx size of x-ray images</td>
<td>Resizing of images is performed to make the number of pixels along the largest image dimension equal to this value.</td>
</tr>

<tr>
<td></td>
<td>Settings: Image file format</td>
<td>Output format of images. Options are bmp, jpeg and png.</td>
</tr>

<tr>
<td></td>
<td>Settings: Smooth images</td>
<td>Turn on image smoothing. <code>PIL.ImageFilter.SMOOTH</code> is used to perform the smoothing.</td>
</tr>

<tr>
<td></td>
<td>Settings: Manual scaling of images</td>
<td>pyvXRAY scales the mapped bone density values when creating the virtual x-ray images. The image files are 24-bit (or 8-bit per channel), so the gray scale range is essentially 0-255. 
The scale factor used ensures that this range is fully utilised and that none of the images in the series are over-exposed. Activating this option reports the scale factors used and gives 
the user the ability to change these values. This may be desirable when comparing virtual x-rays from different models; an equal comparison is possible only if the same scale factors are 
used for both. </td>
</tr>

</table>

## Outputs

pyvXRAY outputs a series of virtual x-rays correponding to the bone density results in a list of specified analysis steps. The bone density is mapped from the Finite Element Model to a 
overlapping regular grid of points and then projected onto each of the three Cartesian coordinate planes. If the model has an implant, then this can also be shown. The virtual x-ray images
are saved in common image formats (bmp, jpeg, and png) and can be opened in any graphics package. These images can then be analysed to determine changes in the grey scale values, which
can be related to the change in Bone Mineral Density (BMD) over time.

The recommended package for analysing these images is [BMDanalyse](https://github.com/mhogg/BMDanalyse), which is available free under the MIT license. BMDanalyse can be used to create 
regions of interest (ROIs) and determine the change in the average grey scale value within each ROI for all images in the series.

## Help
 
For help create an Issue or a Pull Request on Github.
