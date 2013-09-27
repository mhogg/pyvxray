# pyvXRAY

**An ABAQUS plug-in for the creation of virtual x-rays from 3D finite element bone/implant models**

**Developed together with [bonemapy](https://github.com/mhogg/bonemapy) and [BMDanalyse](https://github.com/mhogg/BMDanalyse) to provide tools for preparation and post-processing of bone/implant computer models.**

Copyright 2013, Michael Hogg (michael.christopher.hogg@gmail.com)

MIT license - See pyvxray/LICENSE.txt for details on usage and distribution

## Requirements

### Software requirements

* ABAQUS >= 6.11
* Python Image Library (PIL) >= 1.1.6

For building cython modules from source (e.g. if not using releases with pre-built modules):
* A C compiler. Using ABAQUS Python on Windows requires Microsoft C++. Can use other compilers if you have a separate Python installation.
* Cython >= 0.17. This is optional, as .c files generated by Cython are provided

**NOTES:**

1.  ABAQUS is a commerical software package and requires a license from [Simulia](http://www.3ds.com/products-services/simulia/overview/)
2.  The authors of pyvXRAY are not associated with ABAQUS/Simulia 
3.  Python and numpy are heavily used by pyvXRAY. These are built in to ABAQUS. All of the last few releases (v6.11-v6.13) use Python 2.6 and numpy 1.4

### Model setup requirements

* The model must contain only linear or quadrilateral tetrahedral elements (ABAQUS element types C3D4 and C3D10, respectively)
* The model must have a scalar fieldoutput variable representing bone density. This is typically a state variable such as SDV1
* This scalar variable must be available in the last frame of each step to be analysed, as only the last frame is used.

**LIMITATIONS:**

* Virtual x-rays can only be created for element sets within part instances. Assembly element sets (which may contain elements from more than a single part) are currently not supported

## Installation

pyvXRAY is an ABAQUS plug-in. ABAQUS plug-ins may be installed in several ways. Only one of the ways is discussed here. For other options the user is referred to the ABAQUS user manuals.

The ABAQUS GUI is built on Python, and has its own Python installation. This Python installation is not the typically Python setup, so some guidance is provided here on how to install pyvXRAY's dependencies.

####1. Installation of pyvXRAY plug-in

* _Releases with pre-built modules_

  + Download the latest pyvXRAY release with pre-built modules. This is only made available for a single platform, which is 64-bit Windows
  + Unpack this to a convenient location
  + Copy the pyvXRAY sub-folder to the `abaqus_plugins` directory within your ABAQUS installation. The location of this directory depends on your ABAQUS version. Some possible locations are:

      v6.11-1: `C:\SIMULIA\Abaqus\6.11-1\abaqus_plugins`

      v6.12-1: `C:\SIMULIA\Abaqus\6.12-1\code\python\lib\abaqus_plugins`

      v6.13-1: `C:\SIMULIA\Abaqus\6.13-1\code\python\lib\abaqus_plugins`

* _Installation from source_

  + Download the latest pyvXRAY source, typically called `pyvXRAY-x.x.x.zip` or `pyvXRAY-x.x.x.tar.gz`, where x.x.x is the version number
  + Unpack this to a convenient location
  + Open a command prompt and browse to directory pyvXRAY-x.x.x (containing file setup.py)
  + Run the following command:

            abaqus python setup.py build_ext --inplace

      which will build the cython modules
  + Copy the pyvXRAY sub-folder to the `abaqus_plugins` directory within your ABAQUS installation, following the instructions above for pre-built distribution 

####2. Installation of pyvXRAY dependencies

Currently pyvXRAY has only one dependency that is not part of the ABAQUS Python, which is PIL. On Windows it is easist to install PIL using a binary installer. However, binary installers detect and list the available Python installations in which to install the Python package. ABAQUS Python does not appear in the list of available Python installations, so installation will not work. There are several solutions to this:

* _Use a separate Python installation and binary installer_

  The SIMULIA support site suggests that a separate Python installation be setup. The dependencies can be installed easily into this separate Python installation, which can then be used by ABAQUS Python. This separate Python installation version must match the ABAQUS Python version, which has been 2.6.x for the following few ABAQUS versions. This method is often the easist solution if you need to use several Python packages within ABAQUS, i.e. PIL, scipy etc., especially if you have an older version of EPD or similar Python distibution that comes with all these packages.

  + Download and install Python 2.6 (standard python, or EPD distribution etc)
  + If your Python distribution doesn't come with PIL, then download and run the PIL binary installer. (Note: Binary installers for Windows 64-bit can be downloaded from [here](http://www.lfd.uci.edu/~gohlke/pythonlibs/)).
  + Create an environment variable `PYTHONPATH=C:\Python26\Lib\site-packages` that tells ABAQUS Python where these packages are installed, assuming that `C:\Python26` is the installation directory.

* _Install PIL from source (on Windows this requires Microsoft C++ to be installed)_

  + Download the PIL source, typically called `Imaging-x.x.x.tar.gz` where x.x.x is the version number
  + Unpack this to a convenient location
  + Open a command prompt and browse to folder `Imaging-x.x.x` (containing the setup.py file)
  + At the command prompt enter:

            abaqus python setup.py install

* _Edit the Windows registry and use a binary installer (Windows only)_

  By editing the Windows registry the binary installers will be able to detect the ABAQUS Python version and install as usual. Please use caution when editing the Windows registry or backup your registry before hand.

  You need to create registry entry `HKEY_CURRENT_USER\Software\Python\Pythoncore\2.6\InstallPath` and make its value that of your ABAQUS Python directory location. This location depends on the ABAQUS version. For the default ABAQUS installation location, possible locations are:

      v6.11-1: `C:\\SIMULIA\\Abaqus\\6.11-1\\External\\Python`
     
      v6.12-1: `C:\\SIMULIA\\Abaqus\\6.12-1\\tools\\SMApy`
     
      v6.13-1: `C:\\SIMULIA\\Abaqus\\6.13-1\\tools\\SMApy\\python2.6`

  Editing the Windows registry can be done using the regedit utility.

## Usage

* Open ABAQUS/CAE
* Open the ABAQUS odb file within the current viewport
* To launch the pyvXRAY GUI, go to the menubar at the top of the screen and select:

        Plug-ins --> pyvXRAY --> Create virtual x-rays

* Complete the required inputs in the GUI to suit the current model. More information is given below about the inputs. Error messages will be issued if any information is entered incorrectly.
* Click OK to run pyvXRAY
* Look at the message area at the bottom of the screen for messages. On completion 'Finished' will be shown.

## Required inputs

A basic description of each of the inputs required by pyvXRAY is listed here.

| GUI tab           | Input name                  | Input description     
| ----------------- | --------------------------- | ----------------- 
| Select regions    | Bone part instance          | The name of the part instance containing the bone                                                                                                                                        
|                   | Bone Set name               | The element set in the part instance representing bone. If the entire part instance is bone, then an element set containing all the elements in the part instance is needed.             
|                   | Bone Density variable       | A scalar fieldoutput variable representing bone density. This is most often a state variable i.e. SDV1                                                                                   
|                   | Show implant on x-rays      | Option to include implant on the virtual x-rays                                                                                                                                                              
|                   | Implant part instance name  | The name of the part instance containing the implant                                                                                                                                     
|                   | Implant set name            | The element set in the part instance representing the implant. If the entire part instance is an implant, then an element set containing all the elements in the part instance is needed.
|                   | Implant density (kg/m^3)    | The density of the implant material in kg/m^3 i.e. 4500 for Titanium Alloy                                                                                                               
|                   |                             |                       
| Inputs            | Step list                   | A list of steps to be analysed i.e. 1,2,3. A virtual x-ray is created for the last frame of each step in this list. If only a single step is required, must still contain a comma.
|                   | Coordinate system           | The name of the coordinate system used to create the projections. By default this is the global coordinate system. However, the views can be changed by creating a new coordinate system in ABAQUS and using it instead.
|                   | Mapping resolution (mm)     | pyvXRAY works by mapping the results of the bone density variable onto a regular grid. The mapping resolution is the cell spacing of this regular grid. Decreasing this number increases the accuracy of the mapping, but also increases the calculation time. As a first pass, a value of around 2mm is recommended to ensure that all the inputs are correct.
|                   |                             |                                                                                                                                                                                                                                                                                   
| X-ray settings    | Base name of xray file(s)   | This is the base or root name of the virtual x-ray image files. These are labelled `basename_stepnumber_projection` i.e. `basename_1_XY` for Step 1 projected onto the X-Y plane
|                   | Approx size of x-ray images | Some scaling of images is performed to make the number of pixels along the largest image dimension equal to this value
|                   | Image file format           | Output format of images. Options are png, jpeg, bmp 
|                   | Smooth images               | Turn on image smoothing
|                   | Manual scaling of images    | pyvXRAY scales the mapped bone density values when creating the virtual x-ray images. The image files are 24-bit (or 8-bit for each RGB channel), so the gray scale range is essentially 0-255. The scale factor used ensures that this range is fully utilised and that none of the images in the series are over-exposed. Activating this option reports the scale factors used and gives the user the ability to change these values. This may be desirable when comparing virtual x-rays from different models; an equal comparison is possible only if the same scale factors are used for both. 

<table>
<th align="left">GUI tab</th><th>Input name </th><th>Input description</th>
<tr>
<td width="75">Select regions</td>
<td width="200">Bone: Part instance name</td>
<td>The name of the part instance containing the bone</td>
</tr>
<tr>
<td></td>
<td>Bone: Set name</td>
<td>The element set in the part instance representing bone. If the entire part instance is bone, then an element set containing all the elements in the part instance is needed.</td></tr>
<tr>
<td></td>
<td>Bone: Density variable</td>
<td>A scalar fieldoutput variable representing bone density.<br>This is most often a state variable i.e. SDV1</td>
</tr>
<tr>
<td></td>
<td>Show implant on x-rays</td>
<td>Option to include implant on the virtual x-rays </td>
</tr>
<tr>
<td></td>
<td>Implant: Part instance name</td>
<td>The name of the part instance containing the implant</td>
</tr>
<tr>
<td></td>
<td>Implant: Set name</td>
<td>The element set in the part instance representing the implant. If the entire part instance is an implant, then an element set containing all the elements in the part instance is needed.</td>
</tr>
<tr>
<td></td>
<td>Implant: Density (kg/m^3)</td>
<td>The density of the implant material in kg/m^3 i.e. 4500 for Titanium Alloy</td>
</tr>
<tr>
<td>Inputs</td>
<td>Step list</td>
<td>A list of steps to be analysed i.e. 1,2,3. A virtual x-ray is created for the last frame of each step in this list. If only a single step is required, must still contain a comma.</td>
</tr>
<tr>
<td></td>
<td>Coordinate system</td>
<td>The name of the coordinate system used to create the projections. By default this is the global coordinate system. However, the views can be changed by creating a new coordinate system in ABAQUS and using it instead.</td>
</tr>
<tr>
<td></td>
<td>Mapping resolution (mm)</td>
<td>pyvXRAY works by mapping the results of the bone density variable onto a regular grid. The mapping resolution is the cell spacing of this regular grid. Decreasing this number increases the accuracy of the mapping, but also increases the calculation time. As a first pass, a value of around 2mm is recommended to ensure that all the inputs are correct.</td>
</tr>
</table>

## Outputs

pyvXRAY outputs a series of virtual x-rays correponding to the analysis steps in the ABAQUS odb file. The virtual x-ray images are saved in common image formats (png, jpeg, bmp) and can be opened in any imaging package.

The recommended package for analysing these images is [BMDanalyse](https://github.com/mhogg/BMDanalyse), which is available free under the MIT license. BMDanalyse can be used to create regions of interest (ROIs) and output a graph of the change in average grey scale value in each ROI in the image series. This can be related to the change in time of Bone Mineral Density (BMD).

## Help
 
For help post a question on the [project support page](https://groups.google.com/forum/#!forum/pyvxray) or create an Issue on Github.
