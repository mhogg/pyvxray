# -*- coding: utf-8 -*-

# Copyright (C) 2013 Michael Hogg

# This file is part of pyvXRAY - See LICENSE.txt for information on usage and redistribution

from abaqusGui import *
from abaqusConstants import ALL
import os
from version import version as __version__

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class PyvXRAY_plugin(AFXForm):

    def __init__(self, owner):
        
        # Construct the base class 
        AFXForm.__init__(self, owner)
        
        self.imageFormats=['bmp','jpeg','png']
        
        # Keyword definitions
        self.radioButtonGroups   = {}
        self.cmd                 = AFXGuiCommand(mode=self, method='createVirtualXrays',objectName='virtualXrays', registerQuery=False)
        self.bPartNameKw         = AFXStringKeyword(self.cmd, 'bPartName', True, 'PART-1-1')
        self.bSetNameKw          = AFXStringKeyword(self.cmd, 'bSetName', True, 'BONE')
        self.BMDfonameKw         = AFXStringKeyword(self.cmd, 'BMDfoname', True, 'SDV1')
        self.showImplantKw       = AFXBoolKeyword(self.cmd,   'showImplant', AFXBoolKeyword.TRUE_FALSE, True, False)
        self.iPartNameKw         = AFXStringKeyword(self.cmd, 'iPartName', True, 'PART-1-1')
        self.iSetNameKw          = AFXStringKeyword(self.cmd, 'iSetName', True, 'IMPLANT')
        self.iDensityKw          = AFXFloatKeyword(self.cmd,  'iDensity', True, 4500)
        self.stepListKw          = AFXStringKeyword(self.cmd, 'stepList', True, '1')
        self.csysNameKw          = AFXStringKeyword(self.cmd, 'csysName', True, 'CSYS-1')
        self.resGridKw           = AFXFloatKeyword(self.cmd,  'resGrid', True, 2)
        self.imageNameBaseKw     = AFXStringKeyword(self.cmd, 'imageNameBase', True, 'vxray')
        self.preferredXraySizeKw = AFXIntKeyword(self.cmd,    'preferredXraySize', True, 800)
        self.imageFormatKw       = AFXStringKeyword(self.cmd, 'imageFormat', True, self.imageFormats[-1])
        self.smoothKw            = AFXBoolKeyword(self.cmd,   'smooth', AFXBoolKeyword.TRUE_FALSE, True, True)
        self.manualScalingKw     = AFXBoolKeyword(self.cmd,   'manualImageScaling', AFXBoolKeyword.TRUE_FALSE, True, False)
        self.showImplant         = None

    def getFirstDialog(self):

        import pyvXRAYDB
        return pyvXRAYDB.PyvXRAYDB(self)

    def doCustomChecks(self):
    
        # Check that object in current viewport is an odb file
        displayedType = getDisplayedObjectType()
        if displayedType!=ODB:
            showAFXErrorDialog(self.getCurrentDialog(), 'Error: Object in current viewport is not an odb object')
            return False    
        odb = session.viewports[session.currentViewportName].displayedObject  

        # Check that the selected bone region and element set exists
        bPartName = self.bPartNameKw.getValue()
        bSetName  = self.bSetNameKw.getValue()
        if bPartName not in odb.rootAssembly.instances.keys():
            showAFXErrorDialog(self.getCurrentDialog(), 'Error: %s is not a part instance in the current odb' % bPartName)
            return False
        if bSetName not in odb.rootAssembly.instances[bPartName].elementSets.keys():
            showAFXErrorDialog(self.getCurrentDialog(), 'Error: %s is not an element set in part instance %s' % (bSetName,bPartName))
            return False
            
        # If user has requested the implant to be displayed, then check the selected implant region and element set exists
        if self.showImplant:
            iPartName = self.iPartNameKw.getValue()
            iSetName  = self.iSetNameKw.getValue()
            iDensity  = self.iDensityKw.getValue()
            if iPartName not in odb.rootAssembly.instances.keys():
                showAFXErrorDialog(self.getCurrentDialog(), 'Error: %s is not a part instance in the current odb' % iPartName)
                return False 
            if iSetName not in odb.rootAssembly.instances[iPartName].elementSets.keys():
                showAFXErrorDialog(self.getCurrentDialog(), 'Error: %s is not an element set in part instance %s' % (iSetName,iPartName))
                return False
            try: float(iDensity)
            except: 
                showAFXErrorDialog(self.getCurrentDialog(), 'Error: Implant density value is not a number')
                return False               
            if iDensity<0:
                showAFXErrorDialog(self.getCurrentDialog(), 'Error: Implant density must be greater than 0')
                return False    

        # Check that values in stepList are valid. Also check that the values are in increasing order
        stepList = self.stepListKw.getValue()
        try: 
            stepList = stepList.split(',')
            stepList = [int(s) for s in stepList]
        except: 
            showAFXErrorDialog(self.getCurrentDialog(), 'Error: Cannot convert step list values to integers')
            return False
        diff = [True for i in range(len(stepList)-1) if stepList[i+1]-stepList[i]<0]
        if len(diff)>0: 
            showAFXErrorDialog(self.getCurrentDialog(), 'Error: Step numbers in step list not in increasing order' )
            return False 

        # Check that coordinate system exists
        # Do this is in the kernel - If it doesn't exist, we just use the global coordinate system. A message will be written from the kernel to the Message Area.   
        # NOTE: This should really be replaced with a list of available CSYSs, so the user knows what can be selected.        
    
        # Check mapping resolution, resGrid
        resGrid = self.resGridKw.getValue()
        try: resGrid = float(resGrid)
        except: 
            showAFXErrorDialog(self.getCurrentDialog(), 'Error: "Inputs: Mapping resolution" value not valid')
            return False 

        # Check preferred size of images
        preferredXraySize = self.preferredXraySizeKw.getValue()
        try: preferredXraySize = int(preferredXraySize)
        except: 
            showAFXErrorDialog(self.getCurrentDialog(), 'Error: "X-ray Settings: Approx size of x-ray" value not valid')
            return False
        minXraySize = 100
        if preferredXraySize < minXraySize:
            showAFXErrorDialog(self.getCurrentDialog(), 'Error: Minimum virtual x-ray image size is %i pixels' % minXraySize)
            return False   

        # Check that all steps in step list exist and that density variable exists in all steps (last frame)
        # NOTE: Do this last because it is the most time consuming
        BMDfoname = self.BMDfonameKw.getValue()
        stepInfo  = {}
        for stepName,step in odb.steps.items():
            stepInfo[step.number] = step.frames[-1].fieldOutputs.keys()
        for stepNumber in stepList:
            if stepNumber not in stepInfo:
                showAFXErrorDialog(self.getCurrentDialog(), 'Error: Step number %i is not available in odb' % stepNumber)
                return False
            if BMDfoname not in stepInfo[stepNumber]:
                showAFXErrorDialog(self.getCurrentDialog(), 'Error: Density variable %s is not available in Step number %i' % (BMDfoname,stepNumber))
                return False
                
        # Check for Abaqus version >= 6.11 
        majorNumber, minorNumber, updateNumber = getAFXApp().getVersionNumbers()
        if majorNumber==6 and minorNumber < 11:    
            showAFXErrorDialog( self.getCurrentDialog(), 'Error: ABAQUS 6.11 and above is required' )
            return False
        
        # Check for numpy
        try: import numpy
        except: 
            showAFXErrorDialog( self.getCurrentDialog(), 'Error: Required module numpy cannot be found' )
            return False   

        # Check for PIL or Pillow
        try: from PIL import Image
        except: 
            showAFXErrorDialog( self.getCurrentDialog(), 'Error: Required module PIL / Pillow cannot be found' )
            return False                
                
        return True

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Register the plug-in

desc = 'An ABAQUS plugin used to generate a time series of virtual x-rays from the output of a bone remodelling analysis.\n\n'                  + \
       'The resulting images can be used to analyse the change in bone density over time in a number of regions of interest or "Gruen Zones", ' + \
       'typically due to changes in loading following insertion of an orthopaedic implant. Associated tool BMDanalyse, available on PyPi, was ' + \
       'created for this sole purpose.\n\nRequires an odb file to be open within the current viewport which has a fieldoutput representing '    + \
       'bone density. Works by mapping the bone density fieldoutput onto a regular 3D grid and then projecting the values onto the three '      + \
       'orthogonal planes. Currently only models with C3D4 or C3D10 elements are supported.'

toolset = getAFXApp().getAFXMainWindow().getPluginToolset()
toolset.registerGuiMenuButton(
    buttonText='pyvXRAY|Create virtual x-rays', 
    object=PyvXRAY_plugin(toolset),
    messageId=AFXMode.ID_ACTIVATE,
    icon=None,
    kernelInitString='import virtualXrays',
    applicableModules=ALL,
    version=__version__,
    author='Michael Hogg',
    description=desc,
    helpUrl='https://code.google.com/p/pyvxray/'
)
