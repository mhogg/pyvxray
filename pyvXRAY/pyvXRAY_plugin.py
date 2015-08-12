# -*- coding: utf-8 -*-

# Copyright (C) 2015 Michael Hogg

# This file is part of pyvXRAY - See LICENSE.txt for information on usage and redistribution

from abaqusGui import *
from abaqusConstants import ALL, CARTESIAN, SCALAR, INTEGRATION_POINT, CENTROID, ELEMENT_NODAL
from version import version as __version__
# Required to ensure the CSYS list is up to date
from kernelAccess import session  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class PyvXRAY_plugin(AFXForm):

    def __init__(self, owner):
        
        # Construct the base class 
        AFXForm.__init__(self, owner)

        self.odb          = None    
        self.csysList     = None         
        self.stepList     = None
        self.elementSets  = None
        self.scalarList   = None
        self.imageFormats = ['bmp','jpeg','png']
        self.suppLocs     = [INTEGRATION_POINT, CENTROID, ELEMENT_NODAL]
        
        # Keyword definitions
        self.cmd                 = AFXGuiCommand(mode=self, method='createVirtualXrays',objectName='virtualXrays', registerQuery=False)
        self.odbNameKw           = AFXStringKeyword(self.cmd, 'odbName', True, '')
        self.bSetNameKw          = AFXStringKeyword(self.cmd, 'bRegionSetName', True, '')
        self.BMDfonameKw         = AFXStringKeyword(self.cmd, 'BMDfoname', True, '')
        self.showImplantKw       = AFXBoolKeyword(self.cmd,   'showImplant', AFXBoolKeyword.TRUE_FALSE, True, False)
        self.iSetNameKw          = AFXStringKeyword(self.cmd, 'iRegionSetName', True, '')
        self.iDensityKw          = AFXFloatKeyword(self.cmd,  'iDensity', True, 4500)
        self.stepListKw          = AFXStringKeyword(self.cmd, 'stepList', True, '')
        self.csysNameKw          = AFXStringKeyword(self.cmd, 'csysName', True, '')
        self.resGridKw           = AFXFloatKeyword(self.cmd,  'resGrid', True, 2)
        self.imageNameBaseKw     = AFXStringKeyword(self.cmd, 'imageNameBase', True, 'vxray')
        self.preferredXraySizeKw = AFXIntKeyword(self.cmd,    'preferredXraySize', True, 800)
        self.imageFormatKw       = AFXStringKeyword(self.cmd, 'imageFormat', True, self.imageFormats[-1])
        self.smoothKw            = AFXBoolKeyword(self.cmd,   'smooth', AFXBoolKeyword.TRUE_FALSE, True, True)
        self.manualScalingKw     = AFXBoolKeyword(self.cmd,   'manualImageScaling', AFXBoolKeyword.TRUE_FALSE, True, False)
                   
    def getOdbList(self):
        """Get a list of all available odbs in session"""
        self.odbList = session.odbs.keys()
        
    def getFirstOdb(self):
        """Set first odb in first Dialog Box"""
        if self.odbList==None or len(self.odbList)==0: return
        self.setOdb(self.odbList[0])

    def setOdb(self,odbName):
        """Set odb from name"""
        if odbName=='': return
        self.odb = session.odbs[odbName]   

    def getElementSetList(self):
        """Get list of all element sets in the current odb"""
        self.elementSets=[]
        if self.odb==None: return
        # Check part instances
        for instName,inst in self.odb.rootAssembly.instances.items():
            self.elementSets.append('.'.join([instName,'ALL']))
            for setName in inst.elementSets.keys():
                self.elementSets.append('.'.join([instName,setName]))
        # Check assembly            
        for setName in self.odb.rootAssembly.elementSets.keys():
            self.elementSets.append(setName)
        self.elementSets.sort()         
        
    def getCsyses(self): 
        """Get list of all available csyses"""
        # Check scratch odb csyses
        self.csysList = {'Session':[], 'ODB':[]}
        for k,v in session.scratchOdbs.items():
            for csysName,csys in v.rootAssembly.datumCsyses.items():
                if csys.type==CARTESIAN:
                    self.csysList['Session'].append(csysName)
        # Check odb csyses if an odb is open in the current viewport                
        if self.odb != None:
            for csysName,csys in self.odb.rootAssembly.datumCsyses.items():
                if csys.type==CARTESIAN: 
                    self.csysList['ODB'].append(csysName)  

    def getSteps(self):
        """Get list of all available steps"""
        self.stepList=[]
        if self.odb==None: return
        for stepName in self.odb.steps.keys():
            stepNumber = stepName.split('-')[-1]
            self.stepList.append(stepNumber)
        return      
        
    def getScalarList(self):
        """Get list of available scalars. Check last frame in all steps"""
        self.scalarList=[]       
        if self.odb==None: return
        includeList={}; excludeList={}
        for step in self.odb.steps.values():
            frame = step.frames[-1]
            for k in frame.fieldOutputs.keys():
                if excludeList.has_key(k) or includeList.has_key(k): continue
                v   = frame.fieldOutputs[k]               
                loc = [True for loc in v.locations if loc.position in self.suppLocs]                
                if any(loc) and v.type==SCALAR: includeList[k]=1
                else: excludeList[k]=1
        self.scalarList = includeList.keys()
        self.scalarList.sort()

    def getFirstDialog(self):
        """Create the dialog box"""
        # Get odb information to populate the dialog box
        self.getOdbList()
        self.getFirstOdb()
        self.getElementSetList()
        self.getCsyses()
        self.getSteps()
        self.getScalarList()
        # Create dialog box
        import pyvXRAYDB
        return pyvXRAYDB.PyvXRAYDB(self)

    def doCustomChecks(self):
        """Define empty class function doCustomChecks to check user inputs"""
    
        # Check that odb exists
        self.getOdbList()
        if self.odbNameKw.getValue() not in self.odbList:
            showAFXErrorDialog(self.getCurrentDialog(), 'Error: Odb %s does not exist' % self.modelNameKw.getValue())
            return False    
    
        # Check that bone region exists in model
        self.getElementSetList() 
        if self.bSetNameKw.getValue() not in self.elementSets:
            showAFXErrorDialog(self.getCurrentDialog(), 'Error: Bone region %s does not exist' % self.bSetNameKw.getValue())
            return False 
            
        # If implant is requested, check implant inputs
        if self.showImplantKw.getValue():   
            # Check implant region
            self.getElementSetList() 
            if self.iSetNameKw.getValue() not in self.elementSets:
                showAFXErrorDialog(self.getCurrentDialog(), 'Error: Implant region %s does not exist' % self.iSetNameKw.getValue())
                return False 
            # Check input density
            iDensity  = self.iDensityKw.getValue()
            try: float(iDensity)
            except: 
                showAFXErrorDialog(self.getCurrentDialog(), 'Error: Implant density value is not a valid number')
                return False               
            if iDensity<0:
                showAFXErrorDialog(self.getCurrentDialog(), 'Error: Implant density must be greater than 0')
                return False              

        # Check that values in stepList are valid
        stepList = self.stepListKw.getValue()
        try: 
            stepList = [int(s) for s in stepList.replace(',',' ').split()]
        except: 
            showAFXErrorDialog(self.getCurrentDialog(), 'Error: Cannot convert step list values to integers')
            return False
        stepList.sort()
            
        # Check that all steps in step list exist and that density variable exists in all steps (last frame)
        BMDfoname = self.BMDfonameKw.getValue()
        stepInfo  = {}
        for stepName,step in self.odb.steps.items():
            stepNumber = int(stepName.split('-')[-1])
            stepInfo[stepNumber] = step.frames[-1].fieldOutputs.keys()
        for stepNumber in stepList:
            if stepNumber not in stepInfo:
                showAFXErrorDialog(self.getCurrentDialog(), 'Error: Step number %i is not available in odb' % stepNumber)
                return False
            if BMDfoname not in stepInfo[stepNumber]:
                showAFXErrorDialog(self.getCurrentDialog(), 'Error: Density variable %s is not available in Step number %i' % (BMDfoname,stepNumber))
                return False                 
    
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
                
        # Check for Abaqus version >= 6.11 
        majorNumber, minorNumber, updateNumber = getAFXApp().getVersionNumbers()
        if majorNumber==6 and minorNumber < 11:    
            showAFXErrorDialog( self.getCurrentDialog(), 'Error: ABAQUS 6.11 and above is required' )
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
    helpUrl='https://github.com/mhogg/pyvxray'
)
