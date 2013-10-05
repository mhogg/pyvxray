# -*- coding: utf-8 -*-

# Copyright (C) 2013 Michael Hogg

# This file is part of pyvXRAY - See LICENSE.txt for information on usage and redistribution

import os
from abaqus import session, getInputs
from abaqusConstants import ELEMENT_NODAL
from cythonMods import createElementMap, LinearTetInterpFunc, QuadTetInterpFunc

# Use try to prevent error importing missing modules when pyvXRAY plug-in is launched
try:
    import numpy
    from PIL import Image, ImageFilter
except: pass

# ~~~~~~~~~~

def convert3Dto1Dindex(i,j,k,NX,NY,NZ):
    """Converts 3D array index to 1D array index"""
    index = i+j*NX+k*NX*NY
    return index
    
# ~~~~~~~~~~
  
def convert1Dto3Dindex(index,NX,NY,NZ):
    """Converts 1D array index to 1D array index"""
    k = index / (NX*NY)
    j = (index - k*NX*NY) / NX
    i = index - k*NX*NY - j*NX
    return [i,j,k]
    
# ~~~~~~~~~~   

def transformPoint(TM,point):
    """Transforms point using supplied transform"""
    point = numpy.append(point,1.0)
    return numpy.dot(TM,point)[:3]
    
# ~~~~~~~~~~      

def createTransformationMatrix(Ma,Mb,Vab,rel='a'):
    """
    Creates a transformation matrix that can be used to transform a point from csys a to csys b.
    Ma  = 3x3 matrix containing unit vectors of orthogonal coordinate directions for csys a
    Mb  = 3x3 matrix containing unit vectors of orthogonal coordinate directions for csys b
    Vab = 3x1 vector from origin of csys a to csys b
    rel = 'a' or 'b' = Character to indicate if Vab is relative to csys a or csys b
    """
    if rel!='a' and rel!='b': return None
    a1,a2,a3 = Ma
    b1,b2,b3 = Mb
    # Rotation matrix
    R = numpy.identity(4,numpy.float)
    R[0,0:3] = [numpy.dot(b1,a1), numpy.dot(b1,a2), numpy.dot(b1,a3)]
    R[1,0:3] = [numpy.dot(b2,a1), numpy.dot(b2,a2), numpy.dot(b2,a3)]
    R[2,0:3] = [numpy.dot(b3,a1), numpy.dot(b3,a2), numpy.dot(b3,a3)]    
    # Transformation matrix
    if rel=='b':
        Vab = numpy.append(Vab,1.0)
        Vab = numpy.dot(R.T,Vab)[0:3]
    T = numpy.identity(4,numpy.float)     
    T[0:3,3] = -Vab       
    # Transformation matrix
    return numpy.dot(R,T)
    
# ~~~~~~~~~~ 

def getTMfromCsys(odb,csysName):
    if csysName=='GLOBAL': return None
    # Parse coordinate system name
    csysName = csysName.split(r'(')[0].strip()
    # Get ABAQUS datumCsys
    lcsys = None
    # Check odb csyses
    if csysName in odb.rootAssembly.datumCsyses.keys(): 
        lcsys = odb.rootAssembly.datumCsyses[csysName]
    # Check scratch odb csyses
    if odb.path in session.scratchOdbs.keys():
        if csysName in session.scratchOdbs[odb.path].rootAssembly.datumCsyses.keys():
            lcsys = session.scratchOdbs[odb.path].rootAssembly.datumCsyses[csysName]
    if lcsys==None: return None
    # Global coordinate system
    Og = numpy.zeros(3)
    Mg = numpy.identity(3)
    # Local coordinate system
    Ol    = lcsys.origin
    Ml    = numpy.zeros((3,3))
    Ml[0] = lcsys.xAxis/numpy.linalg.norm(lcsys.xAxis) # NOTE: This should already be a unit vector
    Ml[1] = lcsys.yAxis/numpy.linalg.norm(lcsys.yAxis) #       Shouldn't need to normalise
    Ml[2] = lcsys.zAxis/numpy.linalg.norm(lcsys.zAxis)
    # Create transformation matrix
    Vgl = Ol-Og
    TM  = createTransformationMatrix(Mg,Ml,Vgl,rel='a')
    return TM
        
# ~~~~~~~~~~            

def projectXrayPlane(spaceArray3D,whichPlane):
    """Project 3D BMD data onto the specified plane"""
    # Perform the projection by summing along orthogonal axis    
    if whichPlane=='xy':
        projected = numpy.sum(spaceArray3D,axis=2,dtype=spaceArray3D.dtype)
    elif whichPlane=='yz':
        projected = numpy.sum(spaceArray3D,axis=0,dtype=spaceArray3D.dtype)
    elif whichPlane=='xz':
        projected = numpy.sum(spaceArray3D,axis=1,dtype=spaceArray3D.dtype)   
    return projected

# ~~~~~~~~~~  

def writeImageFile(xrayImageFilename,BMDprojected,imageSize,imageFormat='bmp',smooth=True):
    """Create an image from array and write to file"""
    
    # Convert to 8-bit
    xray = BMDprojected.copy()
    xray = numpy.asarray(xray,dtype=numpy.int8)

    # Create image from array.
    xray = xray[:,::-1]       
    xrayImage = Image.fromarray(xray.transpose(),mode='L')
    
    # Resize image
    xsize,ysize = xrayImage.size    
    bigSide = numpy.argmax(xrayImage.size)
    if bigSide==0: scale = float(imageSize)/xsize
    else:          scale = float(imageSize)/ysize
    xsize = int(numpy.rint(scale*xsize))
    ysize = int(numpy.rint(scale*ysize))  
    xrayImage = xrayImage.resize((xsize,ysize),Image.BILINEAR)
    if smooth: xrayImage = xrayImage.filter(ImageFilter.SMOOTH)
    
    # Save xray image to file
    if xrayImageFilename.split('.')[-1]!=imageFormat: xrayImageFilename+='.%s' % imageFormat
    xrayImage.save(xrayImageFilename,imageFormat) 

    return

# ~~~~~~~~~~   

def getPartData(odb,partName,setName,TM,numNodesPerElem):

    """Get part data based on original (undeformed) coordinates"""
    
    p = odb.rootAssembly.instances[partName] 
    pNodes    = p.nodes
    setRegion = p.elementSets[setName]
    pElems    = setRegion.elements
    
    # Create a list of element connectivities (list of nodes connected to each element)    
    setNodeLabs     = {}
    numElems        = len(pElems)
    elemConnect     = numpy.zeros(numElems,dtype=[('label','|i4'),('connectivity','|i4',(numNodesPerElem,))])
    for e in xrange(numElems):
        elem = pElems[e]
        conn = elem.connectivity
        elemConnect[e] = (elem.label, conn)
        for n in conn:        
            setNodeLabs[n] = 1
    
    # Create a dictionary of node labels and node coordinates for the entire part instance
    numNodes    = len(pNodes)
    numSetNodes = len(setNodeLabs) 
    nodeCount   = 0
    setNodes    = numpy.zeros(numSetNodes,dtype=[('label','|i4'),('coordinates','|f4',(3,))])
    for n in xrange(numNodes):
        node  = pNodes[n]
        label = node.label
        if label in setNodeLabs:
            setNodes[nodeCount] = (label, node.coordinates) 
            nodeCount += 1
        
    # Transform the coordinates from the global csys to the local csys
    if TM is not None:
        for i in xrange(numSetNodes):
            setNodes['coordinates'][i] = transformPoint(TM,setNodes['coordinates'][i])
        
    # Get bounding box
    low  = numpy.min(setNodes['coordinates'],axis=0)
    upp  = numpy.max(setNodes['coordinates'],axis=0) 
    bbox = (low,upp)

    # Convert setNodes to a dictionary for fast indexing by node label
    setNodeList = dict(zip(setNodes['label'],setNodes['coordinates']))       
    
    return setRegion,setNodeList,elemConnect,bbox

# ~~~~~~~~~~ 

def getElementInfo(odb,partName,setName):
    """Get element type and number of nodes per element"""
    
    elements = odb.rootAssembly.instances[partName].elementSets[setName].elements 
    eTypes   = {}
    for e in elements:
        eTypes[e.type]=1
    eTypes=[str(eType) for eType in eTypes.keys()]  
    
    return checkElementTypes(eTypes,partName)
    
# ~~~~~~~~~~ 

def checkElementTypes(eTypes,partName):
    """Check element type. Only supported elements and a single element type per part are allowed"""
    
    # Perform checks of element types    
    supportedElements={}
    supportedElements['C3D4']  =  4
    supportedElements['C3D10'] = 10     
    
    # Check 1: Check for unsupported element types   
    usTypes=[]
    for eType in eTypes:
        if not any([True for seType in supportedElements.keys() if seType in eType]):
            usTypes.append(eType)
    if len(usTypes)>0:
        print 'Element types %s in part %s are not supported' % (', '.join(usTypes),partName)   
        return None

    # Check 2: Check that the part consists of single element type only
    sTypes={}
    for eType in eTypes:
        for sType in [seType for seType in supportedElements.keys() if seType in eType]:
            sTypes[sType]=1
    sTypes = sTypes.keys()
    if len(sTypes)==1:
        partElemType = sTypes[0] 
    else:
        print 'Multiple element types in part %s not supported' % partName 
        return None
           
    return partElemType,supportedElements[partElemType] 
    
# ~~~~~~~~~~      

def createVirtualXrays(odbName,bRegionSetName,BMDfoname,showImplant,iRegionSetName,
                       iDensity,stepList,csysName,resGrid,imageNameBase,preferredXraySize,
                       imageFormat,smooth=False,manualImageScaling=False):
    """Creates virtual x-rays from an ABAQUS odb file. The odb file should contain \n""" + \
    """a number of steps with a fieldoutput variable representing bone mineral density (BMD)"""
        
    # User message
    print '\npyvXRAY: Create virtual x-rays plugin'
    
    # Process inputs    
    resGrid           = float(resGrid)
    stepList          = [int(s) for s in stepList.replace(',',' ').split()]
    preferredXraySize = int(preferredXraySize)
        
    # Set variables
    dx,dy,dz  = (resGrid,)*3
    iDensity /= 1000.    
    odb       = session.odbs[odbName]

    # Get element information for each part. An error will be raised if 
    # unsupported elements are detected or a part is made up of more than one 
    # element type. Currently the only supported elements are C3D4 and C3D10 (and
    # variants of these such as C3D10M etc). The bone and implant parts may use
    # different element types from each other. 
    bElementInfo = getElementInfo(odb,bPartName,bSetName)
    if bElementInfo==None: return
    bElemType,bNumNodesPerElem = bElementInfo
    if bElemType=='C3D4' : tetInterpFunc = LinearTetInterpFunc
    if bElemType=='C3D10': tetInterpFunc = QuadTetInterpFunc
    if showImplant: 
        iElementInfo = getElementInfo(odb,iPartName,iSetName) 
        if iElementInfo==None: return
        iElemType,iNumNodesPerElem = iElementInfo
    
    # Get transformation matrix to convert from global to local coordinate system
    TM = getTMfromCsys(odb,csysName)
    print 'X-ray views will be relative to %s' % csysName

    # Get part data and create a bounding box. The bounding box should include the implant if specified
    bRegion,bNodeList,bElemConnect,bBBox = getPartData(odb,bPartName,bSetName,TM,bNumNodesPerElem)
    if showImplant:    
        iRegion,iNodeList,iElemConnect,iBBox = getPartData(odb,iPartName,iSetName,TM,iNumNodesPerElem)
        bbLow = numpy.min((bBBox[0],iBBox[0]),axis=0)
        bbUpp = numpy.max((bBBox[1],iBBox[1]),axis=0)
    else:
        bbLow,bbUpp = bBBox 
       
    border   = 0.05*(bbUpp-bbLow)  
    bbLow    = bbLow - border
    bbUpp    = bbUpp + border 
    bbSides  = bbUpp - bbLow
    x0,y0,z0 = bbLow 
    xN,yN,zN = bbUpp 
    lx,ly,lz = bbSides

    # Generate Xray grid
    NX = int(numpy.ceil(lx/dx+1))
    x  = numpy.linspace(x0,xN,NX)
    NY = int(numpy.ceil(ly/dy+1))
    y  = numpy.linspace(y0,yN,NY)
    NZ = int(numpy.ceil(lz/dz+1))
    z  = numpy.linspace(z0,zN,NZ)  
    
    # Create element map for the implant, map tp 3D space array and then project onto 3 planes 
    if showImplant: 
        # Get element map       
        iElementMap  = createElementMap(iNodeList,iElemConnect['label'],iElemConnect['connectivity'],iNumNodesPerElem,x,y,z)        
        # Mask 3D array
        iMask = numpy.zeros((NX,NY,NZ),dtype=numpy.float64)   
        for gpi in xrange(iElementMap.size):
            gridPoint = iElementMap[gpi]
            if gridPoint['cte'] > 0:
                i,j,k = convert1Dto3Dindex(gpi,NX,NY,NZ)
                iMask[i,j,k] = iDensity
        # Create projections of 3D space array onto planes 
        iProjectedXY = projectXrayPlane(iMask,'xy')
        iProjectedYZ = projectXrayPlane(iMask,'yz')
        iProjectedXZ = projectXrayPlane(iMask,'xz')
        # Create xrays of implant without bone        
        iprojXY = iProjectedXY.copy()
        iprojYZ = iProjectedYZ.copy()
        iprojXZ = iProjectedXZ.copy()        
        prXY    = [numpy.min(iprojXY),numpy.max(iprojXY)]
        prYZ    = [numpy.min(iprojYZ),numpy.max(iprojYZ)]
        prXZ    = [numpy.min(iprojXZ),numpy.max(iprojXZ)]        
        iprojXY[:,:] = (iprojXY[:,:]-prXY[0])/(prXY[1]-prXY[0])*255.
        iprojYZ[:,:] = (iprojYZ[:,:]-prYZ[0])/(prYZ[1]-prYZ[0])*255. 
        iprojXZ[:,:] = (iprojXZ[:,:]-prXZ[0])/(prXZ[1]-prXZ[0])*255.         
        writeImageFile('implant_XY',iprojXY,preferredXraySize,imageFormat,smooth)
        writeImageFile('implant_YZ',iprojYZ,preferredXraySize,imageFormat,smooth)
        writeImageFile('implant_XZ',iprojXZ,preferredXraySize,imageFormat,smooth)

    # Create the element map for the bone
    bElementMap = createElementMap(bNodeList,bElemConnect['label'],bElemConnect['connectivity'],bNumNodesPerElem,x,y,z)
    
    # Interpolate HU values from tet mesh onto grid using quadratic tet shape function
    # (a) Get HU values from frame
    numSteps  = len(stepList)
    xraysXY   = numpy.zeros((numSteps,NX,NY),dtype=numpy.float64)
    xraysYZ   = numpy.zeros((numSteps,NY,NZ),dtype=numpy.float64)
    xraysXZ   = numpy.zeros((numSteps,NX,NZ),dtype=numpy.float64)
    mappedBMD = numpy.zeros((NX,NY,NZ),dtype=numpy.float64)
    
    for s in xrange(numSteps):
        # Step details
        stepId   = stepList[s]               
        stepName = "Step-%i" % (stepId)
        frame    = odb.steps[stepName].frames[-1]
        # Get BMD data for bRegion in current frame
        BMDfov    = frame.fieldOutputs[BMDfoname].getSubset(region=bRegion, position=ELEMENT_NODAL).values
        BMDvalues = {}
        for i in xrange(len(BMDfov)):
            val = BMDfov[i]            
            elementLabel = val.elementLabel            
            if elementLabel not in BMDvalues:
                BMDvalues[elementLabel] = numpy.zeros(bNumNodesPerElem)
                count=0
            else:
                count+=1
            BMDvalues[elementLabel][count] = val.data        
        # Perform the interpolation from elementMap to 3D space array
        for gpi in xrange(bElementMap.size):
            gridPoint = bElementMap[gpi]
            cte = gridPoint['cte']
            if cte > 0:
                nv  = BMDvalues[cte] 
                ipc = numpy.array([gridPoint['g'],gridPoint['h'],gridPoint['r']])
                i,j,k = convert1Dto3Dindex(gpi,NX,NY,NZ)
                mappedBMD[i,j,k] = tetInterpFunc(nv,ipc)
        # Project onto orthogonal planes    
        xraysXY[s] = projectXrayPlane(mappedBMD,'xy')
        xraysYZ[s] = projectXrayPlane(mappedBMD,'yz')
        xraysXZ[s] = projectXrayPlane(mappedBMD,'xz')
        
    # Get min/max pixel values. Use zero for lower limit (corresponding to background)
    prXY = [0.,numpy.max(xraysXY)]
    prYZ = [0.,numpy.max(xraysYZ)]
    prXZ = [0.,numpy.max(xraysXZ)]
    # Allow user to change scale factors if desired
    if manualImageScaling:
        fields = (('X-Y','%.6f'%prXY[1]),('Y-Z','%.6f'%prYZ[1]),('X-Z','%.6f'%prXZ[1]))
        usf = getInputs(fields=fields,label='X-ray image scale factors:')
        if usf[0] != None:
            try:    usf = [float(sf) for sf in usf]
            except: print 'Error in user supplied X-ray image scale factors. Using pyvXRAY values'
            else:   prXY[1],prYZ[1],prXZ[1] = usf      
    # Add projected implant to projected bone
    if showImplant:
        xraysXY[:] += iProjectedXY
        xraysYZ[:] += iProjectedYZ
        xraysXZ[:] += iProjectedXZ
    # Scale each projection using pixel range from bone
    xraysXY[:,:,:] = (xraysXY[:,:,:]-prXY[0])/(prXY[1]-prXY[0])*255.      
    xraysYZ[:,:,:] = (xraysYZ[:,:,:]-prYZ[0])/(prYZ[1]-prYZ[0])*255.           
    xraysXZ[:,:,:] = (xraysXZ[:,:,:]-prXZ[0])/(prXZ[1]-prXZ[0])*255.                
    xraysXY[numpy.where(xraysXY<0.)]   = 0.
    xraysYZ[numpy.where(xraysYZ<0.)]   = 0.
    xraysXZ[numpy.where(xraysXZ<0.)]   = 0.
    xraysXY[numpy.where(xraysXY>255.)] = 255.
    xraysYZ[numpy.where(xraysYZ>255.)] = 255.
    xraysXZ[numpy.where(xraysXZ>255.)] = 255.
    
    # Create images
    for s in xrange(numSteps):
        stepId = stepList[s]
        writeImageFile(('%s_XY_%i' % (imageNameBase,stepId)),xraysXY[s,:,:],preferredXraySize,imageFormat,smooth)
        writeImageFile(('%s_YZ_%i' % (imageNameBase,stepId)),xraysYZ[s,:,:],preferredXraySize,imageFormat,smooth) 
        writeImageFile(('%s_XZ_%i' % (imageNameBase,stepId)),xraysXZ[s,:,:],preferredXraySize,imageFormat,smooth)
        
    # User message
    print 'Virtual x-rays have been created in %s\n' % os.getcwd()
    