# -*- coding: utf-8 -*-

# Copyright (C) 2015 Michael Hogg

# This file is part of bonemapy - See LICENSE.txt for information on usage and redistribution

import numpy as np
from abaqusConstants import C3D4, C3D4H, C3D10, C3D10H, C3D10I, C3D10M, C3D10MH

# ~~~~~~~~~~ 

class elementC3D4():
    
    def __init__(self):
        self.name     = 'C3D4'
        self.desc     = 'Linear tetrahedral element'
        self.numNodes = 4
        self.N  = np.zeros(self.numNodes)
        self.nv = np.zeros(self.numNodes)
                                            
    def evalN(self,ipc):
        g,h,r = ipc
        self.N[0] = (1.0-g-h-r)
        self.N[1] = g
        self.N[2] = h
        self.N[3] = r
    
    def interp(self,ipc,nv=None):
        self.evalN(ipc)
        if nv==None: return np.dot(self.N,self.nv)
        else:        return np.dot(self.N,nv)
        
    def setNodalValueByIndex(self,indx,val):
        self.nv[indx]=val
    
# ~~~~~~~~~~                              
                                                    
class elementC3D4H(elementC3D4):
    
    def __init__(self):
        elementC3D4.__init__(self)
        self.name = 'C3D4H' 
        self.desc = 'Linear tetrahedral element with hybrid formulation'                              
                            
# ~~~~~~~~~~        
        
class elementC3D10():

    def __init__(self):
        self.name     = 'C3D10'
        self.desc     = 'Quadratic tetrahedral element'
        self.numNodes = 10
        self.N  = np.zeros(self.numNodes)
        self.nv = np.zeros(self.numNodes)   
                              
    def evalN(self,ipc):
        g,h,r = ipc
        self.N[0] = (2.0*(1.0-g-h-r)-1.0)*(1.0-g-h-r)
        self.N[1] = (2.0*g-1.0)*g
        self.N[2] = (2.0*h-1.0)*h
        self.N[3] = (2.0*r-1.0)*r
        self.N[4] =  4.0*(1.0-g-h-r)*g
        self.N[5] =  4.0*g*h
        self.N[6] =  4.0*(1.0-g-h-r)*h
        self.N[7] =  4.0*(1.0-g-h-r)*r
        self.N[8] =  4.0*g*r
        self.N[9] =  4.0*h*r
        
    def interp(self,ipc,nv=None):
        self.evalN(ipc)
        if nv==None: return np.dot(self.N,self.nv)
        else:        return np.dot(self.N,nv)
        
    def setNodalValueByIndex(self,indx,val):
        self.nv[indx]=val
        
# ~~~~~~~~~~                             

class elementC3D10M(elementC3D10):
    
    def __init__(self):
        elementC3D10.__init__(self)
        self.name = 'C3D10M' 
        self.desc = 'Quadratic tetrahedral element with modified formulation'     
                                                  
# ~~~~~~~~~~                              
                                                    
class elementC3D10H(elementC3D10):
    
    def __init__(self):
        elementC3D10.__init__(self)
        self.name = 'C3D10H' 
        self.desc = 'Quadratic tetrahedral element with hybrid formulation'                                                       
                                                                                                                                  
# ~~~~~~~~~~                              
                                                       
class elementC3D10MH(elementC3D10M):
    
    def __init__(self):
        elementC3D10M.__init__(self)
        self.name = 'C3D10MH' 
        self.desc = 'Quadratic tetrahedral element with modified hybrid formulation'                                 

# ~~~~~~~~~~  
                              
class elementC3D10I(elementC3D10):
    
    def __init__(self):
        elementC3D10.__init__(self)
        self.name = 'C3D10I' 
        self.desc = 'Quadratic tetrahedral element with imporved surface stress formulation'
        
# ~~~~~~~~~~         

# Supported element types
seTypes = {}
seTypes['C3D4']    = elementC3D4
seTypes['C3D4H']   = elementC3D4H
seTypes['C3D10']   = elementC3D10
seTypes['C3D10H']  = elementC3D10H
seTypes['C3D10I']  = elementC3D10I
seTypes['C3D10M']  = elementC3D10M
seTypes['C3D10MH'] = elementC3D10MH

