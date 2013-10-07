# -*- coding: utf-8 -*-

# Copyright (C) 2013 Michael Hogg

# This file is part of pyvXRAY - See LICENSE.txt for information on usage and redistribution

# cython: profile=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: infer_types=True
# cython: nonecheck=False
# cython: cdivision=False

import numpy as np
cimport numpy as np
from libc.math cimport fabs, sqrt
from libc.stdlib cimport malloc, free
#ctypedef np.float64_t float64
#ctypedef np.int32_t int32

cdef class elemC3D4:
    
    cdef:
        int numNodes
        double[::1] G
        double[:,::1] nc, JM, dNdG   
        double N[4]

    def __init__(self):
        self.numNodes = 4
        self.G    = np.zeros(3)
        self.nc   = np.empty((3,self.numNodes))
        self.JM   = np.empty((3,3))
        self.dNdG = np.empty((self.numNodes,3))

    cdef int testPointInElem(self,double[::1] X2):
        
        cdef double tol,lowLim,uppLim,X1[3],dX[3]
        cdef int i,result
        tol=1.0e-6; lowLim=0.0-tol; uppLim=1.0+tol 
        
        for i in range(3):
            self.G[i]=0.0
        
        self.evalN()   
        result = MatVecMult(self.nc,self.N,4,X1,3) 
        for i in range(3):
            dX[i] = X2[i]-X1[i]  
                  
        self.evaldNdG()
        result = MatMult(self.nc,self.dNdG,self.JM)

        result = SolveLinearEquations(self.JM,dX)
        for i in range(3): 
            self.G[i] += dX[i] 

        if((self.G[0]+self.G[1]+self.G[2])<=uppLim     and \
           (self.G[0]>= lowLim and self.G[0]<=uppLim)  and \
           (self.G[1]>= lowLim and self.G[1]<=uppLim)  and \
           (self.G[2]>= lowLim and self.G[2]<=uppLim)):            
            return 1
        else:
            return 0
        
    cdef void evalN(self):
        cdef double g,h,r
        g=self.G[0]; h=self.G[1]; r=self.G[2]
        self.N[0] = (1.0-g-h-r)
        self.N[1] = g
        self.N[2] = h
        self.N[3] = r 
    
    cdef void evaldNdG(self): 
        # dNdg
        self.dNdG[0,0] = -1.0
        self.dNdG[1,0] =  1.0
        self.dNdG[2,0] =  0.0
        self.dNdG[3,0] =  0.0           
        # dNdh
        self.dNdG[0,1] = -1.0
        self.dNdG[1,1] =  0.0 
        self.dNdG[2,1] =  1.0
        self.dNdG[3,1] =  0.0         
        # dNdr    
        self.dNdG[0,2] = -1.0
        self.dNdG[1,2] =  0.0
        self.dNdG[2,2] =  0.0
        self.dNdG[3,2] =  1.0
        
    #def getdNdG(self):
    #    return np.array(self.dNdG)
        
    #def setG(self,ipc):
    #    for i in range(3):
    #        self.G[i]=ipc[i]

    cdef void setNodeCoords(self):
        self.nc[:] = np.transpose(np.array([[ 0.0, 0.0, 0.0],
                                            [ 1.0, 0.0, 0.0],
                                            [ 0.0, 1.0, 0.0],
                                            [ 0.0, 0.0, 1.0]]))
                            
    #def getN(self):
    #    self.evalN()  
    #    N = np.zeros(4)
    #    for i in range(4):
    #        N[i]=self.N[i]
    #    return N
                            
    def testPoint(self,pc):
        pc=np.array(pc)
        result = self.testPointInElem(pc)
        if result==True:
            return np.array(self.G)
        else:
            return None
            
cdef int MatMult(double[:,::1] A, double [:,::1] B, double[:,::1] C):
    """Matrix multiplication: A(l,m) x B(m,n) = C(l,n)"""
    cdef int i,j,k,l,m,n
    cdef double csum
    l = A.shape[0]; m=A.shape[1]; n=B.shape[1]
    if m!=B.shape[0] or l!=C.shape[0] or n!=C.shape[1]: return 1
    for i in range(l):
        for j in range(n):
            csum = 0.0
            for k in range(m):
                csum += A[i,k]*B[k,j]
            C[i,j] = csum
    return 0    

                 
cdef int MatVecMult(double[:,::1] M, double v[], int vsize, double y[], int ysize):
    """Multiplies a matrix by a vector, y = M.v"""
    cdef int i,j,m,n
    cdef double csum
    m = M.shape[0]; n=M.shape[1];
    if m!=ysize or n!=vsize: return 1
    for i in range(m):
        csum = 0.0
        for j in range(n):
            csum += M[i,j]*v[j]
        y[i] = csum
    return 0 
    
    
cdef double vectDot(double[::1] v1, double v2[], int v2size):
    """Vector dot product"""
    cdef double vdot
    cdef int i
    vdot = 0.0
    for i in range(v1.shape[0]):
        vdot += v1[i]*v2[i]
    return vdot       

cdef int SolveLinearEquations(double[:,::1] A, double b[]):

    """Solves a system of linear equations, Ax=b, using Gaussian elimination
    with partial pivoting. A must be a square matrix. Both A and b are modified"""

    cdef int i,j,k,m,n,pvtStore
    cdef double pvt,temp
    m = A.shape[1]; n = A.shape[0]
    
    # Allocate memory for dynamic arrays
    cdef double *x  = <double *>malloc(n*sizeof(double))
    cdef int *pivot = <int *>malloc(n*sizeof(int))
    
    # Solve equations
    for j in range(n-1):
        pvt = fabs(A[j,j])
        pvtStore = j
        pivot[j] = pvtStore
        # Find pivot row
        for i in range(j+1,n):
            temp = fabs(A[i,j]) 
            if (temp > pvt):
                pvt = temp
                pvtStore = i
        # Switch rows if necessary
        if (pivot[j] != pvtStore):
            pivot[j] = pvtStore
            pivot[pvtStore] = j   
            for k in range(n):
                temp = A[j,k]
                A[j,k] = A[pivot[j],k]
                A[pivot[j],k] = temp
            temp = b[j]
            b[j] = b[pivot[j]]
            b[pivot[j]] = temp
        # Store multipliers 
        for i in range(j+1,n):
            A[i,j] = A[i,j]/A[j,j]
        # Create zeros below the main diagonal
        for i in range(j+1,n):
            for k in range(j+1,n):
                A[i,k] = A[i,k]-A[i,j]*A[j,k]
            b[i] = b[i]-A[i,j]*b[j]
            
    # Back substitution
    x[n-1] = b[n-1]/A[n-1,n-1]
    for j in range(n-2,-1,-1):
        x[j] = b[j]
        for k in range(n-1,j,-1):
            x[j] = x[j]-x[k]*A[j,k]
        x[j] = x[j]/A[j,j] 
    
    # Return x as b
    for i in range(n):
        b[i] = x[i]
    
    # Free memory
    free(x); free(pivot)
    
    return 0