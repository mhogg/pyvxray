# -*- coding: utf-8 -*-

# Copyright (C) 2015 Michael Hogg

# This file is part of pyvXRAY - See LICENSE.txt for information on usage and redistribution

# cython: profile=False
# cython: boundscheck=False
# cython: wraparound=False
# cython: infer_types=True
# cython: nonecheck=False
# cython: cdivision=False

import numpy as np
cimport numpy as np
from libc.math cimport fabs, sqrt
from libc.stdlib cimport malloc, free
ctypedef np.float64_t float64
ctypedef np.int32_t int32


# Note: fmin and fmax are not included in math.h for MSVC! (Although they are included
#       in gcc). Work-around is to write own functions
cdef double fmax(double a, double b):
    if (a>=b): return a
    else:      return b


cdef double fmin(double a, double b):
    if (a<=b): return a
    else:      return b


# Packed structure for use in numpy record array 
cdef packed struct mappedPoint:
    int32 cte
    float64 g,h,r  
       

#cpdef double LinearTetInterpFunc(double[::1] nv, double[::1] ipc):
#    """Shape function for first order tetrahedral (C3D4) element"""    
#    cdef double N[4], U
#    CLinearTetShapeFuncMatrix(ipc,N)
#    U = vectDot(nv,N,4)
#    return U

        
#cpdef double QuadTetInterpFunc(double[::1] nv, double[::1] ipc):
#    """Shape function for second order tetrahedral (C3D10) element"""    
#    cdef double N[10], U
#    CQuadTetShapeFuncMatrix(ipc,N)
#    U = vectDot(nv,N,10)
#    return U    
    
    
cdef int convert3Dto1Dindex(int i, int j, int k, int NX, int NY, int NZ):
    """Converts 3D array index to 1D array index"""
    return i+j*NX+k*NX*NY  
    

cdef int getMinVals(double[:,::1] arr, double minVals[]):
    """Get minimum values in each column of array"""
    cdef int i, numvals, dim
    dim=arr.shape[0]; numvals=arr.shape[1]
    for i in range(dim):
        minVals[i] = arr[i,0]
    for i in range(dim):
        for j in range(1,numvals):
            minVals[i] = fmin(minVals[i],arr[i,j])
    return 0  
    
    
cdef int getMaxVals(double[:,::1] arr, double maxVals[]):
    """Get maximum values in each column of array"""
    cdef int i, numvals, dim
    dim=arr.shape[0]; numvals=arr.shape[1]
    for i in range(dim):
        maxVals[i] = arr[i,0]
    for i in range(dim):
        for j in range(1,numvals):
            maxVals[i] = fmax(maxVals[i],arr[i,j])
    return 0    
    

cdef int getNearest(double[:] arr, double val, int side):
    """Get nearest index to value in list"""    
    cdef int indx, i, numvals = arr.shape[0]
    if val<=arr[0]: return 0
    if val> arr[numvals-1]: return numvals-1
    for i in range(numvals-1):
        if val>arr[i] and val<=arr[i+1]:
            if side==0: return i  
            if side==1: return i+1     

                        
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
 
       
def createElementMap(dict nodeList, np.ndarray[int32,ndim=1] nConnect_labels, 
                     np.ndarray[int32,ndim=2] nConnect_connectivity, int numNodesPerElem,                    
                     double[:] x, double[:] y, double[:] z):
                                        
    """Creates a map between a list of points and a list of solid tetrahedral elements"""
       
    cdef:
        int i,j,k,e,nlabel,NX,NY,NZ,iLow,jLow,kLow,iUpp,jUpp,kUpp,numElems,elemLabel,numGridPoints
        double xLow,yLow,zLow,xUpp,yUpp,zUpp
        double[::1]   gridPointCoords=np.empty(3), ipc=np.zeros(3)  
        double[:,::1] tetNodeCoords=np.empty((3,numNodesPerElem)), JM=np.empty((3,3))
        double[:,::1] dNdG=np.empty((numNodesPerElem,3))
        double tetCoordsLow[3], tetCoordsUpp[3]

    NX=x.shape[0]; NY=y.shape[0]; NZ=z.shape[0]
    numGridPoints = NX*NY*NZ
    dtype=np.dtype([('cte',np.int32),('g',np.float64),('h',np.float64),('r',np.float64)])
    cdef np.ndarray[mappedPoint,ndim=1] elementMap = np.zeros(numGridPoints,dtype=dtype)

    # Select correct function depending on linear or quadratic element
    if numNodesPerElem == 4:  testPointInElement = TestPointInLinearTetElem
    if numNodesPerElem == 10: testPointInElement = TestPointInQuadTetElem
    
    # Create the element map    
    numElems = nConnect_labels.shape[0]    
    for e in range(numElems): 
        
        # Tet element label
        elemLabel = nConnect_labels[e]

        # Tet node coordinates
        for j in range(numNodesPerElem):            
            nlabel = nConnect_connectivity[e,j]
            for i in range(3):                
                tetNodeCoords[i,j] = nodeList[nlabel][i]

        # Get bounding box around tet to limit grid points searched
        getMinVals(tetNodeCoords,tetCoordsLow)
        getMaxVals(tetNodeCoords,tetCoordsUpp)
        iLow = getNearest(x,tetCoordsLow[0],0); iUpp = getNearest(x,tetCoordsUpp[0],1) 
        jLow = getNearest(y,tetCoordsLow[1],0); jUpp = getNearest(y,tetCoordsUpp[1],1)
        kLow = getNearest(z,tetCoordsLow[2],0); kUpp = getNearest(z,tetCoordsUpp[2],1)
                
        # Find intersections between tet and grid points
        for k in range(kLow,kUpp+1):
            for j in range(jLow,jUpp+1):
                for i in range(iLow,iUpp+1):
                    gridPointCoords[0] = x[i] 
                    gridPointCoords[1] = y[j] 
                    gridPointCoords[2] = z[k]
                    foundIntersection  = testPointInElement(gridPointCoords,ipc,tetNodeCoords,dNdG,JM) 
                    if foundIntersection:
                        gridPointIndex = convert3Dto1Dindex(i,j,k,NX,NY,NZ)
                        elementMap[gridPointIndex].cte = elemLabel                        
                        elementMap[gridPointIndex].g   = ipc[0]
                        elementMap[gridPointIndex].h   = ipc[1]
                        elementMap[gridPointIndex].r   = ipc[2]
    
    return elementMap

        
cdef int TestPointInLinearTetElem(double[::1] X2, double[::1] G, double[:,::1] nv, 
                                  double[:,::1] dNdG, double[:,::1] JM):

    """Tests if a point lies within a first order tetrahedral (C3D4) element. 
    This is a direct calculation, performed using the Jacobian matrix"""  

    cdef:
        double tol,lowLim,uppLim,dX[3],N[4],X1[3]
        int result
    
    # Initialise variables
    tol=1.0e-4; lowLim=0.0-tol; uppLim=1.0+tol 
    for i in range(3): G[i]=0.0
          
    CLinearTetShapeFuncMatrix(G,N)
    # nv(3x4) x N(4,1) = X1(3x1)
    result = MatVecMult(nv,N,4,X1,3)  
    for i in range(3):
        dX[i] = X2[i]-X1[i]  
                  
    CLinearTetShapeFuncDerivMatrix(G,dNdG)
    # nv(3x4) x dNdG(4,3) = JM(3x3)    
    result = MatMult(nv,dNdG,JM)

    # Solve system of linear equations, Dx = J Dg
    result = SolveLinearEquations(JM,dX)
    for i in range(3): 
        G[i] += dX[i] 
        
    # Test if point lies within tet element                    
    if((G[0]+G[1]+G[2])<=uppLim          and \
       (G[0]>= lowLim and G[0]<=uppLim)  and \
       (G[1]>= lowLim and G[1]<=uppLim)  and \
       (G[2]>= lowLim and G[2]<=uppLim)):
        return 1
    else:
        return 0                           


cdef int CLinearTetShapeFuncMatrix(double[::1] ipc, double N[]):
       
    cdef double g,h,r

    # Unpack isoparametric coordinates
    g = ipc[0]; h=ipc[1]; r=ipc[2]

    # Element shape functions
    N[0] = (1.0-g-h-r)
    N[1] = g
    N[2] = h
    N[3] = r

    return 0  
    
        
cdef int CLinearTetShapeFuncDerivMatrix(double[::1] ipc, double[:,::1] dNdG):
    
    cdef double g,h,r
    
    # Unpack isoparametric coordinates
    g = ipc[0]; h=ipc[1]; r=ipc[2]
    
    # Partial derivates of shape functions
    # dNdg
    dNdG[0,0] = -1.0
    dNdG[1,0] =  1.0
    dNdG[2,0] =  0.0
    dNdG[3,0] =  0.0           
    # dNdh
    dNdG[0,1] = -1.0
    dNdG[1,1] =  0.0 
    dNdG[2,1] =  1.0
    dNdG[3,1] =  0.0         
    # dNdr    
    dNdG[0,2] = -1.0
    dNdG[1,2] =  0.0
    dNdG[2,2] =  0.0
    dNdG[3,2] =  1.0    
     
    return 0     
   
     
cdef int TestPointInQuadTetElem(double[::1] X2, double[::1] G, double[:,::1] nv, 
                                double[:,::1] dNdG, double[:,::1] JM):
                                                                     
    """Tests if a point lies within a second order tetrahedral (C3D10) element.
     This is an interative process performed using the Newton-Raphson method"""
    
    cdef:
        int maxIter,numIter,result,converged
        double tol,lowLim,uppLim,err1,err2,f[3],X1[3],N[10]

    # Solver parameters
    maxIter=50; tol=1.0e-6; lowLim=0.0-tol; uppLim=1.0+tol
    
    # Set initial values
    numIter=1; converged=0; G[0]=0.0; G[1]=0.0; G[2]=0.0

    # Run iterative loop to find iso-parametric coordinates G=(g,h,r) corresponding to point X2.
    while (numIter<=maxIter): 
        
        CQuadTetShapeFuncMatrix(G,N)
        # nv(3x10) x N(10,1) = X1(3x1)
        result = MatVecMult(nv,N,10,X1,3)
        for i in range(3):
            f[i]= X2[i]-X1[i]
            
        err1=0.0
        for i in range(3):
            err1 += f[i]**2.0
        err1 = sqrt(err1)          
            
        CQuadTetShapeFuncDerivMatrix(G,dNdG)
        # nv(3x10) x dNdG(10,3) = JM(3x3)        
        result = MatMult(nv,dNdG,JM)

        # Solve system of linear equations, -f = J.dX
        result = SolveLinearEquations(JM,f)
        for i in range(3): 
            G[i] += f[i] 

        err2=0.0
        for i in range(3):
            err2 += f[i]**2.0
        err2 = sqrt(err2)

        # Break if error is within tolerance        
        if (err1<=tol and err2<=tol): 
            converged=1
            break 
        
        # Increment loop counter
        numIter+=1                                                                                          

    # Test if point lies within tet element
    if converged and ((G[0]+G[1]+G[2])<=uppLim         and \
                      (G[0]>=lowLim and G[0]<=uppLim)  and \
                      (G[1]>=lowLim and G[1]<=uppLim)  and \
                      (G[2]>=lowLim and G[2]<=uppLim)):
        return 1
    else:
        return 0
  

cdef int CQuadTetShapeFuncMatrix(double[::1] ipc, double N[]):
    
    cdef double g,h,r
    
    # Unpack isoparametric coordinates
    g = ipc[0]; h=ipc[1]; r=ipc[2]
    
    # Element shape functions
    N[0] = (2.0*(1.0-g-h-r)-1.0)*(1.0-g-h-r)
    N[1] = (2.0*g-1.0)*g
    N[2] = (2.0*h-1.0)*h
    N[3] = (2.0*r-1.0)*r
    N[4] = 4.0*(1.0-g-h-r)*g
    N[5] = 4.0*g*h
    N[6] = 4.0*(1.0-g-h-r)*h
    N[7] = 4.0*(1.0-g-h-r)*r
    N[8] = 4.0*g*r 
    N[9] = 4.0*h*r                
    
    return 0   

        
cdef int CQuadTetShapeFuncDerivMatrix(double[::1] ipc, double[:,::1] dNdG):

    cdef double g,h,r
    
    # Unpack isoparametric coordinates
    g = ipc[0]; h=ipc[1]; r=ipc[2]

    # Partial derivates of shape functions
    # dNdg
    dNdG[0,0] =  4.0*(g+h+r-1.0)+1.0
    dNdG[1,0] =  4.0*g-1.0
    dNdG[2,0] =  0.0
    dNdG[3,0] =  0.0
    dNdG[4,0] =  4.0*(1.0-2.0*g-h-r)
    dNdG[5,0] =  4.0*h
    dNdG[6,0] = -4.0*h
    dNdG[7,0] = -4.0*r
    dNdG[8,0] =  4.0*r
    dNdG[9,0] =  0.0                
    # dNdh
    dNdG[0,1] =  4.0*(g+h+r-1.0)+1.0
    dNdG[1,1] =  0.0
    dNdG[2,1] =  4.0*h-1.0
    dNdG[3,1] =  0.0
    dNdG[4,1] = -4.0*g
    dNdG[5,1] =  4.0*g
    dNdG[6,1] =  4.0*(1.0-g-2.0*h-r)
    dNdG[7,1] = -4.0*r
    dNdG[8,1] =  0.0
    dNdG[9,1] =  4.0*r       
    # dNdr    
    dNdG[0,2] =  4.0*(g+h+r-1.0)+1.0
    dNdG[1,2] =  0.0
    dNdG[2,2] =  0.0
    dNdG[3,2] =  4.0*r-1.0
    dNdG[4,2] = -4.0*g
    dNdG[5,2] =  0.0
    dNdG[6,2] = -4.0*h
    dNdG[7,2] =  4.0*(1.0-g-h-2.0*r)
    dNdG[8,2] =  4.0*g
    dNdG[9,2] =  4.0*h
      
    return 0   
    
            
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
