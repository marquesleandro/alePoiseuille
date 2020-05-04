# ======================================
# Code created by Leandro Marques
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ======================================

# This code is used to assemble the elementary arrays


# ------------------------------------------------------------------
# Use:
# linear = trielem.Linear(mesh.x, mesh.y, mesh.IEN)
# 
# for e in range(0,mesh.nelem):
#  linear.numerical(e)
#
#  for i in range(0,3):
#   ii = mesh.IEN[e][i]
#
#   for j in range(0,3):
#    jj = mesh.IEN[e][j]
#
#    K[ii][jj] += linear.kxx[i][j] + linear.kyy[i][j]
# ------------------------------------------------------------------


# Reference Analytic: Fundamentals of the Finite
#                     Element Method for Heat Transfer
#                     and Fluid Flow - Lewis, Nithiarasu,
#                     Seetharamu - pg. 196-200
#                     For Q_elem pg. 126
#                     For 1D pg. 193


import sys
import numpy as np

class Element1D:
 def __init__(_self, _x, _IEN, _GAUSSPOINTS):
  _self.x = _x
  _self.IEN = _IEN


  if _GAUSSPOINTS == 3:
   _self.NUMGAUSS = 3  #Number of Gauss Points


   #                                l1     
   _self.GQPoints = np.array([[-0.774596669], 
                              [ 0.000000000], 
                              [ 0.774596669]])


   #                                 w
   _self.GQWeights = np.array([[0.555555556], 
                               [0.888888889], 
                               [0.555555556]])



 
  elif _GAUSSPOINTS == 4:
   _self.NUMGAUSS = 4  #Number of Gauss Points


   #                               l1     
   _self.GQPoints = np.array([[-0.861136], 
                              [-0.339981], 
                              [ 0.339981], 
                              [ 0.861136]])


   #                                w
   _self.GQWeights = np.array([[0.347855], 
                               [0.652145], 
                               [0.652145], 
                               [0.347855]])




  elif _GAUSSPOINTS == 5:
   _self.NUMGAUSS = 5  #Number of Gauss Points


   #                               l1     
   _self.GQPoints = np.array([[-0.906180], 
                              [-0.538469], 
                              [ 0.000000], 
                              [ 0.538469], 
                              [ 0.906180]])


   #                                w
   _self.GQWeights = np.array([[0.236927], 
                               [0.478629], 
                               [0.568889], 
                               [0.478629], 
                               [0.236927]])


  elif _GAUSSPOINTS == 10:
   _self.NUMGAUSS = 10  #Number of Gauss Points


   #                               l1     
   _self.GQPoints = np.array([[-0.148874], 
                              [ 0.148874], 
                              [-0.433953], 
                              [ 0.433953], 
                              [-0.679409], 
                              [ 0.679409], 
                              [-0.865063], 
                              [ 0.865063], 
                              [-0.973906], 
                              [ 0.973906]])


   #                                w     
   _self.GQWeights = np.array([[0.295524], 
                               [0.295524], 
                               [0.269266], 
                               [0.269266], 
                               [0.219086], 
                               [0.219086], 
                               [0.149451], 
                               [0.149451], 
                               [0.066671], 
                               [0.066671]])





  else:
   print ""
   print " Error: Gauss Points not found"
   print ""
   sys.exit()




 def linear(_self,_e):
  _self.NUMNODE = 2  #Linear One-dimensional Element - 2 Nodes

  N = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  dNdl1 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  
  dxdl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)

  J = np.zeros([1,1], dtype = float)
  jacobian = np.zeros([_self.NUMGAUSS,1], dtype = float)
  
  dNdx = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)


  # A loop is required for each pair of coordinates (l1,l2)
  for k in range(0,_self.NUMGAUSS):
    
   # Area Coordinates
   L1 = (1.0/2.0)*(1.0 - _self.GQPoints[k][0])      #L1 = (1/2)*(1 - l1)
   L2 = (1.0/2.0)*(1.0 + _self.GQPoints[k][0])      #L2 = (1/2)*(1 + l1)

   # Shape Functions
   N[k][0] = L1  #N1 = L1
   N[k][1] = L2  #N2 = L2

   # Shape Functions Derivatives in respect to l1
   dNdl1[k][0] = -0.5   #dN1/dl1
   dNdl1[k][1] =  0.5   #dN2/dl1


   # Coordinate Transfomation
   # Lewis pag. 64 Eq. 3.108 for 1D
   for i in range(0,_self.NUMNODE):
    ii = _self.IEN[_e][i]
    
    dxdl1[k] += _self.x[ii]*dNdl1[k][i]   # dx/dl1 = x1*dN1/dl1 + x2*dN2/dl1 + x3*dN3/dl1 ...

   # Jacobian Matrix
   # Lewis pag. 64 Eq. 3.108
   J[0][0] = dxdl1[k]

   jacobian[k] = np.linalg.det(J)


   # Shape Functions Derivatives in respect to x
   # Lewis pag. 63 Eq. 3.107
   for i in range(0,_self.NUMNODE):
    dNdx[k][i] = (1.0/jacobian[k])*dNdl1[k][i]

 
  _self.mass = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  
  # Elementary Matrices 
  for k in range(0,_self.NUMGAUSS): 
   for i in range(0,_self.NUMNODE):
    for j in range(0,_self.NUMNODE):
    
     _self.mass[i][j] += N[k][i]*N[k][j]*jacobian[k]*_self.GQWeights[k]
     _self.kx[i][j] += dNdx[k][i]*dNdx[k][j]*jacobian[k]*_self.GQWeights[k]
     _self.gx[i][j] += dNdx[k][i]*N[k][j]*jacobian[k]*_self.GQWeights[k]
     _self.dx[i][j] += dNdx[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]
   

 def quadratic(_self,_e):
  _self.NUMNODE = 3  #Quadratic One-dimensional Element - 3 Nodes

  N = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  dNdl1 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  
  dxdl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)

  J = np.zeros([1,1], dtype = float)
  jacobian = np.zeros([_self.NUMGAUSS,1], dtype = float)
  
  dNdx = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)


  # A loop is required for each pair of coordinates (l1,l2)
  for k in range(0,_self.NUMGAUSS):
    
   # Area Coordinates
   L1 = (1.0/2.0)*(1.0 - _self.GQPoints[k][0])      #L1 = (1/2)*(1 - l1)
   L2 = (1.0/2.0)*(1.0 + _self.GQPoints[k][0])      #L2 = (1/2)*(1 + l1)

   # Shape Functions
   # Lewis pag. 63 Eq. 3.104
   # N3 is quadratic node
   N[k][0] = L1*(-_self.GQPoints[k][0])  #N1 = (1/2)*(1 - l1)*(-l1)
   N[k][1] = L2*( _self.GQPoints[k][0])  #N2 = (1/2)*(1 + l1)*l1
   N[k][2] = 4.0*L1*L2                   #N3 = (1 - l1**2)

   # Shape Functions Derivatives in respect to l1
   #dN3/dl1 is quadratic node
   dNdl1[k][0] = (1.0/2.0)*(-1.0 + 2.0*_self.GQPoints[k][0])   #dN1/dl1
   dNdl1[k][1] = (1.0/2.0)*( 1.0 + 2.0*_self.GQPoints[k][0])   #dN2/dl1
   dNdl1[k][2] = -2.0*_self.GQPoints[k][0]                     #dN3/dl1


   # Coordinate Transfomation
   # Lewis pag. 64 Eq. 3.108 for 1D
   # IEN[e][2] is quadratic node, so N3 and dN3/dl1 is quadratic 
   for i in range(0,_self.NUMNODE):
    ii = _self.IEN[_e][i]
    
    dxdl1[k] += _self.x[ii]*dNdl1[k][i]   # dx/dl1 = x1*dN1/dl1 + x2*dN2/dl1 + x3*dN3/dl1 ...

   # Jacobian Matrix
   # Lewis pag. 64 Eq. 3.108
   J[0][0] = dxdl1[k]

   jacobian[k] = np.linalg.det(J)


   # Shape Functions Derivatives in respect to x
   # Lewis pag. 63 Eq. 3.107
   for i in range(0,_self.NUMNODE):
    dNdx[k][i] = (1.0/jacobian[k])*dNdl1[k][i]

 
  _self.mass = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  
  # Elementary Matrices 
  for k in range(0,_self.NUMGAUSS): 
   for i in range(0,_self.NUMNODE):
    for j in range(0,_self.NUMNODE):
    
     _self.mass[i][j] += N[k][i]*N[k][j]*jacobian[k]*_self.GQWeights[k]
     _self.kx[i][j] += dNdx[k][i]*dNdx[k][j]*jacobian[k]*_self.GQWeights[k]
     _self.gx[i][j] += dNdx[k][i]*N[k][j]*jacobian[k]*_self.GQWeights[k]
     _self.dx[i][j] += dNdx[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]
   


 # Reference Analytic: Fundamentals of the Finite
 #                     Element Method for Heat Transfer
 #                     and Fluid Flow - Lewis, Nithiarasu,
 #                     Seetharamu - pg. 196-200
 #                     For Q_elem pg. 126
 #                     For 1D pg. 193
 def linear_analytic(_self, _e):
  v1 = _self.IEN[_e][0]
  v2 = _self.IEN[_e][1]
  
  dx = _self.x[v2] - _self.x[v1]

  _self.kx = (1.0/dx)*np.array([[1,-1],[-1,1]])
  _self.mass = (dx/6.0)*np.array([[2,1],[1,2]])
  _self.gx = (1.0/2.0)*np.array([[-1,1],[-1,1]])



 # http://kis.tu.kielce.pl//mo/COLORADO_FEM/colorado/IFEM.Ch32.pdf
 # pag. 32-11
 # Eq. 32.24
 def quad_analytic(_self, _e):
  v1 = _self.IEN[_e][0]
  v2 = _self.IEN[_e][1]
  v3 = _self.IEN[_e][2]
  
  dx = _self.x[v2] - _self.x[v1]

  _self.kx = (1.0/(6.0*dx))*np.array([[14,2,-16],[2,14,-16],[-16,-16,32]])
  _self.mass = (dx/90.0)*np.array([[12,-3,6],[-3,12,6],[6,6,48]])
  # _self.gx = (1.0/2.0)*np.array([[-1,1],[-1,1]]) # verify 



class Element2D:
 def __init__(_self, _x, _y, _IEN, _GAUSSPOINTS):
  _self.x = _x
  _self.y = _y
  _self.IEN = _IEN


  if _GAUSSPOINTS == 3:
   _self.NUMGAUSS = 3  #Number of Gauss Points


   #                                 l1                 l2
   _self.GQPoints = np.array([[0.16666666666667, 0.16666666666667], 
                              [0.16666666666667, 0.66666666666667], 
                              [0.66666666666667, 0.16666666666667]])


   #                                    w
   _self.GQWeights = np.array([[0.333333333333333], 
                               [0.333333333333333], 
                               [0.333333333333333]])





  elif _GAUSSPOINTS == 4:
   _self.NUMGAUSS = 4  #Number of Gauss Points


   #                                 l1                 l2
   _self.GQPoints = np.array([[0.33333333333333, 0.33333333333333], 
                              [0.20000000000000, 0.20000000000000], 
                              [0.20000000000000, 0.60000000000000], 
                              [0.60000000000000, 0.20000000000000]])


   #                                    w
   _self.GQWeights = np.array([[-0.56250000000000], 
                               [0.520833333333333], 
                               [0.520833333333333], 
                               [0.520833333333333]])



  
  elif _GAUSSPOINTS == 6:
   _self.NUMGAUSS = 6  #Number of Gauss Points


   #                                 l1                 l2
   _self.GQPoints = np.array([[0.44594849091597, 0.44594849091597],
                              [0.44594849091597, 0.10810301816807],
                              [0.10810301816807, 0.44594849091597],
                              [0.09157621350977, 0.09157621350977],
                              [0.09157621350977, 0.81684757298046],
                              [0.81684757298046, 0.09157621350977]])


   #                                    w
   _self.GQWeights = np.array([[0.22338158967801],
                               [0.22338158967801],
                               [0.22338158967801],
                               [0.10995174365532],
                               [0.10995174365532],
                               [0.10995174365532]])


  elif _GAUSSPOINTS == 12:
   _self.NUMGAUSS = 12  #Number of Gauss Points

   #                                 l1                 l2
   _self.GQPoints = np.array([[0.24928674517091, 0.24928674517091],
                              [0.24928674517091, 0.50142650965818],
                              [0.50142650965818, 0.24928674517091],
                              [0.06308901449150, 0.06308901449150],
                              [0.06308901449150, 0.87382197101700],
                              [0.87382197101700, 0.06308901449150],
                              [0.31035245103378, 0.63650249912140],
                              [0.63650249912140, 0.05314504984482],
                              [0.05314504984482, 0.31035245103378],
                              [0.63650249912140, 0.31035245103378],
                              [0.31035245103378, 0.05314504984482],
                              [0.05314504984482, 0.63650249912140]])

   #                                    w
   _self.GQWeights = np.array([[0.11678627572638],
                               [0.11678627572638],
                               [0.11678627572638],
                               [0.05084490637021],
                               [0.05084490637021],
                               [0.05084490637021],
                               [0.08285107561837],
                               [0.08285107561837],
                               [0.08285107561837],
                               [0.08285107561837],
                               [0.08285107561837],
                               [0.08285107561837]])
   '''
 
 
    #                                 l1                 l2
   _self.GQPoints = np.array([[0.219429982550000, 0.561140034900000],
                              [0.561140034900000, 0.219429982550000],
                              [0.219429982550000, 0.219429982550000],
                              [0.480137964112000, 0.039724071775600],
                              [0.039724071775600, 0.480137964112000],
                              [0.480137964112000, 0.480137964112000],
                              [0.019371724361200, 0.839009259715000],
                              [0.141619015924000, 0.019371724361200],
                              [0.839009259715000, 0.141619015924000],
                              [0.141619015924000, 0.839009259715000],
                              [0.019371724361200, 0.141619015924000],
                              [0.839009259715000, 0.019371724361200]])


    #                                    w
   _self.GQWeights = np.array([[0.171333124153000],
                               [0.171333124153000],
                               [0.171333124153000],
                               [0.080731089593000],
                               [0.080731089593000],
                               [0.080731089593000],
                               [0.040634559793700],
                               [0.040634559793700],
                               [0.040634559793700],
                               [0.040634559793700],
                               [0.040634559793700],
                               [0.040634559793700]])
 
 
  '''
  else:
   print ""
   print " Error: Gauss Points not found"
   print ""
   sys.exit()


 def linear(_self,_e):
  _self.NUMNODE = 3  #Linear Triangular Element - 3 Nodes
  
  N = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  dNdl1 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  dNdl2 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  
  dxdl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dxdl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dydl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dydl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)

  dNdx = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  dNdy = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  J = np.zeros([2,2], dtype = float)
  jacobian = np.zeros([_self.NUMGAUSS,1], dtype = float)


  # A loop is required for each pair of coordinates (l1,l2)
  for k in range(0,_self.NUMGAUSS):
    
   # Area Coordinates
   # Lewis pag. 67 Eq. 3.129
   L1 = 1.0 - _self.GQPoints[k][0] - _self.GQPoints[k][1]   #L1 = 1 - l1 - l2
   L2 = _self.GQPoints[k][0]                                #L2 = l1
   L3 = _self.GQPoints[k][1]                                #L3 = l2

   # Shape Functions
   # Lewis pag. 67 Eq. 3.129
   N[k][0] = L1  #N1 = L1
   N[k][1] = L2  #N2 = L2
   N[k][2] = L3  #N3 = L3

   # Shape Functions Derivatives in respect to l1
   dNdl1[k][0] = -1.0   #dN1/dl1
   dNdl1[k][1] =  1.0   #dN2/dl1
   dNdl1[k][2] =  0.0   #dN3/dl1

   # Shape Functions Derivatives in respect to l2
   dNdl2[k][0] = -1.0   #dN1/dl2
   dNdl2[k][1] =  0.0   #dN2/dl2
   dNdl2[k][2] =  1.0   #dN3/dl2

   # Coordinate Transfomation
   # Lewis pag. 64 Eq. 3.108 for 1D
   # Lewis pag. 66 Eq. 3.121 for 2D
   for i in range(0,_self.NUMNODE):
    ii = _self.IEN[_e][i]
    
    dxdl1[k] += _self.x[ii]*dNdl1[k][i]   # dx/dl1 = x1*dN1/dl1 + x2*dN2/dl1 + x3*dN3/dl1 ...
    dydl1[k] += _self.y[ii]*dNdl1[k][i]   # dy/dl1 = y1*dN1/dl1 + y2*dN2/dl1 + y3*dN3/dl1 ...
    dxdl2[k] += _self.x[ii]*dNdl2[k][i]   # dx/dl2 = x1*dN1/dl2 + x2*dN2/dl2 + x3*dN3/dl2 ...
    dydl2[k] += _self.y[ii]*dNdl2[k][i]   # dy/dl2 = y1*dN1/dl2 + y2*dN2/dl2 + y3*dN3/dl2 ...


   # Jacobian Matrix
   # Lewis pag. 66 Eq. 3.121
   J[0][0] = dxdl1[k]
   J[0][1] = dydl1[k]
   J[1][0] = dxdl2[k]
   J[1][1] = dydl2[k]

   jacobian[k] = np.linalg.det(J)


   # Shape Functions Derivatives in respect to x and y
   # Lewis pag. 65 Eq. 3.116
   for i in range(0,_self.NUMNODE):
    dNdx[k][i] = (1.0/jacobian[k])*( dNdl1[k][i]*dydl2[k] - dNdl2[k][i]*dydl1[k])
    dNdy[k][i] = (1.0/jacobian[k])*(-dNdl1[k][i]*dxdl2[k] + dNdl2[k][i]*dxdl1[k])



  _self.mass = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kxx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kxy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kyx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kyy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  

  # Elementary Matrices
  # P.S: It is divided by 2 due to relation 
  # of parallelogram and triangle areas --> DxDy = (jacobian*Dl1*Dl2)/2
  for k in range(0,_self.NUMGAUSS): 
   for i in range(0,_self.NUMNODE):
    for j in range(0,_self.NUMNODE):
    
     _self.mass[i][j] += N[k][i]*N[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.kxx[i][j] += dNdx[k][i]*dNdx[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kxy[i][j] += dNdx[k][i]*dNdy[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kyx[i][j] += dNdy[k][i]*dNdx[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kyy[i][j] += dNdy[k][i]*dNdy[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.gx[i][j] += dNdx[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.gy[i][j] += dNdy[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.dx[i][j] += dNdx[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.dy[i][j] += dNdy[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
  

 def mini(_self,_e):
  _self.NUMNODE = 4  #Mini Triangular Element - 4 Nodes
  
  N = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  dNdl1 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  dNdl2 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  
  dxdl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dxdl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dydl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dydl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)

  dNdx = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  dNdy = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  J = np.zeros([2,2], dtype = float)
  jacobian = np.zeros([_self.NUMGAUSS,1], dtype = float)



  # A loop is required for each pair of coordinates (l1,l2)
  for k in range(0,_self.NUMGAUSS):
    
   # Area Coordinates
   L1 = _self.GQPoints[k][0]                                #L1 = l1
   L2 = _self.GQPoints[k][1]                                #L2 = l2
   L3 = 1.0 - _self.GQPoints[k][0] - _self.GQPoints[k][1]   #L3 = 1 - l1 - l2
   
   # Shape Functions
   # Mini2D matlab Gustavo
   N[k][0] = L1 - 9.0*L1*L2*L3    #N1 = L1-9*L1*L2*L3
   N[k][1] = L2 - 9.0*L1*L2*L3    #N2 = L2-9*L1*L2*L3
   N[k][2] = L3 - 9.0*L1*L2*L3    #N3 = L3-9*L1*L2*L3
   N[k][3] = 27.0*L1*L2*L3        #N4 = 27*L1*L2*L3

   # Shape Functions Derivatives in respect to l1
   #dN1/dl1
   dNdl1[k][0] =  18.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] + 9.0*_self.GQPoints[k][1]**2\
                 - 9.0*_self.GQPoints[k][1] + 1.0 

   #dN2/dl1
   dNdl1[k][1] =  18.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] + 9.0*_self.GQPoints[k][1]**2\
                 - 9.0*_self.GQPoints[k][1] 

   #dN3/dl1
   dNdl1[k][2] =  18.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] + 9.0*_self.GQPoints[k][1]**2\
                 - 9.0*_self.GQPoints[k][1] - 1.0

   #dN4/dl1
   dNdl1[k][3] =  27.0*(- 2.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] - _self.GQPoints[k][1]**2\
                 + _self.GQPoints[k][1])


   # Shape Functions Derivatives in respect to l2
   #dN1/dl2
   dNdl2[k][0] =  18.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] + 9.0*_self.GQPoints[k][0]**2\
                 - 9.0*_self.GQPoints[k][0] 

   #dN2/dl2
   dNdl2[k][1] =  18.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] + 9.0*_self.GQPoints[k][0]**2\
                 - 9.0*_self.GQPoints[k][0] + 1.0 

   #dN3/dl2
   dNdl2[k][2] =  18.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] + 9.0*_self.GQPoints[k][0]**2\
                 - 9.0*_self.GQPoints[k][0] - 1.0

   #dN4/dl2
   dNdl2[k][3] =  27.0*(- 2.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] - _self.GQPoints[k][0]**2\
                 + _self.GQPoints[k][0])


   # Coordinate Transfomation
   # Lewis pag. 64 Eq. 3.108 for 1D
   # Lewis pag. 66 Eq. 3.121 for 2D
   for i in range(0,_self.NUMNODE):
    ii = _self.IEN[_e][i]
    
    dxdl1[k] += _self.x[ii]*dNdl1[k][i]   # dx/dl1 = x1*dN1/dl1 + x2*dN2/dl1 + x3*dN3/dl1 ...
    dydl1[k] += _self.y[ii]*dNdl1[k][i]   # dy/dl1 = y1*dN1/dl1 + y2*dN2/dl1 + y3*dN3/dl1 ...
    dxdl2[k] += _self.x[ii]*dNdl2[k][i]   # dx/dl2 = x1*dN1/dl2 + x2*dN2/dl2 + x3*dN3/dl2 ...
    dydl2[k] += _self.y[ii]*dNdl2[k][i]   # dy/dl2 = y1*dN1/dl2 + y2*dN2/dl2 + y3*dN3/dl2 ...

   # Jacobian Matrix
   # Lewis pag. 66 Eq. 3.121
   J[0][0] = dxdl1[k]
   J[0][1] = dydl1[k]
   J[1][0] = dxdl2[k]
   J[1][1] = dydl2[k]

   jacobian[k] = np.linalg.det(J)


   # Shape Functions Derivatives in respect to x and y
   # Lewis pag. 65 Eq. 3.116
   for i in range(0,_self.NUMNODE):
    dNdx[k][i] = (1.0/jacobian[k])*( dNdl1[k][i]*dydl2[k] - dNdl2[k][i]*dydl1[k])
    dNdy[k][i] = (1.0/jacobian[k])*(-dNdl1[k][i]*dxdl2[k] + dNdl2[k][i]*dxdl1[k])



  _self.mass = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kxx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kxy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kyx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kyy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  

  # Elementary Matrices
  # P.S: It is divided by 2 due to relation 
  # of parallelogram and triangle areas --> DxDy = (jacobian*Dl1*Dl2)/2
  for k in range(0,_self.NUMGAUSS): 
   for i in range(0,_self.NUMNODE):
    for j in range(0,_self.NUMNODE):
    
     _self.mass[i][j] += N[k][i]*N[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.kxx[i][j] += dNdx[k][i]*dNdx[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kxy[i][j] += dNdx[k][i]*dNdy[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kyx[i][j] += dNdy[k][i]*dNdx[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kyy[i][j] += dNdy[k][i]*dNdy[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.gx[i][j] += dNdx[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.gy[i][j] += dNdy[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.dx[i][j] += dNdx[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.dy[i][j] += dNdy[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
  





 def quadratic(_self,_e):
  _self.NUMNODE = 6  #Quadratic Triangular Element - 6 Nodes
  
#  _self.N = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

#  _self.dNdl1 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
#  _self.dNdl2 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  
  dxdl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dxdl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dydl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dydl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)

  dNdx = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  dNdy = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  J = np.zeros([2,2], dtype = float)
  jacobian = np.zeros([_self.NUMGAUSS,1], dtype = float)


  _self.N = np.array([[-0.124998982535, -0.124998982535, 0.001430579518, 0.248575525272, 0.499995930140, 0.499995930140],
                      [-0.124998982535, 0.001430579518, -0.124998982535, 0.499995930140, 0.499995930140, 0.248575525272],
                      [0.001430579518, -0.124998982535, -0.124998982535, 0.499995930140, 0.248575525272, 0.499995930140],
                      [-0.055128566992, -0.055128566992, 0.653307703047, 0.015920894998, 0.220514267970, 0.220514267970],
                      [-0.055128566992, 0.653307703047, -0.055128566992, 0.220514267970, 0.220514267970, 0.015920894998],
                      [0.653307703047, -0.055128566992, -0.055128566992, 0.220514267970, 0.015920894998, 0.220514267970],
                      [-0.117715163308, 0.173768363654, -0.047496257199, 0.790160442766, 0.135307828169, 0.065974785919],
                      [0.173768363654, -0.047496257199, -0.117715163308, 0.135307828169, 0.065974785919, 0.790160442766],
                      [-0.047496257199, -0.117715163308, 0.173768363654, 0.065974785919, 0.790160442766, 0.135307828169],
                      [0.173768363654, -0.117715163308, -0.047496257199, 0.790160442766, 0.065974785919, 0.135307828169],
                      [-0.117715163308, -0.047496257199, 0.173768363654, 0.065974785919, 0.135307828169, 0.790160442766],
                      [-0.047496257199, 0.173768363654, -0.117715163308, 0.135307828169, 0.790160442766, 0.065974785919]])


  _self.dNdl1 = np.array([[-0.002853019316, 0.000000000000, -1.005706038632, 0.997146980684, -0.997146980684, 1.008559057948],
                          [-0.002853019316, 0.000000000000, 0.002853019316, 2.005706038632, -2.005706038632, 0.000000000000 ],
                          [1.005706038632, 0.000000000000, 0.002853019316, 0.997146980684,  -0.997146980684, -1.008559057948],
                          [-0.747643942034, 0.000000000000, -2.495287884068, 0.252356057966, -0.252356057966, 3.242931826102],
                          [-0.747643942034, 0.000000000000, 0.747643942034, 3.495287884068,  -3.495287884068, -0.00000000000],
                          [2.495287884068, 0.000000000000, 0.747643942034, 0.252356057966, -0.252356057966, -3.242931826102 ],
                          [0.241409804136, 0.000000000000, 0.787419800620, 2.546009996484, -2.546009996484, -1.028829604756 ],
                          [1.546009996484, 0.000000000000, -0.241409804137, 0.212580199379, -0.212580199379, -1.304600192347],
                          [-0.787419800621, 0.000000000000, -1.546009996485, 1.241409804136, -1.241409804136, 2.333429797106],
                          [1.546009996484, 0.000000000000, 0.787419800620, 1.241409804136, -1.241409804136, -2.333429797104 ],
                          [0.241409804136, 0.000000000000, -1.546009996485, 0.212580199379, -0.212580199379, 1.304600192349 ],
                          [-0.787419800621, 0.000000000000, -0.241409804137, 2.546009996484, -2.546009996484, 1.028829604758]])


  _self.dNdl2 = np.array([[0.000000000000, -0.002853019316, -1.005706038632, 0.997146980684, 1.008559057948, -0.997146980684],
                          [0.000000000000, 1.005706038632, 0.002853019316, 0.997146980684, -1.008559057948, -0.997146980684 ],
                          [0.000000000000, -0.002853019316, 0.002853019316, 2.005706038632, -0.000000000000, -2.005706038632],
                          [0.000000000000, -0.747643942034, -2.495287884068, 0.252356057966, 3.242931826102, -0.252356057966],
                          [0.000000000000, 2.495287884068, 0.747643942034, 0.252356057966, -3.242931826102, -0.252356057966 ],
                          [0.000000000000, -0.747643942034, 0.747643942034, 3.495287884068, -0.000000000000, -3.495287884068],
                          [0.000000000000, 1.546009996484, 0.787419800620, 1.241409804136, -2.333429797104, -1.241409804136 ],
                          [0.000000000000, -0.787419800621, -0.241409804137, 2.546009996484, 1.028829604758, -2.546009996484],
                          [0.000000000000, 0.241409804136, -1.546009996485, 0.212580199379, 1.304600192349, -0.212580199379 ],
                          [0.000000000000, 0.241409804136, 0.787419800620, 2.546009996484, -1.028829604756, -2.546009996484 ],
                          [0.000000000000, -0.787419800621, -1.546009996485, 1.241409804136, 2.333429797106, -1.241409804136],
                          [0.000000000000, 1.546009996484, -0.241409804137, 0.212580199379, -1.304600192347, -0.212580199379]])

  # A loop is required for each pair of coordinates (l1,l2)
  for k in range(0,_self.NUMGAUSS):
    
   # Area Coordinates
   # Lewis pag. 67 Eq. 3.129
#   L1 = 1.0 - _self.GQPoints[k][0] - _self.GQPoints[k][1]   #L1 = 1 - l1 - l2
#   L2 = _self.GQPoints[k][0]                                #L2 = l1
#   L3 = _self.GQPoints[k][1]                                #L3 = l2

   # Shape Functions
   # Lewis pag. 67 Eq. 3.130
#   _self.N[k][0] = L1*(2.0*L1 - 1.0)  #N1 = L1*(2*L1-1)
#   _self.N[k][1] = L2*(2.0*L2 - 1.0)  #N2 = L2*(2*L2-1)
#   _self.N[k][2] = L3*(2.0*L3 - 1.0)  #N3 = L3*(2*L3-1)
#   _self.N[k][3] = 4.0*L1*L2          #N4 = 4*L1*L2
#   _self.N[k][4] = 4.0*L2*L3          #N5 = 4*L2*L3
#   _self.N[k][5] = 4.0*L3*L1          #N6 = 4*L3*L1

   # Shape Functions Derivatives in respect to l1
#   _self.dNdl1[k][0] = -3.0 + 4.0*_self.GQPoints[k][0] + 4.0*_self.GQPoints[k][1]   #dN1/dl1
#   _self.dNdl1[k][1] =  4.0*_self.GQPoints[k][0] - 1.0                              #dN2/dl1
#   _self.dNdl1[k][2] =  0.0                                                         #dN3/dl1
#   _self.dNdl1[k][3] =  4.0 - 8.0*_self.GQPoints[k][0] - 4.0*_self.GQPoints[k][1]   #dN4/dl1
#   _self.dNdl1[k][4] =  4.0*_self.GQPoints[k][1]                                    #dN5/dl1
#   _self.dNdl1[k][5] = -4.0*_self.GQPoints[k][1]                                    #dN6/dl1


   # Shape Functions Derivatives in respect to l2
#   _self.dNdl2[k][0] = -3.0 + 4.0*_self.GQPoints[k][0] + 4.0*_self.GQPoints[k][1]   #dN1/dl2
#   _self.dNdl2[k][1] =  0.0                                                         #dN2/dl2
#   _self.dNdl2[k][2] =  4.0*_self.GQPoints[k][1] - 1.0                              #dN3/dl2
#   _self.dNdl2[k][3] = -4.0*_self.GQPoints[k][0]                                    #dN4/dl2
#   _self.dNdl2[k][4] =  4.0*_self.GQPoints[k][0]                                    #dN5/dl2
#   _self.dNdl2[k][5] =  4.0 - 4.0*_self.GQPoints[k][0] - 8.0*_self.GQPoints[k][1]   #dN6/dl2

   # Coordinate Transfomation
   # Lewis pag. 64 Eq. 3.108 for 1D
   # Lewis pag. 66 Eq. 3.121 for 2D
   for i in range(0,_self.NUMNODE):
    ii = _self.IEN[_e][i]
    
    dxdl1[k] += _self.x[ii]*_self.dNdl1[k][i]   # dx/dl1 = x1*dN1/dl1 + x2*dN2/dl1 + x3*dN3/dl1 ...
    dydl1[k] += _self.y[ii]*_self.dNdl1[k][i]   # dy/dl1 = y1*dN1/dl1 + y2*dN2/dl1 + y3*dN3/dl1 ...
    dxdl2[k] += _self.x[ii]*_self.dNdl2[k][i]   # dx/dl2 = x1*dN1/dl2 + x2*dN2/dl2 + x3*dN3/dl2 ...
    dydl2[k] += _self.y[ii]*_self.dNdl2[k][i]   # dy/dl2 = y1*dN1/dl2 + y2*dN2/dl2 + y3*dN3/dl2 ...


   # Jacobian Matrix
   # Lewis pag. 66 Eq. 3.121
   J[0][0] = dxdl1[k]
   J[0][1] = dydl1[k]
   J[1][0] = dxdl2[k]
   J[1][1] = dydl2[k]

   jacobian[k] = np.linalg.det(J)


   # Shape Functions Derivatives in respect to x and y
   # Lewis pag. 65 Eq. 3.116
   for i in range(0,_self.NUMNODE):
    dNdx[k][i] = (1.0/jacobian[k])*( _self.dNdl1[k][i]*dydl2[k] - _self.dNdl2[k][i]*dydl1[k])
    dNdy[k][i] = (1.0/jacobian[k])*(-_self.dNdl1[k][i]*dxdl2[k] + _self.dNdl2[k][i]*dxdl1[k])



  _self.mass = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kxx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kxy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kyx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kyy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  

  # Elementary Matrices
  # P.S: It is divided by 2 due to relation 
  # of parallelogram and triangle areas --> DxDy = (jacobian*Dl1*Dl2)/2
  for k in range(0,_self.NUMGAUSS): 
   for i in range(0,_self.NUMNODE):
    for j in range(0,_self.NUMNODE):
    
     _self.mass[i][j] += _self.N[k][i]*_self.N[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.kxx[i][j] += dNdx[k][i]*dNdx[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kxy[i][j] += dNdx[k][i]*dNdy[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kyx[i][j] += dNdy[k][i]*dNdx[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kyy[i][j] += dNdy[k][i]*dNdy[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.gx[i][j] += dNdx[k][j]*_self.N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.gy[i][j] += dNdy[k][j]*_self.N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.dx[i][j] += dNdx[k][j]*_self.N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.dy[i][j] += dNdy[k][j]*_self.N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
  



 def cubic(_self,_e):
  _self.NUMNODE = 10  #Cubic Triangular Element - 10 Nodes
  
  N = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  dNdl1 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  dNdl2 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  
  dxdl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dxdl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dydl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dydl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)

  dNdx = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  dNdy = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  J = np.zeros([2,2], dtype = float)
  jacobian = np.zeros([_self.NUMGAUSS,1], dtype = float)




  # A loop is required for each pair of coordinates (l1,l2)
  for k in range(0,_self.NUMGAUSS):
    
   # Area Coordinates
   L1 = _self.GQPoints[k][0]                                #L1 = l1
   L2 = _self.GQPoints[k][1]                                #L2 = l2
   L3 = 1.0 - _self.GQPoints[k][0] - _self.GQPoints[k][1]   #L3 = 1 - l1 - l2

   # Shape Functions
   # Lewis pag. 56 and 57 Eq. 3.81/3.82/3.83/3.84/3.85
   N[k][0] = 0.5*L1*(3.0*L1 - 1.0)*(3.0*L1 - 2.0)   #N1  = 1/2*L1*(3*L1-1)*(3*L1-2)
   N[k][1] = 0.5*L2*(3.0*L2 - 1.0)*(3.0*L2 - 2.0)   #N2  = 1/2*L2*(3*L2-1)*(3*L2-2)
   N[k][2] = 0.5*L3*(3.0*L3 - 1.0)*(3.0*L3 - 2.0)   #N3  = 1/2*L3*(3*L3-1)*(3*L3-2)
   N[k][3] = 4.5*L1*L2*(3.0*L1 - 1.0)               #N4  = 9/2*L1*L2(3*L1-1)
   N[k][4] = 4.5*L1*L2*(3.0*L2 - 1.0)               #N5  = 9/2*L1*L2(3*L2-1)
   N[k][5] = 4.5*L2*L3*(3.0*L2 - 1.0)               #N6  = 9/2*L2*L3(3*L2-1)
   N[k][6] = 4.5*L2*L3*(3.0*L3 - 1.0)               #N7  = 9/2*L2*L3(3*L3-1)
   N[k][7] = 4.5*L3*L1*(3.0*L3 - 1.0)               #N8  = 9/2*L3*L1(3*L3-1)
   N[k][8] = 4.5*L3*L1*(3.0*L1 - 1.0)               #N9  = 9/2*L3*L1(3*L1-1)
   N[k][9] = 27.0*L1*L2*L3                          #N10 = 27*L1*L2*L3


   # Shape Functions Derivatives in respect to l1
   #dN1/dl1
   dNdl1[k][0] =  0.5*( 27.0*_self.GQPoints[k][0]**2 - 18.0*_self.GQPoints[k][0] + 2.0)

   #dN2/dl1
   dNdl1[k][1] =  0.0

   #dN3/dl1
   dNdl1[k][2] =  0.5*(-11.0 + 36.0*_self.GQPoints[k][1] - 54.0*_self.GQPoints[k][0]*_self.GQPoints[k][1]\
                      - 27.0*_self.GQPoints[k][1]**2 + 36.0*_self.GQPoints[k][0]\
                      - 27.0*_self.GQPoints[k][0]**2)
   #dN4/dl1
   dNdl1[k][3] = 4.5*( 6.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] - _self.GQPoints[k][1])

   #dN5/dl1
   dNdl1[k][4] = 4.5*( 3.0*_self.GQPoints[k][1]**2 - _self.GQPoints[k][1])

   #dN6/dl1
   dNdl1[k][5] = 4.5*(-3.0*_self.GQPoints[k][1]**2 + _self.GQPoints[k][1])

   #dN7/dl1
   dNdl1[k][6] = 4.5*(-5.0*_self.GQPoints[k][1] + 6.0*_self.GQPoints[k][0]*_self.GQPoints[k][1]\
                     + 6.0*_self.GQPoints[k][1]**2)

   #dN8/dl1
   dNdl1[k][7] = 4.5*( 2.0 - 5.0*_self.GQPoints[k][1] - 10.0*_self.GQPoints[k][0]\
                     + 3.0*_self.GQPoints[k][1]**2 + 12.0*_self.GQPoints[k][0]*_self.GQPoints[k][1]\
                     + 9.0*_self.GQPoints[k][0]**2)

   #dN9/dl1
   dNdl1[k][8] = 4.5*(-1.0 + _self.GQPoints[k][1] + 8.0*_self.GQPoints[k][0]\
                     - 6.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] - 9.0*_self.GQPoints[k][0]**2)

   #dN10/dl1
   dNdl1[k][9] = 27.0*(_self.GQPoints[k][1] - 2.0*_self.GQPoints[k][0]*_self.GQPoints[k][1]\
                     - _self.GQPoints[k][1]**2)



   # Shape Functions Derivatives in respect to l2
   #dN1/dl2
   dNdl2[k][0] =  0.0

   #dN2/dl2
   dNdl2[k][1] =  0.5*( 27.0*_self.GQPoints[k][1]**2 - 18.0*_self.GQPoints[k][1] + 2.0)

   #dN3/dl2
   dNdl2[k][2] =  0.5*(-11.0 + 36.0*_self.GQPoints[k][0] - 54.0*_self.GQPoints[k][0]*_self.GQPoints[k][1]\
                      - 27.0*_self.GQPoints[k][0]**2 + 36.0*_self.GQPoints[k][1]\
                      - 27.0*_self.GQPoints[k][1]**2)

   #dN4/dl2
   dNdl2[k][3] = 4.5*( 3.0*_self.GQPoints[k][0]**2 - _self.GQPoints[k][0])

   #dN5/dl2
   dNdl2[k][4] = 4.5*( 6.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] - _self.GQPoints[k][0])

   #dN6/dl2
   dNdl2[k][5] = 4.5*(-1.0 + _self.GQPoints[k][0] + 8.0*_self.GQPoints[k][1]\
                     - 6.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] - 9.0*_self.GQPoints[k][1]**2)

   #dN7/dl2
   dNdl2[k][6] = 4.5*( 2.0 - 5.0*_self.GQPoints[k][0] - 10.0*_self.GQPoints[k][1]\
                     + 3.0*_self.GQPoints[k][0]**2 + 12.0*_self.GQPoints[k][0]*_self.GQPoints[k][1]\
                     + 9.0*_self.GQPoints[k][1]**2)

   #dN8/dl2
   dNdl2[k][7] = 4.5*(-5.0*_self.GQPoints[k][0] + 6.0*_self.GQPoints[k][0]*_self.GQPoints[k][1]\
                     + 6.0*_self.GQPoints[k][0]**2)

   #dN9/dl2
   dNdl2[k][8] = 4.5*(-3.0*_self.GQPoints[k][0]**2 + _self.GQPoints[k][0])

   #dN10/dl2
   dNdl2[k][9] = 27.0*(_self.GQPoints[k][0] - 2.0*_self.GQPoints[k][0]*_self.GQPoints[k][1]\
                     - _self.GQPoints[k][0]**2)



   # Coordinate Transfomation
   # Lewis pag. 64 Eq. 3.108 for 1D
   # Lewis pag. 66 Eq. 3.121 for 2D
   for i in range(0,_self.NUMNODE):
    ii = _self.IEN[_e][i]
    
    dxdl1[k] += _self.x[ii]*dNdl1[k][i]   # dx/dl1 = x1*dN1/dl1 + x2*dN2/dl1 + x3*dN3/dl1 ...
    dydl1[k] += _self.y[ii]*dNdl1[k][i]   # dy/dl1 = y1*dN1/dl1 + y2*dN2/dl1 + y3*dN3/dl1 ...
    dxdl2[k] += _self.x[ii]*dNdl2[k][i]   # dx/dl2 = x1*dN1/dl2 + x2*dN2/dl2 + x3*dN3/dl2 ...
    dydl2[k] += _self.y[ii]*dNdl2[k][i]   # dy/dl2 = y1*dN1/dl2 + y2*dN2/dl2 + y3*dN3/dl2 ...



   # Jacobian Matrix
   # Lewis pag. 66 Eq. 3.121
   J[0][0] = dxdl1[k]
   J[0][1] = dydl1[k]
   J[1][0] = dxdl2[k]
   J[1][1] = dydl2[k]

   jacobian[k] = np.linalg.det(J)


   # Shape Functions Derivatives in respect to x and y
   # Lewis pag. 65 Eq. 3.116
   for i in range(0,_self.NUMNODE):
    dNdx[k][i] = (1.0/jacobian[k])*( dNdl1[k][i]*dydl2[k] - dNdl2[k][i]*dydl1[k])
    dNdy[k][i] = (1.0/jacobian[k])*(-dNdl1[k][i]*dxdl2[k] + dNdl2[k][i]*dxdl1[k])



  _self.mass = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kxx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kxy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kyx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kyy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dx = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dy = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  

  # Elementary Matrices
  # P.S: It is divided by 2 due to relation 
  # of parallelogram and triangle areas --> DxDy = (jacobian*Dl1*Dl2)/2
  for k in range(0,_self.NUMGAUSS): 
   for i in range(0,_self.NUMNODE):
    for j in range(0,_self.NUMNODE):
    
     _self.mass[i][j] += N[k][i]*N[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.kxx[i][j] += dNdx[k][i]*dNdx[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kxy[i][j] += dNdx[k][i]*dNdy[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kyx[i][j] += dNdy[k][i]*dNdx[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kyy[i][j] += dNdy[k][i]*dNdy[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.gx[i][j] += dNdx[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.gy[i][j] += dNdy[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.dx[i][j] += dNdx[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.dy[i][j] += dNdy[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0





 # Reference Analytic: Fundamentals of the Finite
 #                     Element Method for Heat Transfer
 #                     and Fluid Flow - Lewis, Nithiarasu,
 #                     Seetharamu - pg. 196-200
 #                     For Q_elem pg. 126
 #                     For 1D pg. 193
 def analytic(_self, _e):

  i = _self.IEN[_e][0]
  j = _self.IEN[_e][1]
  k = _self.IEN[_e][2]

  bi = _self.y[j]-_self.y[k]
  bj = _self.y[k]-_self.y[i]
  bk = _self.y[i]-_self.y[j]
  ci = _self.x[k]-_self.x[j]
  cj = _self.x[i]-_self.x[k]
  ck = _self.x[j]-_self.x[i]


  A = 0.5*np.linalg.det(np.array([[1, _self.x[i], _self.y[i]],
 				  [1, _self.x[j], _self.y[j]],
				  [1, _self.x[k], _self.y[k]]]))


  _self.mass = (A/12.)*np.array([[2.,1.,1.],
                                 [1.,2.,1.],
                                 [1.,1.,2.]])

  _self.q = (A/3.)*np.ones([3,1], dtype = float)
  
  _self.gx = (1./6)*np.array([[bi,bj,bk],
                              [bi,bj,bk],
                              [bi,bj,bk]]) 
   
  _self.gy = (1./6)*np.array([[ci,cj,ck],
                              [ci,cj,ck],
                              [ci,cj,ck]])

  _self.kxx = (1./(4*A))*np.array([[bi*bi,bi*bj,bi*bk],
                                   [bj*bi,bj*bj,bj*bk],
                                   [bk*bi,bk*bj,bk*bk]])

  _self.kyy = (1./(4*A))*np.array([[ci*ci,ci*cj,ci*ck],
                                   [cj*ci,cj*cj,cj*ck],
                                   [ck*ci,ck*cj,ck*ck]])

  _self.kxy = (1./(4*A))*np.array([[bi*ci,bi*cj,bi*ck],
                                   [bj*ci,bj*cj,bj*ck],
                                   [bk*ci,bk*cj,bk*ck]])

  _self.kyx = (1./(4*A))*np.array([[ci*bi,ci*bj,ci*bk],
                                   [cj*bi,cj*bj,cj*bk],
                                   [ck*bi,ck*bj,ck*bk]])



 def axisymmetric(_self, _e):
  _self.r = _self.y
  _self.z = _self.x

  i = _self.IEN[_e][0]
  j = _self.IEN[_e][1]
  k = _self.IEN[_e][2]

  bi = _self.z[j] - _self.z[k]
  bj = _self.z[k] - _self.z[i]
  bk = _self.z[i] - _self.z[j]
  ci = _self.r[k] - _self.r[j]
  cj = _self.r[i] - _self.r[k]
  ck = _self.r[j] - _self.r[i]

  A = 0.5*np.linalg.det(np.array([[1, _self.r[i], _self.z[i]],
 				  [1, _self.r[j], _self.z[j]],
				  [1, _self.r[k], _self.z[k]]]))

  r = (_self.r[i] + _self.r[j] + _self.r[k])/3.

  r_vec = np.array([[_self.r[i]],
                    [_self.r[j]],
                    [_self.r[k]]])

  _self.M_elem = (A/12.)*np.array([[2.,1.,1.],
				   [1.,2.,1.],
				   [1.,1.,2.]])

  _self.Q_elem = (2*np.pi)*np.dot(_self.M_elem,r_vec)
  
  _self.Gr_elem = (1./6)*np.array([[bi,bj,bk],
                                   [bi,bj,bk],
                                   [bi,bj,bk]]) 
   
  _self.Gz_elem = (1./6)*np.array([[ci,cj,ck],
                                   [ci,cj,ck],
                                   [ci,cj,ck]])

  _self.Kr_elem = ((2*np.pi*r)/(4*A))*np.array([[bi*bi,bj*bi,bk*bi],
                                                [bi*bj,bj*bj,bk*bj],
                                                [bi*bk,bj*bk,bk*bk]])

  _self.Kz_elem = ((2*np.pi*r)/(4*A))*np.array([[ci*ci,cj*ci,ck*ci],
                                                [ci*cj,cj*cj,ck*cj],
                                                [ci*ck,cj*ck,ck*ck]])






