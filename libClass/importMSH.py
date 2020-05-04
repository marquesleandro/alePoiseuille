# ==========================================
# Code created by Leandro Marques at 12/2018
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used to import .msh
# for n governament equations 


import numpy as np


class Linear1D:
 def __init__(_self, _dir, _file, _number_equations):
  _self.name = _dir + '/' + _file
  _self.neq = _number_equations
  _self.gmsh = []
  _self.neumann_lines = {}
  _self.neumann_pts = {}
  _self.dirichlet_lines = {}
  _self.dirichlet_pts = {}
  _self.neighbors_nodes = {}
  _self.neighbors_elements = {}
  _self.far_neighbors_nodes = {}
  _self.far_neighbors_elements = {}

  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.gmsh.append(row[:])

  
  _self.nphysical = int(_self.gmsh[4][0])
  _self.npoints = int(_self.gmsh[_self.nphysical+7][0])
  
  for i in range(1, _self.neq + 1):
   _self.neumann_lines[i] = []
   _self.neumann_pts[i] = []
   _self.dirichlet_lines[i] = []
   _self.dirichlet_pts[i] = []

  # Lines classification in neumann or dirichlet
  for i in range(0,_self.nphysical):
   for j in range(1,_self.neq + 1):
    if _self.gmsh[5+i][2] == '"neumann' + str(j):
     _self.neumann_lines[j].append(int(_self.gmsh[5+i][1]))

    elif _self.gmsh[5+i][2] == '"dirichlet' + str(j):
     _self.dirichlet_lines[j].append(int(_self.gmsh[5+i][1]))


  jj = _self.nphysical + _self.npoints 
  ii = 1
  element_type = int(_self.gmsh[(jj + 10 + ii)][1])

  # Assembly of the neumann edges and dirichlet points
  while element_type == 15:
    a_1 = int(_self.gmsh[(jj + 10 + ii)][3])
    a_2 = int(_self.gmsh[(jj + 10 + ii)][5])
    a_3 = [a_1,a_2]
   
    for i in range(1, _self.neq + 1): 
     if a_1 in _self.neumann_lines[i]:
      _self.neumann_pts[i].append(a_3)
 
      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])
 
     elif a_1 in _self.dirichlet_lines[i]:
      _self.dirichlet_pts[i].append(a_3)

      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])

   
  for i in range(1, _self.neq + 1):
   _self.neumann_pts[i] = np.array(_self.neumann_pts[i]) 
   _self.dirichlet_pts[i] = np.array(_self.dirichlet_pts[i]) 


  _self.nelem = int(_self.gmsh[jj + 10][0]) - ii + 1

  for i in range(0, _self.npoints):
   _self.neighbors_nodes[i] = []
   _self.neighbors_elements[i] = []
   _self.far_neighbors_nodes[i] = []
   _self.far_neighbors_elements[i] = []


  _self.ii = ii
  _self.jj = jj


 def coord(_self):
  _self.x = np.zeros([_self.npoints,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.npoints):  
   _self.x[i] = _self.gmsh[_self.nphysical + 8 + i][1]
   _self.npts.append(i)


 def ien(_self):
  _self.IEN = np.zeros([_self.nelem,2], dtype = int)
  _self.GL = len(_self.IEN[0,:])
  length = [] 

  for e in range(0, _self.nelem):
   v1 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][5]) - 1
   v2 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][6]) - 1
  
   _self.IEN[e] = [v1,v2]
  
   _self.neighbors_nodes[v1].extend(_self.IEN[e])  
   _self.neighbors_nodes[v2].extend(_self.IEN[e])  
   
   _self.neighbors_nodes[v1] = list(set(_self.neighbors_nodes[v1]))
   _self.neighbors_nodes[v2] = list(set(_self.neighbors_nodes[v2]))
   
   _self.neighbors_elements[v1].append(e)  
   _self.neighbors_elements[v2].append(e)  

   x_a = _self.x[v1] - _self.x[v2]
   length1 = np.sqrt(x_a**2)
   length.append(length1)

  _self.length_min = min(length)


  for i in range(0, _self.npoints):
   for j in _self.neighbors_nodes[i]:
    _self.far_neighbors_nodes[i].extend(_self.neighbors_nodes[j]) 
    _self.far_neighbors_elements[i].extend(_self.neighbors_elements[j]) 
 
   _self.far_neighbors_nodes[i] = list(set(_self.far_neighbors_nodes[i])\
                                     - set(_self.neighbors_nodes[i]))
   
   _self.far_neighbors_elements[i] = list(set(_self.far_neighbors_elements[i])\
                                        - set(_self.neighbors_elements[i]))



class Quad1D:
 def __init__(_self, _dir, _file, _number_equations):
  _self.name = _dir + '/' + _file
  _self.neq = _number_equations
  _self.gmsh = []
  _self.neumann_lines = {}
  _self.neumann_pts = {}
  _self.dirichlet_lines = {}
  _self.dirichlet_pts = {}
  _self.neighbors_nodes = {}
  _self.neighbors_elements = {}
  _self.far_neighbors_nodes = {}
  _self.far_neighbors_elements = {}

  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.gmsh.append(row[:])

  
  _self.nphysical = int(_self.gmsh[4][0])
  _self.npoints = int(_self.gmsh[_self.nphysical+7][0])
  
  for i in range(1, _self.neq + 1):
   _self.neumann_lines[i] = []
   _self.neumann_pts[i] = []
   _self.dirichlet_lines[i] = []
   _self.dirichlet_pts[i] = []

  # Lines classification in neumann or dirichlet
  for i in range(0,_self.nphysical):
   for j in range(1,_self.neq + 1):
    if _self.gmsh[5+i][2] == '"neumann' + str(j):
     _self.neumann_lines[j].append(int(_self.gmsh[5+i][1]))

    elif _self.gmsh[5+i][2] == '"dirichlet' + str(j):
     _self.dirichlet_lines[j].append(int(_self.gmsh[5+i][1]))


  jj = _self.nphysical + _self.npoints 
  ii = 1
  element_type = int(_self.gmsh[(jj + 10 + ii)][1])

  # Assembly of the neumann edges and dirichlet points
  while element_type == 15:
    a_1 = int(_self.gmsh[(jj + 10 + ii)][3])
    a_2 = int(_self.gmsh[(jj + 10 + ii)][5])
    a_3 = [a_1,a_2]
   
    for i in range(1, _self.neq + 1): 
     if a_1 in _self.neumann_lines[i]:
      _self.neumann_pts[i].append(a_3)
 
      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])
 
     elif a_1 in _self.dirichlet_lines[i]:
      _self.dirichlet_pts[i].append(a_3)

      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])

   
  for i in range(1, _self.neq + 1):
   _self.neumann_pts[i] = np.array(_self.neumann_pts[i]) 
   _self.dirichlet_pts[i] = np.array(_self.dirichlet_pts[i]) 


  _self.nelem = int(_self.gmsh[jj + 10][0]) - ii + 1

  for i in range(0, _self.npoints):
   _self.neighbors_nodes[i] = []
   _self.neighbors_elements[i] = []
   _self.far_neighbors_nodes[i] = []
   _self.far_neighbors_elements[i] = []


  _self.ii = ii
  _self.jj = jj


 def coord(_self):
  _self.x = np.zeros([_self.npoints,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.npoints):  
   _self.x[i] = _self.gmsh[_self.nphysical + 8 + i][1]
   _self.npts.append(i)


 def ien(_self):
  _self.IEN = np.zeros([_self.nelem,3], dtype = int)
  _self.GL = len(_self.IEN[0,:])
  length = [] 

  for e in range(0, _self.nelem):
   v1 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][5]) - 1
   v2 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][6]) - 1
   v3 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][7]) - 1
  
   _self.IEN[e] = [v1,v2,v3]
  
   _self.neighbors_nodes[v1].extend(_self.IEN[e])  
   _self.neighbors_nodes[v2].extend(_self.IEN[e])  
   _self.neighbors_nodes[v3].extend(_self.IEN[e])  
   
   _self.neighbors_nodes[v1] = list(set(_self.neighbors_nodes[v1]))
   _self.neighbors_nodes[v2] = list(set(_self.neighbors_nodes[v2]))
   _self.neighbors_nodes[v3] = list(set(_self.neighbors_nodes[v2]))

   _self.neighbors_elements[v1].append(e)  
   _self.neighbors_elements[v2].append(e)  
   _self.neighbors_elements[v3].append(e)  

   x_a = _self.x[v1] - _self.x[v2]
   length1 = np.sqrt(x_a**2)
   length.append(length1)

  _self.length_min = min(length)


  for i in range(0, _self.npoints):
   for j in _self.neighbors_nodes[i]:
    _self.far_neighbors_nodes[i].extend(_self.neighbors_nodes[j]) 
    _self.far_neighbors_elements[i].extend(_self.neighbors_elements[j]) 
 
   _self.far_neighbors_nodes[i] = list(set(_self.far_neighbors_nodes[i])\
                                     - set(_self.neighbors_nodes[i]))
   
   _self.far_neighbors_elements[i] = list(set(_self.far_neighbors_elements[i])\
                                        - set(_self.neighbors_elements[i]))




class Linear2D:
 def __init__(_self, _dir, _file):
  _self.name = _dir + '/' + _file
  _self.mshFile = []
  _self.boundaryEdges = []
  _self.boundaryNodes = []
  _self.neighborsNodes = {}
  _self.neighborsNodesALE = {}
  _self.neighborsElements = {}


  
  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.mshFile.append(row[:])

  _self.numPhysical = int(_self.mshFile[4][0])
  _self.numNodes = int(_self.mshFile[_self.numPhysical+7][0])
 

  # Boundary edges Assembly
  countLineStart = _self.numPhysical + _self.numNodes 
  countLine = 1
  mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])

  while mshFileElementType == 1:
   line = int(_self.mshFile[(countLineStart + 10 + countLine)][3])
   v1 = int(_self.mshFile[(countLineStart + 10 + countLine)][5])
   v2 = int(_self.mshFile[(countLineStart + 10 + countLine)][6])
   
   _self.boundaryEdges.append([line,v1,v2])
   _self.boundaryNodes.append(v1)
   _self.boundaryNodes.append(v2)

   countLine = countLine + 1
   mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])
 
  _self.boundaryEdges = np.array(_self.boundaryEdges) 
  _self.boundaryNodes = list(set(_self.boundaryNodes))
  _self.numElements = int(_self.mshFile[countLineStart + 10][0]) - countLine + 1



  # Coordinate vectors Assembly
  _self.x = np.zeros([_self.numNodes,1], dtype = float)
  _self.y = np.zeros([_self.numNodes,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.numNodes):  
   _self.x[i] = _self.mshFile[_self.numPhysical + 8 + i][1]
   _self.y[i] = _self.mshFile[_self.numPhysical + 8 + i][2]
   _self.neighborsNodes[i] = []
   _self.neighborsNodesALE[i] = []
   _self.neighborsElements[i] = []
   _self.npts.append(i)



  # IEN matrix and Neighbors list Assembly
  _self.IEN = np.zeros([_self.numElements,3], dtype = int)
  _self.FreedomDegree = len(_self.IEN[0,:])
  length = []

  for e in range(0, _self.numElements):
   v1 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][5]) - 1
   v2 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][6]) - 1
   v3 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][7]) - 1
  
   _self.IEN[e] = [v1,v2,v3]


  
   _self.neighborsNodes[v1].extend(_self.IEN[e])  
   _self.neighborsNodes[v2].extend(_self.IEN[e])  
   _self.neighborsNodes[v3].extend(_self.IEN[e])  

   _self.neighborsNodes[v1] = list(set(_self.neighborsNodes[v1]))
   _self.neighborsNodes[v2] = list(set(_self.neighborsNodes[v2]))
   _self.neighborsNodes[v3] = list(set(_self.neighborsNodes[v3]))

 

   _self.neighborsNodesALE[v1].extend(_self.IEN[e])  
   _self.neighborsNodesALE[v2].extend(_self.IEN[e])  
   _self.neighborsNodesALE[v3].extend(_self.IEN[e])  
     
   _self.neighborsNodesALE[v1] = list(set(_self.neighborsNodesALE[v1]))
   _self.neighborsNodesALE[v2] = list(set(_self.neighborsNodesALE[v2]))
   _self.neighborsNodesALE[v3] = list(set(_self.neighborsNodesALE[v3]))
   


   _self.neighborsElements[v1].append(e)  
   _self.neighborsElements[v2].append(e)  
   _self.neighborsElements[v3].append(e)  



   x_a = _self.x[v1] - _self.x[v2]
   x_b = _self.x[v2] - _self.x[v3]
   x_c = _self.x[v3] - _self.x[v1]
   
   y_a = _self.y[v1] - _self.y[v2]
   y_b = _self.y[v2] - _self.y[v3]
   y_c = _self.y[v3] - _self.y[v1]
   
   length1 = np.sqrt(x_a**2 + y_a**2)
   length2 = np.sqrt(x_b**2 + y_b**2)
   length3 = np.sqrt(x_c**2 + y_c**2)

   length.append(length1)
   length.append(length2)
   length.append(length3)
   
  _self.minLengthMesh = min(length)





class Mini2D:
 def __init__(_self, _dir, _file, _number_equations):
  _self.name = _dir + '/' + _file
  _self.neq = _number_equations
  _self.gmsh = []
  _self.neumann_lines = {}
  _self.dirichlet_lines = {}
  _self.neumann_edges = {}
  _self.dirichlet_pts = {}
  _self.neighbors_nodes = {}
  _self.neighbors_elements = {}
  _self.far_neighbors_nodes = {}
  _self.far_neighbors_elements = {}

  
  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.gmsh.append(row[:])

  
  _self.nphysical = int(_self.gmsh[4][0])
  _self.npoints = int(_self.gmsh[_self.nphysical+7][0])
  
  for i in range(1, _self.neq + 1):
   _self.neumann_lines[i] = []
   _self.neumann_edges[i] = []
   _self.dirichlet_lines[i] = []
   _self.dirichlet_pts[i] = []


  # Lines classification in neumann or dirichlet
  for i in range(0,_self.nphysical):
   for j in range(1,_self.neq + 1):
    if _self.gmsh[5+i][2] == '"neumann' + str(j):
     _self.neumann_lines[j].append(int(_self.gmsh[5+i][1]))

    elif _self.gmsh[5+i][2] == '"dirichlet' + str(j):
     _self.dirichlet_lines[j].append(int(_self.gmsh[5+i][1]))


  jj = _self.nphysical + _self.npoints 
  ii = 1
  element_type = int(_self.gmsh[(jj + 10 + ii)][1])

  # Assembly of the neumann edges and dirichlet points
  while element_type == 1:
    a_1 = int(_self.gmsh[(jj + 10 + ii)][3])
    a_2 = int(_self.gmsh[(jj + 10 + ii)][5])
    a_3 = int(_self.gmsh[(jj + 10 + ii)][6])
    a_4 = [a_1,a_2,a_3]
   
    for i in range(1, _self.neq + 1): 
     if a_1 in _self.neumann_lines[i]:
      _self.neumann_edges[i].append(a_4)
 
      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])
 
     elif a_1 in _self.dirichlet_lines[i]:
      _self.dirichlet_pts[i].append(a_4)

      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])

   
  for i in range(1, _self.neq + 1):
   _self.neumann_edges[i] = np.array(_self.neumann_edges[i]) 
   _self.dirichlet_pts[i] = np.array(_self.dirichlet_pts[i]) 


  _self.nelem = int(_self.gmsh[jj + 10][0]) - ii + 1
  _self.NP = _self.npoints
  _self.npoints = _self.NP + _self.nelem


  for i in range(0, _self.npoints):
   _self.neighbors_nodes[i] = []
   _self.neighbors_elements[i] = []
   _self.far_neighbors_nodes[i] = []
   _self.far_neighbors_elements[i] = []


  _self.ii = ii
  _self.jj = jj


 def coord(_self):
  _self.x = np.zeros([_self.npoints,1], dtype = float)
  _self.y = np.zeros([_self.npoints,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.NP):  
   _self.x[i] = _self.gmsh[_self.nphysical + 8 + i][1]
   _self.y[i] = _self.gmsh[_self.nphysical + 8 + i][2]
   _self.npts.append(i)


 def ien(_self):
  _self.IEN = np.zeros([_self.nelem,4], dtype = int)
  _self.GL = len(_self.IEN[0,:])
  length = []

  for e in range(0, _self.nelem):
   v1 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][5]) - 1
   v2 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][6]) - 1
   v3 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][7]) - 1
   v4 = _self.NP + e
  
   _self.IEN[e] = [v1,v2,v3,v4]

   _self.x[v4] = (_self.x[v1] + _self.x[v2] + _self.x[v3])/3.0
   _self.y[v4] = (_self.y[v1] + _self.y[v2] + _self.y[v3])/3.0
   _self.npts.append(v4)

   _self.neighbors_nodes[v1].extend(_self.IEN[e])  
   _self.neighbors_nodes[v2].extend(_self.IEN[e])  
   _self.neighbors_nodes[v3].extend(_self.IEN[e])  
   _self.neighbors_nodes[v4].extend(_self.IEN[e])  
   
   _self.neighbors_nodes[v1] = list(set(_self.neighbors_nodes[v1]))
   _self.neighbors_nodes[v2] = list(set(_self.neighbors_nodes[v2]))
   _self.neighbors_nodes[v3] = list(set(_self.neighbors_nodes[v3]))
   _self.neighbors_nodes[v4] = list(set(_self.neighbors_nodes[v4]))
   
   _self.neighbors_elements[v1].append(e)  
   _self.neighbors_elements[v2].append(e)  
   _self.neighbors_elements[v3].append(e)  
   _self.neighbors_elements[v4].append(e)  

   x_a = _self.x[v1] - _self.x[v2]
   x_b = _self.x[v2] - _self.x[v3]
   x_c = _self.x[v3] - _self.x[v1]
   
   y_a = _self.y[v1] - _self.y[v2]
   y_b = _self.y[v2] - _self.y[v3]
   y_c = _self.y[v3] - _self.y[v1]
   
   length1 = np.sqrt(x_a**2 + y_a**2)
   length2 = np.sqrt(x_b**2 + y_b**2)
   length3 = np.sqrt(x_c**2 + y_c**2)

   length.append(length1)
   length.append(length2)
   length.append(length3)
   
  _self.length_min = min(length)

  for i in range(0, _self.npoints):
   for j in _self.neighbors_nodes[i]:
    _self.far_neighbors_nodes[i].extend(_self.neighbors_nodes[j]) 
    _self.far_neighbors_elements[i].extend(_self.neighbors_elements[j]) 
 
   _self.far_neighbors_nodes[i] = list(set(_self.far_neighbors_nodes[i])\
                                     - set(_self.neighbors_nodes[i]))
   
   _self.far_neighbors_elements[i] = list(set(_self.far_neighbors_elements[i])\
                                        - set(_self.neighbors_elements[i]))





class Quad2D:
 def __init__(_self, _dir, _file):
  _self.name = _dir + '/' + _file
  _self.mshFile = []
  _self.boundaryEdges = []
  _self.boundaryNodes = []
  _self.neighborsNodes = {}
  _self.neighborsNodesALE = {}
  _self.neighborsElements = {}


  
  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.mshFile.append(row[:])

  _self.numPhysical = int(_self.mshFile[4][0])
  _self.numNodes = int(_self.mshFile[_self.numPhysical+7][0])
 

  # Boundary edges Assembly
  countLineStart = _self.numPhysical + _self.numNodes 
  countLine = 1
  mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])

  while mshFileElementType == 8:
   line = int(_self.mshFile[(countLineStart + 10 + countLine)][3])
   v1 = int(_self.mshFile[(countLineStart + 10 + countLine)][5])
   v2 = int(_self.mshFile[(countLineStart + 10 + countLine)][6])
   v3 = int(_self.mshFile[(countLineStart + 10 + countLine)][7])
   
   _self.boundaryEdges.append([line,v1,v2,v3])
   _self.boundaryNodes.append(v1)
   _self.boundaryNodes.append(v2)
   _self.boundaryNodes.append(v3)

   countLine = countLine + 1
   mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])
 
  _self.boundaryEdges = np.array(_self.boundaryEdges) 
  _self.boundaryNodes = list(set(_self.boundaryNodes))
  _self.numElements = int(_self.mshFile[countLineStart + 10][0]) - countLine + 1



  # Coordinate vectors Assembly
  _self.x = np.zeros([_self.numNodes,1], dtype = float)
  _self.y = np.zeros([_self.numNodes,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.numNodes):  
   _self.x[i] = _self.mshFile[_self.numPhysical + 8 + i][1]
   _self.y[i] = _self.mshFile[_self.numPhysical + 8 + i][2]
   _self.neighborsNodes[i] = []
   _self.neighborsNodesALE[i] = []
   _self.neighborsElements[i] = []
   _self.npts.append(i)



  # IEN matrix and Neighbors list Assembly
  _self.IEN = np.zeros([_self.numElements,6], dtype = int)
  _self.FreedomDegree = len(_self.IEN[0,:])
  length = []

  for e in range(0, _self.numElements):
   v1 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][5]) - 1
   v2 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][6]) - 1
   v3 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][7]) - 1
   v4 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][8]) - 1
   v5 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][9]) - 1
   v6 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][10]) - 1
  
   _self.IEN[e] = [v1,v2,v3,v4,v5,v6]


  
   _self.neighborsNodes[v1].extend(_self.IEN[e])  
   _self.neighborsNodes[v2].extend(_self.IEN[e])  
   _self.neighborsNodes[v3].extend(_self.IEN[e])  
   _self.neighborsNodes[v4].extend(_self.IEN[e])  
   _self.neighborsNodes[v5].extend(_self.IEN[e])  
   _self.neighborsNodes[v6].extend(_self.IEN[e])  

   _self.neighborsNodes[v1] = list(set(_self.neighborsNodes[v1]))
   _self.neighborsNodes[v2] = list(set(_self.neighborsNodes[v2]))
   _self.neighborsNodes[v3] = list(set(_self.neighborsNodes[v3]))
   _self.neighborsNodes[v4] = list(set(_self.neighborsNodes[v4]))
   _self.neighborsNodes[v5] = list(set(_self.neighborsNodes[v5]))
   _self.neighborsNodes[v6] = list(set(_self.neighborsNodes[v6]))

 

   _self.neighborsNodesALE[v1].extend([v1,v2,v3,v4,v6])  
   _self.neighborsNodesALE[v2].extend([v1,v2,v3,v4,v5])  
   _self.neighborsNodesALE[v3].extend([v1,v2,v3,v5,v6])  
   _self.neighborsNodesALE[v4].extend([v1,v2])  
   _self.neighborsNodesALE[v5].extend([v2,v3])  
   _self.neighborsNodesALE[v6].extend([v3,v1])  
   
   _self.neighborsNodesALE[v1] = list(set(_self.neighborsNodesALE[v1]))
   _self.neighborsNodesALE[v2] = list(set(_self.neighborsNodesALE[v2]))
   _self.neighborsNodesALE[v3] = list(set(_self.neighborsNodesALE[v3]))
   _self.neighborsNodesALE[v4] = list(set(_self.neighborsNodesALE[v4]))
   _self.neighborsNodesALE[v5] = list(set(_self.neighborsNodesALE[v5]))
   _self.neighborsNodesALE[v6] = list(set(_self.neighborsNodesALE[v6]))

 
   
   _self.neighborsElements[v1].append(e)  
   _self.neighborsElements[v2].append(e)  
   _self.neighborsElements[v3].append(e)  
   _self.neighborsElements[v4].append(e)  
   _self.neighborsElements[v5].append(e)  
   _self.neighborsElements[v6].append(e)  




   x_a = _self.x[v1] - _self.x[v2]
   x_b = _self.x[v2] - _self.x[v3]
   x_c = _self.x[v3] - _self.x[v1]
   
   y_a = _self.y[v1] - _self.y[v2]
   y_b = _self.y[v2] - _self.y[v3]
   y_c = _self.y[v3] - _self.y[v1]
   
   length1 = np.sqrt(x_a**2 + y_a**2)
   length2 = np.sqrt(x_b**2 + y_b**2)
   length3 = np.sqrt(x_c**2 + y_c**2)

   length.append(length1)
   length.append(length2)
   length.append(length3)
   
  _self.minLengthMesh = min(length)




   
class Cubic2D:
 def __init__(_self, _dir, _file, _number_equations):
  _self.name = _dir + '/' + _file
  _self.neq = _number_equations
  _self.gmsh = []
  _self.neumann_lines = {}
  _self.neumann_edges = {}
  _self.dirichlet_lines = {}
  _self.dirichlet_pts = {}
  _self.neighbors_nodes = {}
  _self.neighbors_elements = {}
  _self.far_neighbors_nodes = {}
  _self.far_neighbors_elements = {}

  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.gmsh.append(row[:])

  
  _self.nphysical = int(_self.gmsh[4][0])
  _self.npoints = int(_self.gmsh[_self.nphysical+7][0])
  
  for i in range(1, _self.neq + 1):
   _self.neumann_lines[i] = []
   _self.neumann_edges[i] = []
   _self.dirichlet_lines[i] = []
   _self.dirichlet_pts[i] = []


  # Lines classification in neumann or dirichlet
  for i in range(0,_self.nphysical):
   for j in range(1,_self.neq + 1):
    if _self.gmsh[5+i][2] == '"neumann' + str(j):
     _self.neumann_lines[j].append(int(_self.gmsh[5+i][1]))

    elif _self.gmsh[5+i][2] == '"dirichlet' + str(j):
     _self.dirichlet_lines[j].append(int(_self.gmsh[5+i][1]))


  jj = _self.nphysical + _self.npoints 
  ii = 1
  element_type = int(_self.gmsh[(jj + 10 + ii)][1])

  # Assembly of the neumann edges and dirichlet points
  while element_type == 26:
    a_1 = int(_self.gmsh[(jj + 10 + ii)][3])
    a_2 = int(_self.gmsh[(jj + 10 + ii)][5])
    a_3 = int(_self.gmsh[(jj + 10 + ii)][6])
    a_4 = int(_self.gmsh[(jj + 10 + ii)][7])
    a_5 = int(_self.gmsh[(jj + 10 + ii)][8])
    a_6 = [a_1,a_2,a_3,a_4,a_5]
   
    for i in range(1, _self.neq + 1): 
     if a_1 in _self.neumann_lines[i]:
      _self.neumann_edges[i].append(a_6)
 
      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])
 
     elif a_1 in _self.dirichlet_lines[i]:
      _self.dirichlet_pts[i].append(a_6)

      ii = ii + 1
      element_type = int(_self.gmsh[(jj + 10 + ii)][1])

   
  for i in range(1, _self.neq + 1):
   _self.neumann_edges[i] = np.array(_self.neumann_edges[i]) 
   _self.dirichlet_pts[i] = np.array(_self.dirichlet_pts[i]) 


  _self.nelem = int(_self.gmsh[jj + 10][0]) - ii + 1

  for i in range(0, _self.npoints):
   _self.neighbors_nodes[i] = []
   _self.neighbors_elements[i] = []
   _self.far_neighbors_nodes[i] = []
   _self.far_neighbors_elements[i] = []


  _self.ii = ii
  _self.jj = jj


 def coord(_self):
  _self.x = np.zeros([_self.npoints,1], dtype = float)
  _self.y = np.zeros([_self.npoints,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.npoints):  
   _self.x[i] = _self.gmsh[_self.nphysical + 8 + i][1]
   _self.y[i] = _self.gmsh[_self.nphysical + 8 + i][2]
   _self.npts.append(i)


 def ien(_self):
  _self.IEN = np.zeros([_self.nelem,10], dtype = int)
  _self.GL = len(_self.IEN[0,:])
  _self.nodes_linear = [] 
  _self.nodes_quad = [] 
  length = [] 
  
  for e in range(0, _self.nelem):
   v1 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][5]) - 1
   v2 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][6]) - 1
   v3 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][7]) - 1
   v4 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][8]) - 1
   v5 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][9]) - 1
   v6 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][10]) - 1
   v7 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][11]) - 1
   v8 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][12]) - 1
   v9 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][13]) - 1
   v10 = int(_self.gmsh[(_self.jj + 10 + _self.ii + e)][14]) - 1
  
   _self.IEN[e] = [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10]
  
   _self.neighbors_nodes[v1].extend(_self.IEN[e])  
   _self.neighbors_nodes[v2].extend(_self.IEN[e])  
   _self.neighbors_nodes[v3].extend(_self.IEN[e])  
   _self.neighbors_nodes[v4].extend(_self.IEN[e])  
   _self.neighbors_nodes[v5].extend(_self.IEN[e])  
   _self.neighbors_nodes[v6].extend(_self.IEN[e])  
   _self.neighbors_nodes[v7].extend(_self.IEN[e])  
   _self.neighbors_nodes[v8].extend(_self.IEN[e])  
   _self.neighbors_nodes[v9].extend(_self.IEN[e])  
   _self.neighbors_nodes[v10].extend(_self.IEN[e])  
   
   _self.neighbors_nodes[v1] = list(set(_self.neighbors_nodes[v1]))
   _self.neighbors_nodes[v2] = list(set(_self.neighbors_nodes[v2]))
   _self.neighbors_nodes[v3] = list(set(_self.neighbors_nodes[v3]))
   _self.neighbors_nodes[v4] = list(set(_self.neighbors_nodes[v4]))
   _self.neighbors_nodes[v5] = list(set(_self.neighbors_nodes[v5]))
   _self.neighbors_nodes[v6] = list(set(_self.neighbors_nodes[v6]))
   _self.neighbors_nodes[v7] = list(set(_self.neighbors_nodes[v7]))
   _self.neighbors_nodes[v8] = list(set(_self.neighbors_nodes[v8]))
   _self.neighbors_nodes[v9] = list(set(_self.neighbors_nodes[v9]))
   _self.neighbors_nodes[v10] = list(set(_self.neighbors_nodes[v10]))
   
   _self.neighbors_elements[v1].append(e)  
   _self.neighbors_elements[v2].append(e)  
   _self.neighbors_elements[v3].append(e)  
   _self.neighbors_elements[v4].append(e)  
   _self.neighbors_elements[v5].append(e)  
   _self.neighbors_elements[v6].append(e)  
   _self.neighbors_elements[v7].append(e)  
   _self.neighbors_elements[v8].append(e)  
   _self.neighbors_elements[v9].append(e)  
   _self.neighbors_elements[v10].append(e)  
   
   _self.nodes_linear.append(v1)  
   _self.nodes_linear.append(v2)  
   _self.nodes_linear.append(v3)  

   _self.nodes_quad.append(v4)  
   _self.nodes_quad.append(v5)  
   _self.nodes_quad.append(v6)  
   _self.nodes_quad.append(v7)  
   _self.nodes_quad.append(v8)  
   _self.nodes_quad.append(v9)  
   _self.nodes_quad.append(v10)  

   x_a = _self.x[v1] - _self.x[v2]
   x_b = _self.x[v2] - _self.x[v3]
   x_c = _self.x[v3] - _self.x[v1]
   
   y_a = _self.y[v1] - _self.y[v2]
   y_b = _self.y[v2] - _self.y[v3]
   y_c = _self.y[v3] - _self.y[v1]
   
   length1 = np.sqrt(x_a**2 + y_a**2)
   length2 = np.sqrt(x_b**2 + y_b**2)
   length3 = np.sqrt(x_c**2 + y_c**2)

   length.append(length1)
   length.append(length2)
   length.append(length3)
   
  _self.length_min = min(length)

  _self.nodes_linear = list(set(_self.nodes_linear))  
  _self.nodes_quad = list(set(_self.nodes_quad))


  for i in range(0, _self.npoints):
   for j in _self.neighbors_nodes[i]:
    _self.far_neighbors_nodes[i].extend(_self.neighbors_nodes[j]) 
    _self.far_neighbors_elements[i].extend(_self.neighbors_elements[j]) 
 
   _self.far_neighbors_nodes[i] = list(set(_self.far_neighbors_nodes[i])\
                                     - set(_self.neighbors_nodes[i]))
   
   _self.far_neighbors_elements[i] = list(set(_self.far_neighbors_elements[i])\
                                        - set(_self.neighbors_elements[i]))

