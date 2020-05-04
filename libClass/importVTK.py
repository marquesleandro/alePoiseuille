# ==========================================
# Code created by Leandro Marques at 03/2020
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used to import .vtk file


# Converting .msh in a python list

import numpy as np

def vtkFile(_file, _polynomial_option): 

 vtkList = [] 
 with open(_file) as vtkFile:
   for line in vtkFile:
    row = line.split()
    vtkList.append(row[:])

 for i in range(0,len(vtkList)):
  for j in range(0,len(vtkList[i])):
   if vtkList[i][j] == "POINTS":
    numNodes = int(vtkList[i][j+1])

    x = np.zeros([numNodes,1], dtype = float)
    y = np.zeros([numNodes,1], dtype = float)
    for k in range(0,numNodes):
     x[k] = float(vtkList[i+k+1][0])
     y[k] = float(vtkList[i+k+1][1])
    continue  

   if vtkList[i][j] == "CELLS":
    numElements = int(vtkList[i][j+1])

    # Linear Element
    if _polynomial_option == 1:
     IEN = np.zeros([numElements,3], dtype = int)
     for e in range(0,numElements):
      IEN[e][0] = int(vtkList[i+e+1][1])
      IEN[e][1] = int(vtkList[i+e+1][2])
      IEN[e][2] = int(vtkList[i+e+1][3])
     continue 

    # Quad Element
    elif _polynomial_option == 3:
     IEN = np.zeros([numElements,6], dtype = int)
     for e in range(0,numElements):
      IEN[e][0] = int(vtkList[i+e+1][1])
      IEN[e][1] = int(vtkList[i+e+1][2])
      IEN[e][2] = int(vtkList[i+e+1][3])
      IEN[e][3] = int(vtkList[i+e+1][4])
      IEN[e][4] = int(vtkList[i+e+1][5])
      IEN[e][5] = int(vtkList[i+e+1][6])
     continue 


   if vtkList[i][j] == "VECTORS":
    vx = np.zeros([numNodes,1], dtype = float)
    vy = np.zeros([numNodes,1], dtype = float)
    for k in range(0,numNodes):
     vx[k] = float(vtkList[i+k+1][0])
     vy[k] = float(vtkList[i+k+1][1])
    continue  

   if vtkList[i][j] == "scalar1":
    scalar1 = np.zeros([numNodes,1], dtype = float)
    for k in range(0,numNodes):
     scalar1[k] = float(vtkList[i+k+2][0])
    continue  

   if vtkList[i][j] == "scalar2":
    scalar2 = np.zeros([numNodes,1], dtype = float)
    for k in range(0,numNodes):
     scalar2[k] = float(vtkList[i+k+2][0])
    continue  

   if vtkList[i][j] == "scalar3":
    scalar3 = np.zeros([numNodes,1], dtype = float)
    for k in range(0,numNodes):
     scalar3[k] = float(vtkList[i+k+2][0])
    continue  

 return numNodes, numElements, IEN, x, y, vx, vy, scalar1, scalar2, scalar3
