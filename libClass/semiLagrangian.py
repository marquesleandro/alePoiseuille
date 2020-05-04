# ==========================================
# Code created by Leandro Marques at 12/2018
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used to use semi-lagrangian scheme

# ------------------------------------------------------------------------------------
# Use:
# scalar_d = semi_lagrangian.Linear2D(
# mesh.npoints, mesh.neighbors_elements, mesh.IEN, mesh.x, mesh.y, vx, vy, dt, scalar_n)
# ------------------------------------------------------------------------------------

import numpy as np


# 1D Semi-Lagrangian using npoints x nelem to find departure node
# correct method
def Linear1D_v2(_npoints, _nelem, _IEN, _xn, _vx, _dt, _scalar):
 xd = _xn - _vx*_dt

 scalar = np.zeros([_npoints,1], dtype = float) 

 for i in range(0,_npoints):
  x = float(xd[i])

  breaking = 0
  length = []

  for e in range(0,_nelem):
   v1 = _IEN[e][0]
   v2 = _IEN[e][1]

   x1 = float(_xn[v1])
   x2 = float(_xn[v2])

   len1 = abs(x2 - x)
   len2 = abs(x1 - x)
   lent = abs(x1 - x2)

   L1 = len1/lent
   L2 = len2/lent

   alpha = [L1,L2]
   alpha = np.array(alpha)

   if np.all(alpha >= 0.0) and np.all(alpha <= 1.0):
    N1 = L1
    N2 = L2

    scalar1 = _scalar[v1]
    scalar2 = _scalar[v2]

    scalar[i] = N1*scalar1 + N2*scalar2
    breaking = 1
    break

   else:
    x_a = x1 - x
    x_b = x2 - x
  
    length1 = np.sqrt(x_a**2)
    length2 = np.sqrt(x_b**2)

    a_1 = [v1,length1]
    a_2 = [v2,length2]
 
    length.append(a_1)
    length.append(a_2)
   
  if breaking == 0:
   length_min = min(length, key=lambda k:k[1])
   node = length_min[0]
   scalar[i] = _scalar[node]
  
 return scalar  

# fast find method
# correct method
def Linear1D(_npoints, _neighbors_elements, _IEN, _xn, _vx, _dt, _scalar):
 xd = _xn - _vx*_dt
 
 scalar = np.zeros([_npoints,1], dtype = float) 

 for i in range(0,_npoints):
  x = float(xd[i])

  node = i
  length = []
  breaking = 0
  while breaking == 0:
   for e in _neighbors_elements[node]:
    v1 = _IEN[e][0]
    v2 = _IEN[e][1]

    x1 = float(_xn[v1])
    x2 = float(_xn[v2])

    len1 = abs(x2 - x)
    len2 = abs(x1 - x)
    lent = abs(x1 - x2)

    L1 = len1/lent
    L2 = len2/lent

    alpha = [L1,L2]
    alpha = np.array(alpha)
 
    if np.all(alpha >= 0.0) and np.all(alpha <= 1.0):
     N1 = L1
     N2 = L2
 
     scalar1 = _scalar[v1]
     scalar2 = _scalar[v2]

     scalar[i] = N1*scalar1 + N2*scalar2
     breaking = 1
     break

    else:
     x_a = x1 - x
     x_b = x2 - x
  
     length1 = np.sqrt(x_a**2)
     length2 = np.sqrt(x_b**2)

     a_1 = [v1,length1]
     a_2 = [v2,length2]
 
     length.append(a_1)
     length.append(a_2)
   
     breaking = 0


   # first neighbor is element found 
   if breaking == 1:
     break
 
 
   # coordinate not found
   else:
    length_min = min(length, key=lambda k:k[1])
    node1 = node
    node = length_min[0]

    # outside domain
    if node == node1 and breaking == 0:
     scalar[i] = _scalar[node]
     
     breaking = 1
     break


 return scalar  






# 1D Semi-Lagrangian using npoints x nelem to find departure node
# correct method
def Quad1D_v2(_npoints, _nelem, _IEN, _xn, _vx, _dt, _scalar):
 xd = _xn - _vx*_dt

 scalar = np.zeros([_npoints,1], dtype = float) 

 for i in range(0,_npoints):
  x = float(xd[i])

  breaking = 0
  length = []

  for e in range(0,_nelem):
   v1 = _IEN[e][0]
   v2 = _IEN[e][1]
   v3 = _IEN[e][2]

   x1 = float(_xn[v1])
   x2 = float(_xn[v2])
   x3 = float(_xn[v3])

   len1 = abs(x2 - x)
   len2 = abs(x1 - x)
   lent = abs(x1 - x2)

   L1 = len1/lent
   L2 = len2/lent

   alpha = [L1,L2]
   alpha = np.array(alpha)

   if np.all(alpha >= 0.0) and np.all(alpha <= 1.0):
    
    # Lewis pag 48 Eq. 3.33
    # x3 is quadratic node
    N1 = ((x - x3)/(x1 - x3))*((x - x2)/(x1 - x2))
    N2 = ((x - x1)/(x2 - x1))*((x - x3)/(x2 - x3))
    N3 = ((x - x1)/(x3 - x1))*((x - x2)/(x3 - x2))

    scalar1 = _scalar[v1]
    scalar2 = _scalar[v2]
    scalar3 = _scalar[v3]

    scalar[i] = N1*scalar1 + N2*scalar2 + N3*scalar3

    # Interpolation limits
    scalar_limits = [scalar1,scalar2,scalar3]
    if scalar[i] < min(scalar_limits):
     scalar[i] = min(scalar_limits)

    elif scalar[i] > max(scalar_limits):
     scalar[i] = max(scalar_limits)


    breaking = 1
    break


   else:
    x_a = x1 - x
    x_b = x2 - x
  
    length1 = np.sqrt(x_a**2)
    length2 = np.sqrt(x_b**2)

    a_1 = [v1,length1]
    a_2 = [v2,length2]
 
    length.append(a_1)
    length.append(a_2)
   
  if breaking == 0:
   length_min = min(length, key=lambda k:k[1])
   node = length_min[0]
   scalar[i] = _scalar[node]
  
 return scalar  





# fast find method with quad points find
# correct method
def Quad1D(_npoints, _neighbors_elements, _IEN, _xn, _vx, _dt, _scalar):
 xd = _xn - _vx*_dt

 scalar = np.zeros([_npoints,1], dtype = float) 

 for i in range(0,_npoints):
  x = float(xd[i])

  node = i
  length = []
  breaking = 0
  while breaking == 0:
   for e in _neighbors_elements[node]:
    v1 = _IEN[e][0]
    v2 = _IEN[e][1]
    v3 = _IEN[e][2]

    x1 = float(_xn[v1])
    x2 = float(_xn[v2])
    x3 = float(_xn[v3])

    len1 = abs(x2 - x)
    len2 = abs(x1 - x)
    lent = abs(x1 - x2)

    L1 = len1/lent
    L2 = len2/lent

    alpha = [L1,L2]
    alpha = np.array(alpha)

    if np.all(alpha >= 0.0) and np.all(alpha <= 1.0):
    
     # Lewis pag 48 Eq. 3.33
     # x3 is quadratic node
     N1 = ((x - x3)/(x1 - x3))*((x - x2)/(x1 - x2))
     N2 = ((x - x1)/(x2 - x1))*((x - x3)/(x2 - x3))
     N3 = ((x - x1)/(x3 - x1))*((x - x2)/(x3 - x2))

     scalar1 = _scalar[v1]
     scalar2 = _scalar[v2]
     scalar3 = _scalar[v3]

     scalar[i] = N1*scalar1 + N2*scalar2 + N3*scalar3

     # Interpolation limits
     scalar_limits = [scalar1,scalar2,scalar3]
     if scalar[i] < min(scalar_limits):
      scalar[i] = min(scalar_limits)

     elif scalar[i] > max(scalar_limits):
      scalar[i] = max(scalar_limits)

     breaking = 1
     break

    else:
     x_a = x1 - x
     x_b = x2 - x
  
     length1 = np.sqrt(x_a**2)
     length2 = np.sqrt(x_b**2)

     a_1 = [v1,length1]
     a_2 = [v2,length2]
 
     length.append(a_1)
     length.append(a_2)
     
     breaking = 0
   
   # first neighbor is element found 
   if breaking == 1:
     break
 
 
   # coordinate not found
   else:
    length_min = min(length, key=lambda k:k[1])
    node1 = node
    node = length_min[0]

    # outside domain
    if node == node1 and breaking == 0:
     scalar[i] = _scalar[node]
     
     breaking = 1
     break


 return scalar  






# 2D Semi-Lagrangian using npoints x nelem to find departure node
def Linear2D_v2(_npoints, _nelem, _IEN, _xn, _yn, _vx, _vy, _dt, _scalar):
 xd = _xn - _vx*_dt
 yd = _yn - _vy*_dt
 
 scalar = np.zeros([_npoints,1], dtype = float) 
 
 for i in range(0,_npoints):
  x = float(xd[i])
  y = float(yd[i])
  
  breaking = 0
  length = []

  for e in range(0,_nelem):
   v1 = _IEN[e][0]
   v2 = _IEN[e][1]
   v3 = _IEN[e][2]

   x1 = float(_xn[v1])
   x2 = float(_xn[v2])
   x3 = float(_xn[v3])

   y1 = float(_yn[v1])
   y2 = float(_yn[v2])
   y3 = float(_yn[v3])

   A = np.array([[x1,x2,x3],
                 [y1,y2,y3],
                 [1.0,1.0,1.0]])

   b = np.array([x,y,1.0])
 
   alpha = np.linalg.solve(A,b)

   
   if np.all(alpha >= 0.0) and np.all(alpha <= 1.0):
    A1 = 0.5*np.linalg.det(np.array([[1, x, y],
                                     [1, x2, y2],
                                     [1, x3, y3]]))
 
    A2 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                     [1, x, y],
                                     [1, x3, y3]]))
 
    A3 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                     [1, x2, y2],
                                     [1, x, y]]))
 
    At = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                     [1, x2, y2],
                                     [1, x3, y3]]))
   
    L1 = A1/At
    L2 = A2/At
    L3 = A3/At
     
    N1 = L1
    N2 = L2
    N3 = L3

    scalar1 = _scalar[v1]
    scalar2 = _scalar[v2]
    scalar3 = _scalar[v3]

    scalar[i] = N1*scalar1 + N2*scalar2 + N3*scalar3
    breaking = 1
    break

   else:
    x_a = x1 - x
    x_b = x2 - x
    x_c = x3 - x
   
    y_a = y1 - y
    y_b = y2 - y
    y_c = y3 - y
  
    length1 = np.sqrt(x_a**2 + y_a**2)
    length2 = np.sqrt(x_b**2 + y_b**2)
    length3 = np.sqrt(x_c**2 + y_c**2)

    a_1 = [v1,length1]
    a_2 = [v2,length2]
    a_3 = [v3,length3]
 
    length.append(a_1)
    length.append(a_2)
    length.append(a_3)
   
  if breaking == 0:
   length_min = min(length, key=lambda k:k[1])
   node = length_min[0]
   scalar[i] = _scalar[node]
     
 return scalar  





# 2D Semi-Lagrangian using npoints x neighbors_elements to find departure node
# semilagrangian test ok
def Linear2D(_npoints, _neighbors_elements, _IEN, _xn, _yn, _vx, _vy, _dt, _scalar):
 xd = _xn - _vx*_dt
 yd = _yn - _vy*_dt
 
 scalar = np.zeros([_npoints,1], dtype = float) 
 
 for i in range(0,_npoints):
  x = float(xd[i])
  y = float(yd[i])

  node = i
  length = []
  breaking = 0

  while breaking == 0:
   for e in _neighbors_elements[node]:
    v1 = _IEN[e][0]
    v2 = _IEN[e][1]
    v3 = _IEN[e][2]

    x1 = float(_xn[v1])
    x2 = float(_xn[v2])
    x3 = float(_xn[v3])

    y1 = float(_yn[v1])
    y2 = float(_yn[v2])
    y3 = float(_yn[v3])
  
    A = np.array([[x1,x2,x3],
                  [y1,y2,y3],
                  [1.0,1.0,1.0]])

    b = np.array([x,y,1.0])
 
    alpha = np.linalg.solve(A,b)

 
    if np.all(alpha >= 0.0) and np.all(alpha <= 1.0):
 
     A1 = 0.5*np.linalg.det(np.array([[1, x, y],
                                      [1, x2, y2],
                                      [1, x3, y3]]))
 
     A2 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x, y],
                                      [1, x3, y3]]))
 
     A3 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x2, y2],
                                      [1, x, y]]))
 
     At = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x2, y2],
                                      [1, x3, y3]]))
   
     L1 = A1/At
     L2 = A2/At
     L3 = A3/At
     
     N1 = L1
     N2 = L2
     N3 = L3

     scalar1 = _scalar[v1]
     scalar2 = _scalar[v2]
     scalar3 = _scalar[v3]

     scalar[i] = N1*scalar1 + N2*scalar2 + N3*scalar3

     breaking = 1
     break


    else:
     x_a = x1 - x
     x_b = x2 - x
     x_c = x3 - x
   
     y_a = y1 - y
     y_b = y2 - y
     y_c = y3 - y
  
     length1 = np.sqrt(x_a**2 + y_a**2)
     length2 = np.sqrt(x_b**2 + y_b**2)
     length3 = np.sqrt(x_c**2 + y_c**2)

     a_1 = [v1,length1]
     a_2 = [v2,length2]
     a_3 = [v3,length3]
 
     length.append(a_1)
     length.append(a_2)
     length.append(a_3)
   
     breaking = 0


   # first neighbor is element found 
   if breaking == 1:
     break
 
 
   # coordinate not found
   else:
    length_min = min(length, key=lambda k:k[1])
    node1 = node
    node = length_min[0]

    # outside domain
    if node == node1 and breaking == 0:
     scalar[i] = _scalar[node]
     
     breaking = 1
     break


 return scalar




# 2D Semi-Lagrangian using npoints x neighbors_elements to find departure node
# neighbors distance
# semilagrangian test OK 
def Linear2D_v3(_npoints, _neighbors_nodes, _neighbors_elements, _IEN, _xn, _yn, _vx, _vy, _dt, _scalar):
 xd = _xn - _vx*_dt
 yd = _yn - _vy*_dt
 
 scalar = np.zeros([_npoints,1], dtype = float) 
 
 for i in range(0,_npoints):
  x = float(xd[i])
  y = float(yd[i])

  node = i
  length = []
  barycentric = []
  element = []
  breaking = 0

  #while breaking == 0: #until find
  for k in range(0,10): #search limited
   for j in _neighbors_nodes[node]:
    xx = float(_xn[j])
    yy = float(_yn[j])
 
    x_a = xx - x
    y_a = yy - y
  
    length1 = np.sqrt(x_a**2 + y_a**2)

    a_1 = [j,length1]
 
    length.append(a_1)
   
   length_min = min(length, key=lambda k:k[1])
   node1 = node
   node = length_min[0]

   # node more short
   if node == node1:
    for e in _neighbors_elements[node]:
     v1 = _IEN[e][0]
     v2 = _IEN[e][1]
     v3 = _IEN[e][2]

     x1 = float(_xn[v1])
     x2 = float(_xn[v2])
     x3 = float(_xn[v3])

     y1 = float(_yn[v1])
     y2 = float(_yn[v2])
     y3 = float(_yn[v3])
  
     A = np.array([[x1,x2,x3],
                   [y1,y2,y3],
                   [1.0,1.0,1.0]])

     b = np.array([x,y,1.0])
 
     alpha = np.linalg.solve(A,b)

     barycentric.append(alpha)
     element.append(e)

    barycentric = np.array(barycentric)
    idx = (barycentric >= 0.0) & (barycentric <= 1.0)
    aa = np.where((idx == True).all(axis=1))

    # inside element
    try:
     aa = int(aa[0])
     e = element[aa]

     v1 = _IEN[e][0]
     v2 = _IEN[e][1]
     v3 = _IEN[e][2]

     x1 = float(_xn[v1])
     x2 = float(_xn[v2])
     x3 = float(_xn[v3])

     y1 = float(_yn[v1])
     y2 = float(_yn[v2])
     y3 = float(_yn[v3])
  

     A1 = 0.5*np.linalg.det(np.array([[1, x, y],
                                      [1, x2, y2],
                                      [1, x3, y3]]))
 
     A2 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x, y],
                                      [1, x3, y3]]))
 
     A3 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x2, y2],
                                      [1, x, y]]))
 
     At = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x2, y2],
                                      [1, x3, y3]]))
   
     L1 = A1/At
     L2 = A2/At
     L3 = A3/At
     
     N1 = L1
     N2 = L2
     N3 = L3

     scalar1 = _scalar[v1]
     scalar2 = _scalar[v2]
     scalar3 = _scalar[v3]

     scalar[i] = N1*scalar1 + N2*scalar2 + N3*scalar3
     
     breaking = 1
     break

    # outside domain
    except TypeError:
     scalar[i] = _scalar[node]
    
     breaking = 1
     break

   # node far yet
   else:
    continue

   break


 return scalar

# 2D Semi-Lagrangian using npoints x neighbors_elements to find departure node
def Mini2D(_npoints, _neighbors_elements, _IEN, _xn, _yn, _vx, _vy, _dt, _scalar):
 xd = _xn - _vx*_dt
 yd = _yn - _vy*_dt
 
 scalar = np.zeros([_npoints,1], dtype = float) 
 
 for i in range(0,_npoints):
  x = float(xd[i])
  y = float(yd[i])

  node = i
  length = []
  breaking = 0

  while breaking == 0:
   for e in _neighbors_elements[node]:
    v1 = _IEN[e][0]
    v2 = _IEN[e][1]
    v3 = _IEN[e][2]

    x1 = float(_xn[v1])
    x2 = float(_xn[v2])
    x3 = float(_xn[v3])

    y1 = float(_yn[v1])
    y2 = float(_yn[v2])
    y3 = float(_yn[v3])
  
    A = np.array([[x1,x2,x3],
                  [y1,y2,y3],
                  [1.0,1.0,1.0]])

    b = np.array([x,y,1.0])
 
    alpha = np.linalg.solve(A,b)

 
    if np.all(alpha >= 0.0) and np.all(alpha <= 1.0):
     v4 = _IEN[e][3]

     A1 = 0.5*np.linalg.det(np.array([[1, x, y],
                                      [1, x2, y2],
                                      [1, x3, y3]]))
 
     A2 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x, y],
                                      [1, x3, y3]]))
 
     A3 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x2, y2],
                                      [1, x, y]]))
 
     At = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x2, y2],
                                      [1, x3, y3]]))
   
     L1 = A1/At
     L2 = A2/At
     L3 = A3/At

     N1 = L1 - 9.0*L1*L2*L3    
     N2 = L2 - 9.0*L1*L2*L3    
     N3 = L3 - 9.0*L1*L2*L3
     N4 = 27.0*L1*L2*L3      

     scalar1 = _scalar[v1]
     scalar2 = _scalar[v2]
     scalar3 = _scalar[v3]
     scalar4 = _scalar[v4]

     scalar[i] = N1*scalar1 + N2*scalar2 + N3*scalar3 + N4*scalar4

     # Interpolation limits
     scalar_limits = [scalar1,scalar2,scalar3,scalar4]
     if scalar[i] < min(scalar_limits):
      scalar[i] = min(scalar_limits)

     elif scalar[i] > max(scalar_limits):
      scalar[i] = max(scalar_limits)

     breaking = 1
     break


    else:
     x_a = x1 - x
     x_b = x2 - x
     x_c = x3 - x
   
     y_a = y1 - y
     y_b = y2 - y
     y_c = y3 - y
  
     length1 = np.sqrt(x_a**2 + y_a**2)
     length2 = np.sqrt(x_b**2 + y_b**2)
     length3 = np.sqrt(x_c**2 + y_c**2)

     a_1 = [v1,length1]
     a_2 = [v2,length2]
     a_3 = [v3,length3]
 
     length.append(a_1)
     length.append(a_2)
     length.append(a_3)
   
     breaking = 0


   # first neighbor is element found 
   if breaking == 1:
     break
 
 
   # coordinate not found
   else:
    length_min = min(length, key=lambda k:k[1])
    node1 = node
    node = length_min[0]

    # outside domain
    if node == node1 and breaking == 0:
     scalar[i] = _scalar[node]
     
     breaking = 1
     break


 return scalar







# 2D Semi-Lagrangian using npoints x neighbors_elements to find departure node
# semilagrangian test ok
def Quad2D(_npoints, _neighbors_elements, _IEN, _xn, _yn, _vx, _vy, _dt, _scalar):
 xd = _xn - _vx*_dt
 yd = _yn - _vy*_dt
 
 scalar = np.zeros([_npoints,1], dtype = float) 
 
 for i in range(0,_npoints):
  x = float(xd[i])
  y = float(yd[i])

  node = i
  length = []
  breaking = 0

  while breaking == 0:
   for e in _neighbors_elements[node]:
    v1 = _IEN[e][0]
    v2 = _IEN[e][1]
    v3 = _IEN[e][2]

    x1 = float(_xn[v1])
    x2 = float(_xn[v2])
    x3 = float(_xn[v3])

    y1 = float(_yn[v1])
    y2 = float(_yn[v2])
    y3 = float(_yn[v3])
  
    A = np.array([[x1,x2,x3],
                  [y1,y2,y3],
                  [1.0,1.0,1.0]])

    b = np.array([x,y,1.0])
 
    alpha = np.linalg.solve(A,b)

 
    if np.all(alpha >= 0.0) and np.all(alpha <= 1.0):
     
     v4 = _IEN[e][3]
     v5 = _IEN[e][4]
     v6 = _IEN[e][5]

     A1 = 0.5*np.linalg.det(np.array([[1, x, y],
                                      [1, x2, y2],
                                      [1, x3, y3]]))
 
     A2 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x, y],
                                      [1, x3, y3]]))
 
     A3 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x2, y2],
                                      [1, x, y]]))
 
     At = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x2, y2],
                                      [1, x3, y3]]))
   
     L1 = A1/At
     L2 = A2/At
     L3 = A3/At
     
     N1 = L1*(2.0*L1 - 1.0)
     N2 = L2*(2.0*L2 - 1.0)
     N3 = L3*(2.0*L3 - 1.0)
     N4 = 4.0*L1*L2
     N5 = 4.0*L2*L3
     N6 = 4.0*L3*L1

     scalar1 = _scalar[v1]
     scalar2 = _scalar[v2]
     scalar3 = _scalar[v3]
     scalar4 = _scalar[v4]
     scalar5 = _scalar[v5]
     scalar6 = _scalar[v6]

     scalar[i] = N1*scalar1 + N2*scalar2 + N3*scalar3 + N4*scalar4 + N5*scalar5 + N6*scalar6

     # Interpolation limits
     scalar_limits = [scalar1,scalar2,scalar3,scalar4,scalar5,scalar6]
     if scalar[i] < min(scalar_limits):
      scalar[i] = min(scalar_limits)

     elif scalar[i] > max(scalar_limits):
      scalar[i] = max(scalar_limits)

     breaking = 1
     break


    else:
     x_a = x1 - x
     x_b = x2 - x
     x_c = x3 - x
   
     y_a = y1 - y
     y_b = y2 - y
     y_c = y3 - y
  
     length1 = np.sqrt(x_a**2 + y_a**2)
     length2 = np.sqrt(x_b**2 + y_b**2)
     length3 = np.sqrt(x_c**2 + y_c**2)

     a_1 = [v1,length1]
     a_2 = [v2,length2]
     a_3 = [v3,length3]
 
     length.append(a_1)
     length.append(a_2)
     length.append(a_3)
   
     breaking = 0


   # first neighbor is element found 
   if breaking == 1:
     break
 
 
   # coordinate not found
   else:
    length_min = min(length, key=lambda k:k[1])
    node1 = node
    node = length_min[0]

    # outside domain
    if node == node1 and breaking == 0:
     scalar[i] = _scalar[node]
     
     breaking = 1
     break


 return scalar




# 2D Semi-Lagrangian using npoints x nelem to find departure node
def Quad2D_v2(_npoints, _nelem, _IEN, _xn, _yn, _vx, _vy, _dt, _scalar):
 xd = _xn - _vx*_dt
 yd = _yn - _vy*_dt
 
 scalar = np.zeros([_npoints,1], dtype = float) 
 
 for i in range(0,_npoints):
  x = float(xd[i])
  y = float(yd[i])
  
  breaking = 0
  length = []

  for e in range(0,_nelem):
   v1 = _IEN[e][0]
   v2 = _IEN[e][1]
   v3 = _IEN[e][2]
   v4 = _IEN[e][3]
   v5 = _IEN[e][4]
   v6 = _IEN[e][5]

   x1 = float(_xn[v1])
   x2 = float(_xn[v2])
   x3 = float(_xn[v3])
   x4 = float(_xn[v4])
   x5 = float(_xn[v5])
   x6 = float(_xn[v6])

   y1 = float(_yn[v1])
   y2 = float(_yn[v2])
   y3 = float(_yn[v3])
   y4 = float(_yn[v4])
   y5 = float(_yn[v5])
   y6 = float(_yn[v6])
  
   A = np.array([[x1,x2,x3],
                 [y1,y2,y3],
                 [1.0,1.0,1.0]])

   b = np.array([x,y,1.0])
 
   alpha = np.linalg.solve(A,b)

 
   if np.all(alpha >= 0.0) and np.all(alpha <= 1.0):
 
    A1 = 0.5*np.linalg.det(np.array([[1, x, y],
                                     [1, x2, y2],
                                     [1, x3, y3]]))
 
    A2 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                     [1, x, y],
                                     [1, x3, y3]]))
 
    A3 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                     [1, x2, y2],
                                     [1, x, y]]))
 
    At = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                     [1, x2, y2],
                                     [1, x3, y3]]))
   
    L1 = A1/At
    L2 = A2/At
    L3 = A3/At
     
    N1 = L1*(2.0*L1 - 1.0)
    N2 = L2*(2.0*L2 - 1.0)
    N3 = L3*(2.0*L3 - 1.0)
    N4 = 4.0*L1*L2
    N5 = 4.0*L2*L3
    N6 = 4.0*L3*L1

    scalar1 = _scalar[v1]
    scalar2 = _scalar[v2]
    scalar3 = _scalar[v3]
    scalar4 = _scalar[v4]
    scalar5 = _scalar[v5]
    scalar6 = _scalar[v6]

    scalar[i] = N1*scalar1 + N2*scalar2 + N3*scalar3 + N4*scalar4 + N5*scalar5 + N6*scalar6

    breaking = 1
    break


   else:
    x_a = x1 - x
    x_b = x2 - x
    x_c = x3 - x
   
    y_a = y1 - y
    y_b = y2 - y
    y_c = y3 - y
  
    length1 = np.sqrt(x_a**2 + y_a**2)
    length2 = np.sqrt(x_b**2 + y_b**2)
    length3 = np.sqrt(x_c**2 + y_c**2)

    a_1 = [v1,length1]
    a_2 = [v2,length2]
    a_3 = [v3,length3]
 
    length.append(a_1)
    length.append(a_2)
    length.append(a_3)
   
    breaking = 0

  if breaking == 0:
   length_min = min(length, key=lambda k:k[1])
   node = length_min[0]
   scalar[i] = _scalar[node]
     
 return scalar  




# 2D Semi-Lagrangian using npoints x neighbors_elements to find departure node
# neighbors distance
# semilagrangian test ok
def Quad2D_v3(_npoints, _neighbors_nodes, _neighbors_elements, _IEN, _xn, _yn, _vx, _vy, _dt, _scalar):
 xd = _xn - _vx*_dt
 yd = _yn - _vy*_dt
 
 scalar = np.zeros([_npoints,1], dtype = float) 
 
 for i in range(0,_npoints):
  x = float(xd[i])
  y = float(yd[i])

  node = i
  length = []
  barycentric = []
  element = []
  breaking = 0

  #while breaking == 0: #until find
  for k in range(0,10): #search limited
   for j in _neighbors_nodes[node]:
    xx = float(_xn[j])
    yy = float(_yn[j])
 
    x_a = xx - x
    y_a = yy - y
  
    length1 = np.sqrt(x_a**2 + y_a**2)

    a_1 = [j,length1]
 
    length.append(a_1)
   
   length_min = min(length, key=lambda k:k[1])
   node1 = node
   node = length_min[0]

   # node more short
   if node == node1:
    for e in _neighbors_elements[node]:
     v1 = _IEN[e][0]
     v2 = _IEN[e][1]
     v3 = _IEN[e][2]

     x1 = float(_xn[v1])
     x2 = float(_xn[v2])
     x3 = float(_xn[v3])

     y1 = float(_yn[v1])
     y2 = float(_yn[v2])
     y3 = float(_yn[v3])
  
     A = np.array([[x1,x2,x3],
                   [y1,y2,y3],
                   [1.0,1.0,1.0]])

     b = np.array([x,y,1.0])
 
     alpha = np.linalg.solve(A,b)

     barycentric.append(alpha)
     element.append(e)

    barycentric = np.array(barycentric)
    idx = (barycentric >= 0.0) & (barycentric <= 1.0)
    aa = np.where((idx == True).all(axis=1))

    # inside element
    try:
     aa = int(aa[0])
     e = element[aa]

     v1 = _IEN[e][0]
     v2 = _IEN[e][1]
     v3 = _IEN[e][2]
     v4 = _IEN[e][3]
     v5 = _IEN[e][4]
     v6 = _IEN[e][5]

     x1 = float(_xn[v1])
     x2 = float(_xn[v2])
     x3 = float(_xn[v3])
     x4 = float(_xn[v4])
     x5 = float(_xn[v5])
     x6 = float(_xn[v6])

     y1 = float(_yn[v1])
     y2 = float(_yn[v2])
     y3 = float(_yn[v3])
     y4 = float(_yn[v4])
     y5 = float(_yn[v5])
     y6 = float(_yn[v6])
  

     A1 = 0.5*np.linalg.det(np.array([[1, x, y],
                                      [1, x2, y2],
                                      [1, x3, y3]]))
 
     A2 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x, y],
                                      [1, x3, y3]]))
 
     A3 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x2, y2],
                                      [1, x, y]]))
 
     At = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x2, y2],
                                      [1, x3, y3]]))
   
     L1 = A1/At
     L2 = A2/At
     L3 = A3/At
     
     N1 = L1*(2.0*L1 - 1.0)
     N2 = L2*(2.0*L2 - 1.0)
     N3 = L3*(2.0*L3 - 1.0)
     N4 = 4.0*L1*L2
     N5 = 4.0*L2*L3
     N6 = 4.0*L3*L1

     scalar1 = _scalar[v1]
     scalar2 = _scalar[v2]
     scalar3 = _scalar[v3]
     scalar4 = _scalar[v4]
     scalar5 = _scalar[v5]
     scalar6 = _scalar[v6]

     scalar[i] = N1*scalar1 + N2*scalar2 + N3*scalar3 + N4*scalar4 + N5*scalar5 + N6*scalar6

     breaking = 1
     break

    # outside domain
    except TypeError:
     scalar[i] = _scalar[node]
    
     breaking = 1
     break

   # node far yet
   else:
    continue

   break


 return scalar



# 2D Semi-Lagrangian using npoints x neighbors_elements to find departure node
def Cubic2D(_npoints, _neighbors_elements, _IEN, _xn, _yn, _vx, _vy, _dt, _scalar):
 xd = _xn - _vx*_dt
 yd = _yn - _vy*_dt
 
 scalar = np.zeros([_npoints,1], dtype = float) 
 
 for i in range(0,_npoints):
  x = float(xd[i])
  y = float(yd[i])

  node = i
  length = []
  breaking = 0

  while breaking == 0:
   for e in _neighbors_elements[node]:
    v1 = _IEN[e][0]
    v2 = _IEN[e][1]
    v3 = _IEN[e][2]

    x1 = float(_xn[v1])
    x2 = float(_xn[v2])
    x3 = float(_xn[v3])

    y1 = float(_yn[v1])
    y2 = float(_yn[v2])
    y3 = float(_yn[v3])
  
    A = np.array([[x1,x2,x3],
                  [y1,y2,y3],
                  [1.0,1.0,1.0]])

    b = np.array([x,y,1.0])
 
    alpha = np.linalg.solve(A,b)

 
    if np.all(alpha >= 0.0) and np.all(alpha <= 1.0):
     
     v4 = _IEN[e][3]
     v5 = _IEN[e][4]
     v6 = _IEN[e][5]
     v7 = _IEN[e][6]
     v8 = _IEN[e][7]
     v9 = _IEN[e][8]
     v10 = _IEN[e][9]

     A1 = 0.5*np.linalg.det(np.array([[1, x, y],
                                      [1, x2, y2],
                                      [1, x3, y3]]))
 
     A2 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x, y],
                                      [1, x3, y3]]))
 
     A3 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x2, y2],
                                      [1, x, y]]))
 
     At = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                      [1, x2, y2],
                                      [1, x3, y3]]))
   
     L1 = A1/At
     L2 = A2/At
     L3 = A3/At
     
     N1 = 0.5*L1*(3.0*L1 - 1.0)*(3.0*L1 - 2.0)   #N1 = 1/2*L1*(3*L1-1)*(3*L1-2)
     N2 = 0.5*L2*(3.0*L2 - 1.0)*(3.0*L2 - 2.0)   #N2 = 1/2*L2*(3*L2-1)*(3*L2-2)
     N3 = 0.5*L3*(3.0*L3 - 1.0)*(3.0*L3 - 2.0)   #N3 = 1/2*L3*(3*L3-1)*(3*L3-2)
     N4 = 4.5*L1*L2*(3.0*L1 - 1.0)               #N4 = 9/2*L1*L2(3*L1-1)
     N5 = 4.5*L1*L2*(3.0*L2 - 1.0)               #N5 = 9/2*L1*L2(3*L2-1)
     N6 = 4.5*L2*L3*(3.0*L2 - 1.0)               #N6 = 9/2*L2*L3(3*L2-1)
     N7 = 4.5*L2*L3*(3.0*L3 - 1.0)               #N7 = 9/2*L2*L3(3*L3-1)
     N8 = 4.5*L3*L1*(3.0*L3 - 1.0)               #N8 = 9/2*L3*L1(3*L3-1)
     N9 = 4.5*L3*L1*(3.0*L1 - 1.0)               #N9 = 9/2*L3*L1(3*L1-1)
     N10 = 27.0*L1*L2*L3                         #N10 = 27*L1*L2*L3

     scalar1 = _scalar[v1]
     scalar2 = _scalar[v2]
     scalar3 = _scalar[v3]
     scalar4 = _scalar[v4]
     scalar5 = _scalar[v5]
     scalar6 = _scalar[v6]
     scalar7 = _scalar[v7]
     scalar8 = _scalar[v8]
     scalar9 = _scalar[v9]
     scalar10 = _scalar[v10]

     scalar[i] = N1*scalar1 + N2*scalar2 + N3*scalar3 + N4*scalar4 + N5*scalar5 + N6*scalar6 + N7*scalar7 + N8*scalar8 + N9*scalar9 + N10*scalar10

     # Interpolation limits
     scalar_limits = [scalar1,scalar2,scalar3,scalar4,scalar5,scalar6,scalar7,scalar8,scalar9,scalar10]
     if scalar[i] < min(scalar_limits):
      scalar[i] = min(scalar_limits)

     elif scalar[i] > max(scalar_limits):
      scalar[i] = max(scalar_limits)
 

     breaking = 1
     break


    else:
     x_a = x1 - x
     x_b = x2 - x
     x_c = x3 - x
   
     y_a = y1 - y
     y_b = y2 - y
     y_c = y3 - y
  
     length1 = np.sqrt(x_a**2 + y_a**2)
     length2 = np.sqrt(x_b**2 + y_b**2)
     length3 = np.sqrt(x_c**2 + y_c**2)

     a_1 = [v1,length1]
     a_2 = [v2,length2]
     a_3 = [v3,length3]
 
     length.append(a_1)
     length.append(a_2)
     length.append(a_3)
   
     breaking = 0


   # first neighbor is element found 
   if breaking == 1:
     break
 
 
   # coordinate not found
   else:
    length_min = min(length, key=lambda k:k[1])
    node1 = node
    node = length_min[0]

    # outside domain
    if node == node1 and breaking == 0:
     scalar[i] = _scalar[node]
     
     breaking = 1
     break


 return scalar




