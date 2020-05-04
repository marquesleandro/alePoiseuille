# ==========================================
# Code created by Leandro Marques at 12/2018
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used to compute boundary condition


import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg

# OBS: para vetor devemos unir eixo x e y no mesmo vetor, logo usar np.row_stack([dirichlet_pts[1],dirichlet_pts[2]])



class linearPoiseuille:

 # ------------------------------------------------------------------------------------------------------
 # Use:

 # # Applying vx condition
 # condition_xvelocity = bc_apply.Half_Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_xvelocity.neumann_condition(mesh.neumann_edges[1])
 # condition_xvelocity.dirichlet_condition(mesh.dirichlet_pts[1])
 # condition_xvelocity.gaussian_elimination(LHS_vx0,mesh.neighbors_nodes)
 # vorticity_ibc = condition_xvelocity.ibc

 # # Applying vy condition
 # condition_yvelocity = bc_apply.Half_Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_yvelocity.neumann_condition(mesh.neumann_edges[2])
 # condition_yvelocity.dirichlet_condition(mesh.dirichlet_pts[2])
 # condition_yvelocity.gaussian_elimination(LHS_vy0,mesh.neighbors_nodes)

 # # Applying psi condition
 # condition_streamfunction = bc_apply.Half_Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_streamfunction.streamfunction_condition(mesh.dirichlet_pts[3],LHS_psi0,mesh.neighbors_nodes)
 # ------------------------------------------------------------------------------------------------------


 def __init__(_self, _numPhysical, _numNodes, _x, _y):
  _self.numPhysical = _numPhysical
  _self.numNodes = _numNodes
  _self.x = _x
  _self.y = _y
  _self.maxVx = 3.0/2.0
  _self.L = 1.0
  _self.benchmark_problem = 'linear Poiseuille'


 def xVelocityCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Noslip 
   if line == 1 or line == 4:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Inflow
   elif line == 3:
    _self.aux1BC[v1] = (4.0*_self.maxVx)*(_self.y[v1]/_self.L**2)*(_self.L - _self.y[v1])
    _self.aux1BC[v2] = (4.0*_self.maxVx)*(_self.y[v2]/_self.L**2)*(_self.L - _self.y[v2])

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for vx
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 



 def yVelocityCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Noslip 
   if line == 1 or line == 4:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Inflow
   elif line == 3:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for vy
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 


 def streamFunctionCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.csr_matrix.copy(_LHS0) #used csr matrix because LHS = lil_matrix + lil_matrix
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Symmetric axis (Bottom Line)
   # psi_bottom can be any value. Because, important is psi_top - psi_bottom.
   # In this case, psi_bottom is zero
   if line == 4:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Noslip (Top Line)
   # Ref: Batchelor 1967 pag. 78 eq. 2.2.12
   # As psi_bottom is zero, so psi_top is:
   elif line == 1:
    _self.aux1BC[v1] = (_self.maxVx*(2.0/3.0))*(_self.L)
    _self.aux1BC[v2] = (_self.maxVx*(2.0/3.0))*(_self.L)

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for psi
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 

