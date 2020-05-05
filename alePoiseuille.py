# =======================
# Importing the libraries
# =======================

import os
initial_path = os.getcwd()

import sys
folderClass = './libClass'
sys.path.insert(0, folderClass)

from tqdm import tqdm
from time import time

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg

import searchMSH
import importMSH
import assembly
import benchmarkProblems
import importVTK
import ALE	
import semiLagrangian
import exportVTK
import relatory



print '''
               COPYRIGHT                    
 ======================================
 Simulator: %s
 created by Leandro Marques at 02/2019
 e-mail: marquesleandro67@gmail.com
 Gesar Search Group
 State University of the Rio de Janeiro
 ======================================
\n''' %sys.argv[0]



print ' ------'
print ' INPUT:'
print ' ------'

print ""


print ' ----------------------------------------------------------------------------'
print ' (1) - Simulation'
print ' (2) - Debug'
simulation_option = int(raw_input("\n enter simulation option above: "))
print' ----------------------------------------------------------------------------\n'


print ' ----------------------------------------------------------------------------'
print ' (1) - Linear Element'
print ' (2) - Mini Element'
print ' (3) - Quadratic Element'
print ' (4) - Cubic Element'
polynomial_option = int(raw_input("\n Enter polynomial degree option above: "))
print' ----------------------------------------------------------------------------\n'


print ' ----------------------------------------------------------------------------'
print ' 3 Gauss Points'
print ' 4 Gauss Points'
print ' 6 Gauss Points'
print ' 12 Gauss Points'
gausspoints = int(raw_input("\n Enter Gauss Points Number option above: "))
print' ----------------------------------------------------------------------------\n'


print ' ----------------------------------------------------------------------------'
print ' (1) - Taylor Galerkin Scheme'
print ' (2) - Semi Lagrangian Scheme'
scheme_option = int(raw_input("\n Enter simulation scheme option above: "))
print' ----------------------------------------------------------------------------\n'


print ' ----------------------------------------------------------------------------'
nt = int(raw_input(" Enter number of time interations (nt): "))
print' ----------------------------------------------------------------------------\n'


if simulation_option == 1:
 print ' ----------------------------------------------------------------------------'
 folderResults = raw_input(" Enter folder name to save simulations: ")
 print' ----------------------------------------------------------------------------\n'

 print ' ----------------------------------------------------------------------------'
 observation = raw_input(" Digit observation: ")
 print' ----------------------------------------------------------------------------\n'


elif simulation_option == 2:
 folderResults  = 'deletar'
 observation = 'debug'



print '\n ------------'
print ' IMPORT MESH:'
print ' ------------'

start_time = time()

# Linear and Mini Elements
if polynomial_option == 1 or polynomial_option == 2:
 mshFileName = 'linearPoiseuille.msh'
 #mshFileName = 'mesh1.msh'
 #mshFileName = 'mesh2.msh'
 #mshFileName = 'mesh3.msh'
 #mshFileName = 'mesh4.msh'
 #mshFileName = 'mesh5.msh'

 pathMSHFile = searchMSH.Find(mshFileName)
 if pathMSHFile == 'File not found':
  sys.exit()

 if polynomial_option == 1:
  mesh = importMSH.Linear2D(pathMSHFile, mshFileName)

  numNodes               = mesh.numNodes
  numElements            = mesh.numElements
  x                      = mesh.x
  y                      = mesh.y
  IEN                    = mesh.IEN
  boundaryEdges          = mesh.boundaryEdges
  boundaryNodes          = mesh.boundaryNodes
  neighborsNodes         = mesh.neighborsNodes
  neighborsNodesALE      = mesh.neighborsNodesALE
  neighborsElements      = mesh.neighborsElements
  minLengthMesh          = mesh.minLengthMesh
  FreedomDegree          = mesh.FreedomDegree
  numPhysical            = mesh.numPhysical 

  Re = 100.0
  Sc = 1.0
  CFL = 0.5
  #dt = float(CFL*minLengthMesh)
  dt = 0.1   #linear ok 


 
 elif polynomial_option == 2:
  mesh = importMSH.Mini2D(pathMSHFile, mshFileName)

  numNodes               = mesh.numNodes
  numElements            = mesh.numElements
  x                      = mesh.x
  y                      = mesh.y
  IEN                    = mesh.IEN
  boundaryEdges          = mesh.boundaryEdges
  boundaryNodes          = mesh.boundaryNodes
  neighborsNodes         = mesh.neighborsNodes
  neighborsNodesALE      = mesh.neighborsNodesALE
  neighborsElements      = mesh.neighborsElements
  minLengthMesh          = mesh.minLengthMesh
  FreedomDegree          = mesh.FreedomDegree
  numPhysical            = mesh.numPhysical 
  Re = 100.0
  Sc = 1.0
  CFL = 0.5
  dt = float(CFL*minLengthMesh)
  #dt = 0.1   #linear result ok 




# Quad Element
elif polynomial_option == 3:
 mshFileName = 'quadPoiseuille.msh'
 
 pathMSHFile = searchMSH.Find(mshFileName)
 if pathMSHFile == 'File not found':
  sys.exit()

 mesh = importMSH.Quad2D(pathMSHFile, mshFileName)

 numNodes               = mesh.numNodes
 numElements            = mesh.numElements
 x                      = mesh.x
 y                      = mesh.y
 IEN                    = mesh.IEN
 boundaryEdges          = mesh.boundaryEdges
 boundaryNodes          = mesh.boundaryNodes
 neighborsNodes         = mesh.neighborsNodes
 neighborsNodesALE      = mesh.neighborsNodesALE
 neighborsElements      = mesh.neighborsElements
 minLengthMesh          = mesh.minLengthMesh
 FreedomDegree          = mesh.FreedomDegree
 numPhysical            = mesh.numPhysical 

 Re = 100.0
 Sc = 1.0
 CFL = 0.5
 dt = float(CFL*minLengthMesh)
 #dt = 0.1  




# Cubic Element
elif polynomial_option == 4:
 mshFileName = 'cubicPoiseuille_cubic.msh'

 pathMSHFile = searchMSH.Find(mshFileName)
 if pathMSHFile == 'File not found':
  sys.exit()

 mesh = importMSH.Cubic2D(pathMSHFile, mshFileName, equation_number)
 mesh.coord()
 mesh.ien()



end_time = time()
import_mesh_time = end_time - start_time
print ' time duration: %.1f seconds \n' %import_mesh_time



print ' ---------'
print ' ASSEMBLY:'
print ' ---------'

start_time = time()
Kxx, Kxy, Kyx, Kyy, K, M, MLump, Gx, Gy, polynomial_order = assembly.Element2D(simulation_option, polynomial_option, FreedomDegree, numNodes, numElements, IEN, x, y, gausspoints)


end_time = time()
assembly_time = end_time - start_time
print ' time duration: %.1f seconds \n' %assembly_time




print ' --------------------------------'
print ' INITIAL AND BOUNDARY CONDITIONS:'
print ' --------------------------------'

start_time = time()


# ------------------------ Boundaries Conditions ----------------------------------

# Linear and Mini Elements
if polynomial_option == 1 or polynomial_option == 2:

 # Applying vx condition
 xVelocityLHS0 = sps.lil_matrix.copy(M)
 xVelocityBC = benchmarkProblems.linearPoiseuille(numPhysical,numNodes,x,y)
 xVelocityBC.xVelocityCondition(boundaryEdges,xVelocityLHS0,neighborsNodes)
 benchmark_problem = xVelocityBC.benchmark_problem

 # Applying vr condition
 yVelocityLHS0 = sps.lil_matrix.copy(M)
 yVelocityBC = benchmarkProblems.linearPoiseuille(numPhysical,numNodes,x,y)
 yVelocityBC.yVelocityCondition(boundaryEdges,yVelocityLHS0,neighborsNodes)
 
 # Applying psi condition
 streamFunctionLHS0 = sps.lil_matrix.copy(Kxx) + sps.lil_matrix.copy(Kyy)
 streamFunctionBC = benchmarkProblems.linearPoiseuille(numPhysical,numNodes,x,y)
 streamFunctionBC.streamFunctionCondition(boundaryEdges,streamFunctionLHS0,neighborsNodes)

 # Applying vorticity condition
 vorticityDirichletNodes = boundaryNodes



# Quad Element
elif polynomial_option == 3:

 # Applying vx condition
 xVelocityLHS0 = sps.lil_matrix.copy(M)
 xVelocityBC = benchmarkProblems.quadPoiseuille(numPhysical,numNodes,x,y)
 xVelocityBC.xVelocityCondition(boundaryEdges,xVelocityLHS0,neighborsNodes)
 benchmark_problem = xVelocityBC.benchmark_problem

 # Applying vr condition
 yVelocityLHS0 = sps.lil_matrix.copy(M)
 yVelocityBC = benchmarkProblems.quadPoiseuille(numPhysical,numNodes,x,y)
 yVelocityBC.yVelocityCondition(boundaryEdges,yVelocityLHS0,neighborsNodes)
 
 # Applying psi condition
 streamFunctionLHS0 = sps.lil_matrix.copy(Kxx) + sps.lil_matrix.copy(Kyy)
 streamFunctionBC = benchmarkProblems.quadPoiseuille(numPhysical,numNodes,x,y)
 streamFunctionBC.streamFunctionCondition(boundaryEdges,streamFunctionLHS0,neighborsNodes)

 # Applying vorticity condition
 vorticityDirichletNodes = boundaryNodes
# ---------------------------------------------------------------------------------



# -------------------------- Initial condition ------------------------------------
vx = np.copy(xVelocityBC.aux1BC)
vy = np.copy(yVelocityBC.aux1BC)
psi = np.copy(streamFunctionBC.aux1BC)
w = np.zeros([numNodes,1], dtype = float)
# ---------------------------------------------------------------------------------




#---------- Step 1 - Compute the vorticity and stream field --------------------
# -----Vorticity initial-----
vorticityRHS = sps.lil_matrix.dot(Gx,vy) - sps.lil_matrix.dot(Gy,vx)
vorticityLHS = sps.lil_matrix.copy(M)
w = scipy.sparse.linalg.cg(vorticityLHS,vorticityRHS,w, maxiter=1.0e+05, tol=1.0e-05)
w = w[0].reshape((len(w[0]),1))


# -----Streamline initial-----
streamFunctionRHS = sps.lil_matrix.dot(M,w)
streamFunctionRHS = np.multiply(streamFunctionRHS,streamFunctionBC.aux2BC)
streamFunctionRHS = streamFunctionRHS + streamFunctionBC.dirichletVector
psi = scipy.sparse.linalg.cg(streamFunctionBC.LHS,streamFunctionRHS,psi, maxiter=1.0e+05, tol=1.0e-05)
psi = psi[0].reshape((len(psi[0]),1))
#----------------------------------------------------------------------------------




# -------------------------- Import VTK File ------------------------------------
#numNodes, numElements, IEN, x, y, vx, vy, w, w, psi = importVTK.vtkFile("/home/marquesleandro/alePoiseuille/libClass/quad499.vtk", polynomial_option)
#----------------------------------------------------------------------------------




end_time = time()
bc_apply_time = end_time - start_time
print ' time duration: %.1f seconds \n' %bc_apply_time




print ' -----------------------------'
print ' PARAMETERS OF THE SIMULATION:'
print ' -----------------------------'

print ' Benchmark Problem: %s' %benchmark_problem
print ' Scheme: %s' %str(scheme_option)
print ' Element Type: %s' %str(polynomial_order)
print ' Gaussian Quadrature (Gauss Points): %s' %str(gausspoints)
print ' Mesh: %s' %mshFileName
print ' Number of nodes: %s' %numNodes
print ' Number of elements: %s' %numElements
print ' Smallest edge length: %f' %minLengthMesh
print ' Time step: %s' %dt
print ' Number of time iteration: %s' %nt
print ' Reynolds number: %s' %Re
print ' Schmidt number: %s' %Sc
print ""


print ' ----------------------------'
print ' SOLVE THE LINEARS EQUATIONS:'
print ' ---------------------------- \n'

print ' Saving simulation in %s \n' %folderResults



solution_start_time = time()
os.chdir(initial_path)



# ------------------------ Export VTK File ---------------------------------------
# Linear and Mini Elements
if polynomial_option == 1 or polynomial_option == 2:   
 save = exportVTK.Linear2D(x,y,IEN,numNodes,numElements,w,w,psi,vx,vy)
 save.create_dir(folderResults)
 save.saveVTK(folderResults + str(0))

# Quad Element
elif polynomial_option == 3:   
 save = exportVTK.Quad2D(x,y,IEN,numNodes,numElements,w,w,psi,vx,vy)
 save.create_dir(folderResults)
 save.saveVTK(folderResults + str(0))
# ---------------------------------------------------------------------------------



vorticityAux1BC = np.zeros([numNodes,1], dtype = float) 
x_old = np.zeros([numNodes,1], dtype = float)
y_old = np.zeros([numNodes,1], dtype = float)
vx_old = np.zeros([numNodes,1], dtype = float)
vy_old = np.zeros([numNodes,1], dtype = float)
end_type = 0
for t in tqdm(range(1, nt)):
 numIteration = t

 try:
  print ""
  print '''
                 COPYRIGHT                    
   ======================================
   Simulator: %s
   created by Leandro Marques at 04/2019
   e-mail: marquesleandro67@gmail.com
   Gesar Search Group
   State University of the Rio de Janeiro
   ======================================
  ''' %sys.argv[0]
 
 
 
  print ' -----------------------------'
  print ' PARAMETERS OF THE SIMULATION:'
  print ' -----------------------------'
 
  print ' Benchmark Problem: %s' %benchmark_problem
  print ' Scheme: %s' %str(scheme_option)
  print ' Element Type: %s' %str(polynomial_order)
  print ' Gaussian Quadrature (Gauss Points): %s' %str(gausspoints)
  print ' Mesh: %s' %mshFileName
  print ' Number of nodes: %s' %numNodes
  print ' Number of elements: %s' %numElements
  print ' Smallest edge length: %f' %minLengthMesh
  print ' Time step: %s' %dt
  print ' Number of time iteration: %s' %numIteration
  print ' Reynolds number: %s' %Re
  print ' Schmidt number: %s' %Sc
  print ""
 
 
 
  # ------------------------- ALE Scheme --------------------------------------------
  xmeshALE_dif = np.linalg.norm(x-x_old)
  ymeshALE_dif = np.linalg.norm(y-y_old)
  if not xmeshALE_dif < 5e-3 and not ymeshALE_dif < 5e-3:
   x_old = np.copy(x)
   y_old = np.copy(y)
  
   print ' ----'
   print ' ALE:'
   print ' ----'
  
  
   start_time = time()
  
  
   kLagrangian = 0.0
   kLaplace = 0.5
   kVelocity = 0.0
   
   vxLaplacianSmooth, vyLaplacianSmooth = ALE.Laplacian_smoothing(neighborsNodesALE, numNodes, x, y, dt)
   #vxLaplacianSmooth, vyLaplacianSmooth = ALE.Laplacian_smoothing_avg(neighborsNodesALE, numNodes, x, y, dt)
   vxVelocitySmooth,  vyVelocitySmooth  = ALE.Velocity_smoothing(neighborsNodesALE, numNodes, vx, vy)
 
   vxALE = kLagrangian*vx + kLaplace*vxLaplacianSmooth + kVelocity*vxVelocitySmooth
   vyALE = kLagrangian*vy + kLaplace*vyLaplacianSmooth + kVelocity*vyVelocitySmooth
  
  
   for i in boundaryNodes:
    node = i-1 
    vxALE[node] = 0.0
    vyALE[node] = 0.0
  
   x = x + vxALE*dt
   y = y + vyALE*dt
  
   vxSL = vx - vxALE
   vySL = vy - vyALE
  
   end_time = time()
   ALE_time_solver = end_time - start_time
   print ' time duration: %.1f seconds' %ALE_time_solver
   print ""
   # ---------------------------------------------------------------------------------
  
 
 
 
   # ------------------------- Assembly --------------------------------------------
   print ' ---------'
   print ' ASSEMBLY:'
   print ' ---------'
 
   Kxx, Kxy, Kyx, Kyy, K, M, MLump, Gx, Gy, polynomial_order = assembly.Element2D(simulation_option, polynomial_option, FreedomDegree, numNodes, numElements, IEN, x, y, gausspoints)
   print ""
   # --------------------------------------------------------------------------------
 
 
 
 
   # ------------------------ Boundaries Conditions ----------------------------------
   print ' --------------------------------'
   print ' INITIAL AND BOUNDARY CONDITIONS:'
   print ' --------------------------------'
   
   start_time = time()
  
   # Linear and Mini Elements
   if polynomial_option == 1 or polynomial_option == 2:
 
    # Applying vx condition
    xVelocityLHS0 = sps.lil_matrix.copy(M)
    xVelocityBC = benchmarkProblems.linearPoiseuille(numPhysical,numNodes,x,y)
    xVelocityBC.xVelocityCondition(boundaryEdges,xVelocityLHS0,neighborsNodes)
    benchmark_problem = xVelocityBC.benchmark_problem
   
    # Applying vr condition
    yVelocityLHS0 = sps.lil_matrix.copy(M)
    yVelocityBC = benchmarkProblems.linearPoiseuille(numPhysical,numNodes,x,y)
    yVelocityBC.yVelocityCondition(boundaryEdges,yVelocityLHS0,neighborsNodes)
    
    # Applying psi condition
    streamFunctionLHS0 = sps.lil_matrix.copy(Kxx) + sps.lil_matrix.copy(Kyy)
    streamFunctionBC = benchmarkProblems.linearPoiseuille(numPhysical,numNodes,x,y)
    streamFunctionBC.streamFunctionCondition(boundaryEdges,streamFunctionLHS0,neighborsNodes)
   
    # Applying vorticity condition
    vorticityDirichletNodes = boundaryNodes
   
   
  
   # Quad Element
   elif polynomial_option == 3:
 
    # Applying vx condition
    xVelocityLHS0 = sps.lil_matrix.copy(M)
    xVelocityBC = benchmarkProblems.quadPoiseuille(numPhysical,numNodes,x,y)
    xVelocityBC.xVelocityCondition(boundaryEdges,xVelocityLHS0,neighborsNodes)
    benchmark_problem = xVelocityBC.benchmark_problem
   
    # Applying vr condition
    yVelocityLHS0 = sps.lil_matrix.copy(M)
    yVelocityBC = benchmarkProblems.quadPoiseuille(numPhysical,numNodes,x,y)
    yVelocityBC.yVelocityCondition(boundaryEdges,yVelocityLHS0,neighborsNodes)
    
    # Applying psi condition
    streamFunctionLHS0 = sps.lil_matrix.copy(Kxx) + sps.lil_matrix.copy(Kyy)
    streamFunctionBC = benchmarkProblems.quadPoiseuille(numPhysical,numNodes,x,y)
    streamFunctionBC.streamFunctionCondition(boundaryEdges,streamFunctionLHS0,neighborsNodes)
   
    # Applying vorticity condition
    vorticityDirichletNodes = boundaryNodes
   
   
   end_time = time()
   bc_apply_time_solver = end_time - start_time
   print ' time duration: %.1f seconds' %bc_apply_time_solver
   print ""
   # ---------------------------------------------------------------------------------
   
   
 
  # ------------------------ SOLVE LINEAR EQUATIONS ----------------------------------
  print ' ----------------------------'
  print ' SOLVE THE LINEARS EQUATIONS:'
  print ' ----------------------------'
  print ""
  print ' Saving simulation in %s' %folderResults
  print ""
 
  start_solver_time = time()
 
 
 
 
  #---------- Step 2 - Compute the boundary conditions for vorticity --------------
  vorticityRHS = sps.lil_matrix.dot(Gx,vy) - sps.lil_matrix.dot(Gy,vx)
  vorticityLHS = sps.lil_matrix.copy(M)
  vorticityAux1BC = scipy.sparse.linalg.cg(vorticityLHS,vorticityRHS,vorticityAux1BC, maxiter=1.0e+05, tol=1.0e-05)
  vorticityAux1BC = vorticityAux1BC[0].reshape((len(vorticityAux1BC[0]),1))
 
  # Gaussian elimination
  vorticityDirichletVector = np.zeros([numNodes,1], dtype = float)
  vorticityNeumannVector = np.zeros([numNodes,1], dtype = float)
  vorticityAux2BC = np.ones([numNodes,1], dtype = float)
 
  vorticityLHS = (np.copy(M)/dt) + (1.0/Re)*np.copy(Kxx) + (1.0/Re)*np.copy(Kyy)
  for mm in vorticityDirichletNodes:
   for nn in neighborsNodes[mm]:
    vorticityDirichletVector[nn] -= float(vorticityLHS[nn,mm]*vorticityAux1BC[mm])
    vorticityLHS[nn,mm] = 0.0
    vorticityLHS[mm,nn] = 0.0
    
   vorticityLHS[mm,mm] = 1.0
   vorticityDirichletVector[mm] = vorticityAux1BC[mm]
   vorticityAux2BC[mm] = 0.0
  #----------------------------------------------------------------------------------
 
 
 
  #---------- Step 3 - Solve the vorticity transport equation ----------------------
  # Taylor Galerkin Scheme
  if scheme_option == 1:
   scheme_name = 'Taylor Galerkin'
   A = np.copy(M)/dt 
   vorticityRHS = sps.lil_matrix.dot(A,w) - np.multiply(vx,sps.lil_matrix.dot(Gx,w))\
         - np.multiply(vy,sps.lil_matrix.dot(Gy,w))\
         - (dt/2.0)*np.multiply(vx,(np.multiply(vx,sps.lil_matrix.dot(Kxx,w)) + np.multiply(vy,sps.lil_matrix.dot(Kyx,w))))\
         - (dt/2.0)*np.multiply(vy,(np.multiply(vx,sps.lil_matrix.dot(Kxy,w)) + np.multiply(vy,sps.lil_matrix.dot(Kyy,w))))
   vorticityRHS = np.multiply(vorticityRHS,vorticityAux2BC)
   vorticityRHS = vorticityRHS + vorticityDirichletVector
   w = scipy.sparse.linalg.cg(vorticityLHS,vorticityRHS,w, maxiter=1.0e+05, tol=1.0e-05)
   w = w[0].reshape((len(w[0]),1))
 
 
 
  # Semi-Lagrangian Scheme
  elif scheme_option == 2:
 
   # Linear Element   
   if polynomial_option == 1:
    scheme_name = 'Semi Lagrangian Linear'
 
    w_d = semiLagrangian.Linear2D(numNodes, neighborsElements, IEN, x, y, vxSL, vySL, dt, w)
 
    A = np.copy(M)/dt
    vorticityRHS = sps.lil_matrix.dot(A,w_d)
 
    vorticityRHS = vorticityRHS + (1.0/Re)*vorticityNeumannVector
    vorticityRHS = np.multiply(vorticityRHS,vorticityAux2BC)
    vorticityRHS = vorticityRHS + vorticityDirichletVector
 
    w = scipy.sparse.linalg.cg(vorticityLHS,vorticityRHS,w, maxiter=1.0e+05, tol=1.0e-05)
    w = w[0].reshape((len(w[0]),1))
 
 
 
   # Mini Element   
   elif polynomial_option == 2:
    scheme_name = 'Semi Lagrangian Mini'
 
    w_d = semiLagrangian.Mini2D(numNodes, neighborsElements, IEN, z, r, vz, vr, dt, w)
 
    A = np.copy(Mr)/dt
    vorticityRHS = sps.lil_matrix.dot(A,w_d)
 
    vorticityRHS = vorticityRHS + (1.0/Re)*vorticityNeumannVector
    vorticityRHS = np.multiply(vorticityRHS,vorticityAux2BC)
    vorticityRHS = vorticityRHS + vorticityDirichletVector
 
    w = scipy.sparse.linalg.cg(vorticityLHS,vorticityRHS,w, maxiter=1.0e+05, tol=1.0e-05)
    w = w[0].reshape((len(w[0]),1))
 
 
 
   # Quad Element   
   elif polynomial_option == 3:
    scheme_name = 'Semi Lagrangian Quad'
 
    w_d = semiLagrangian.Quad2D(numNodes, neighborsElements, IEN, x, y, vxSL, vySL, dt, w)
 
    A = np.copy(M)/dt
    vorticityRHS = sps.lil_matrix.dot(A,w_d)
 
    vorticityRHS = vorticityRHS + (1.0/Re)*vorticityNeumannVector
    vorticityRHS = np.multiply(vorticityRHS,vorticityAux2BC)
    vorticityRHS = vorticityRHS + vorticityDirichletVector
 
    w = scipy.sparse.linalg.cg(vorticityLHS,vorticityRHS,w, maxiter=1.0e+05, tol=1.0e-05)
    w = w[0].reshape((len(w[0]),1)) 
  #----------------------------------------------------------------------------------
 
 
 
  #---------- Step 4 - Solve the streamline equation --------------------------------
  # Solve Streamline
  # psi condition
  streamFunctionRHS = sps.lil_matrix.dot(M,w)
  streamFunctionRHS = np.multiply(streamFunctionRHS,streamFunctionBC.aux2BC)
  streamFunctionRHS = streamFunctionRHS + streamFunctionBC.dirichletVector
  psi = scipy.sparse.linalg.cg(streamFunctionBC.LHS,streamFunctionRHS,psi, maxiter=1.0e+05, tol=1.0e-05)
  psi = psi[0].reshape((len(psi[0]),1))
  #----------------------------------------------------------------------------------
 
 
 
  #---------- Step 5 - Compute the velocity field -----------------------------------
  # Velocity vx
  vx_old = np.copy(vx)
  xVelocityRHS = sps.lil_matrix.dot(Gy,psi)
  xVelocityRHS = np.multiply(xVelocityRHS,xVelocityBC.aux2BC)
  xVelocityRHS = xVelocityRHS + xVelocityBC.dirichletVector
  vx = scipy.sparse.linalg.cg(xVelocityBC.LHS,xVelocityRHS,vx, maxiter=1.0e+05, tol=1.0e-05)
  vx = vx[0].reshape((len(vx[0]),1))
  
  # Velocity vy
  vy_old = np.copy(vy)
  yVelocityRHS = -sps.lil_matrix.dot(Gx,psi)
  yVelocityRHS = np.multiply(yVelocityRHS,yVelocityBC.aux2BC)
  yVelocityRHS = yVelocityRHS + yVelocityBC.dirichletVector
  vy = scipy.sparse.linalg.cg(yVelocityBC.LHS,yVelocityRHS,vy, maxiter=1.0e+05, tol=1.0e-05)
  vy = vy[0].reshape((len(vy[0]),1))
  #----------------------------------------------------------------------------------
 
  end_solver_time = time()
  solver_time = end_solver_time - start_solver_time
  print ' time duration: %.1f seconds' %solver_time
  print ""
  #----------------------------------------------------------------------------------
  
 
 
 
 
  # ------------------------ Export VTK File ---------------------------------------
  # Linear and Mini Elements
  if polynomial_option == 1 or polynomial_option == 2:   
   save = exportVTK.Linear2D(x,y,IEN,numNodes,numElements,w,w,psi,vx,vy)
   save.create_dir(folderResults)
   save.saveVTK(folderResults + str(t))
 
  # Quad Element
  elif polynomial_option == 3:   
   save = exportVTK.Quad2D(x,y,IEN,numNodes,numElements,w,w,psi,vx,vy)
   save.create_dir(folderResults)
   save.saveVTK(folderResults + str(t))
  # ---------------------------------------------------------------------------------
 
 
 
 
  # ------------------------ CHECK STEADY STATE ----------------------------------
  if np.all(vx == vx_old) and np.all(vy == vy_old):
   end_type = 1
   break
  # ---------------------------------------------------------------------------------
 
  # ------------------------ CHECK CONVERGENCE RESULT ----------------------------------
  if np.linalg.norm(vx) > 10e2 or np.linalg.norm(vy) > 10e2:
   end_type = 2
   break
  # ---------------------------------------------------------------------------------
  

 except KeyboardInterrupt:
  end_type = 3
  break 




end_time = time()
solution_time = end_time - solution_start_time
print ' time duration: %.1f seconds \n' %solution_time



print ' ----------------'
print ' SAVING RELATORY:'
print ' ----------------'
print ""

if end_type == 0:
 print ' END SIMULATION. NOT STEADY STATE'
 print ' Relatory saved in %s' %folderResults
 print ""

elif end_type == 1:
 print ' END SIMULATION. STEADY STATE'
 print ' Relatory saved in %s' %folderResults
 print ""

elif end_type == 2:
 print ' END SIMULATION. ERROR CONVERGENCE RESULT'
 print ' Relatory saved in %s' %folderResults
 print ""

elif end_type == 3:
 print ' END SIMULATION. FORCED INTERRUPTION'
 print ' Relatory saved in %s' %folderResults
 print ""




# -------------------------------- Export Relatory ---------------------------------------
relatory.export(save.path, folderResults, sys.argv[0], benchmark_problem, scheme_name, mshFileName, numNodes, numElements, minLengthMesh, dt, numIteration, Re, Sc, import_mesh_time, assembly_time, bc_apply_time, solution_time, polynomial_order, gausspoints, observation)
# ----------------------------------------------------------------------------------------



