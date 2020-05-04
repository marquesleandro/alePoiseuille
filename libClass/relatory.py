# ==========================================
# Code created by Leandro Marques at 02/2019
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used to create a simulation relatory


# ------------------------------------------------------------------
# Use:
# relatory.export(directory_name, sys.argv[0], scheme, mesh_name, equation_number, mesh.npoints, mesh.nelem, mesh.length_min, dt, nt, Re, Sc, import_mesh_time, assembly_time, bc_apply_time, solution_time)
# ------------------------------------------------------------------



import os
import datetime

def export(_save_path, _directory_name, _simulator_name, _benchmark_problem, _scheme, _mesh_name, _npoints, _nelem, _length_min, _dt, _nt, _Re, _Sc, _import_mesh_time, _assembly_time, _bc_apply_time, _solution_time, _polynomial_order, _gausspoints, observation):

 save_path = _save_path
 directory_name = _directory_name
 simulator_name = _simulator_name
 benchmark_problem = _benchmark_problem
 scheme = _scheme

 mesh_name = _mesh_name
 npoints = _npoints
 nelem = _nelem
 length_min = round(_length_min,7)
 dt = round(_dt,7)
 nt = _nt
 Re = _Re
 Sc = _Sc
 polynomial_order = _polynomial_order
 gausspoints = int(_gausspoints)

 import_mesh_time = round(_import_mesh_time,3)
 assembly_time = round(_assembly_time,3)
 bc_apply_time = round(_bc_apply_time,3)
 solution_time = round(_solution_time,3)
 
 path = save_path
 os.chdir(path)

 relatory_name = 'relatory_' + directory_name
 relatory = open(path + '/' + relatory_name + '.txt', 'w')



 today = datetime.datetime.now()
 relatory.write('                                 COPYRIGHT                                 \n')
 relatory.write(' ========================================================================== \n')
 relatory.write(' This relatory was created by Leandro Marques at ' + str(today) + '\n')
 relatory.write(' e-mail: marquesleandro67@gmail.com \n')
 relatory.write(' Gesar Search Group \n')
 relatory.write(' State University of the Rio de Janeiro \n')
 relatory.write(' ========================================================================== \n')
 relatory.write('\n')
 relatory.write('\n')

 relatory.write(' Simulator Name: ' + str(simulator_name) + '\n')
 relatory.write(' Benchmark Problem: ' + str(benchmark_problem) + '\n')
 relatory.write(' Scheme: ' + str(scheme) + '\n')
 relatory.write(' Element Type: ' + str(polynomial_order) + '\n')
 relatory.write(' Gaussian Quadrature (Gauss Points): ' + str(gausspoints) + '\n')
 relatory.write('\n')

 relatory.write(' ------------------------- \n')
 relatory.write(' PARAMETERS OF SIMULATION: \n')
 relatory.write(' ------------------------- \n')
 relatory.write(' Mesh: ' + str(mesh_name) + '\n')
 relatory.write(' Number of nodes: ' + str(npoints) + '\n')
 relatory.write(' Number of elements: ' + str(nelem) + '\n')
 relatory.write(' Smallest edge length: ' + str(length_min) + '\n')
 relatory.write(' Time step: ' + str(dt) + '\n')
 relatory.write(' Number of time iteration: ' + str(nt) + '\n')
 relatory.write(' Reynolds number: ' + str(Re) + '\n')
 relatory.write(' Schmidt number: ' + str(Sc) + '\n')
 relatory.write('\n')
 relatory.write('\n')


 relatory.write(' ------------ \n')
 relatory.write(' IMPORT MESH: \n')
 relatory.write(' ------------ \n')
 relatory.write(' time duration: ' + str(import_mesh_time) + ' seconds')
 relatory.write('\n')
 relatory.write('\n')


 relatory.write(' --------- \n')
 relatory.write(' ASSEMBLY: \n')
 relatory.write(' --------- \n')
 relatory.write(' time duration: ' + str(assembly_time) + ' seconds')
 relatory.write('\n')
 relatory.write('\n')


 relatory.write(' -------------------------------- \n')
 relatory.write(' INITIAL AND BOUNDARY CONDITIONS: \n')
 relatory.write(' -------------------------------- \n')
 relatory.write(' time duration: ' + str(bc_apply_time) + ' seconds')
 relatory.write('\n')
 relatory.write('\n')


 relatory.write(' ---------------------------- \n')
 relatory.write(' SOLVE THE LINEARS EQUATIONS: \n')
 relatory.write(' ---------------------------- \n')
 relatory.write(' time duration: ' + str(solution_time) + ' seconds')
 relatory.write('\n')
 relatory.write('\n')


 relatory.write(' ------------ \n')
 relatory.write(' OBSERVATION: \n')
 relatory.write(' ------------ \n')
 relatory.write(str(observation))
 relatory.write('\n')
 relatory.write('\n')



 
 relatory.close()
