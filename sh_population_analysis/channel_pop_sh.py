import math
import sys
import numpy as np
import random
import time
import os
import string
import scipy
from scipy import signal
#from memory_profiler import profile
##############################################
##############################################
##############################################
# ls -1 movie_* | tr '\n'  ' '


def input_check():
    global molecule,natoms,results_file   # same number and molecule for all movies and geoms
    if len(sys.argv) < 2:
      print("Error: not enought parameters.\nUsage: python ",sys.argv[0]," th/tm/molecule movie.xyz movie2.xyz....")
      sys.exit(1)
    ntrajs = int(sys.argv[2])
    movies  = []
    lines   = []
    geoms   = []
    results_file = sys.argv[1]   # file used for saving final data
    
    ##### MODIFIE HERE FOR SPECIFIC MOLECULE: #####

    mov_files = []
    for ntr in range(1,ntrajs+1):
        mov_file = "TRAJ.{}/movie.xyz".format(ntr)
        cwd = os.getcwd()
        mov_path = os.path.join(cwd,mov_file)
        if os.path.isfile(mov_path):
            movies.append(mov_path)
            gc = geoms_check(movies[-1],lines_per_mol) #no need to specify order
            #print(gc)
            lines.append(gc[0])
            geoms.append(gc[1])
            print(mov_file,' OK, Nlines: ', lines[-1], 'Ngeoms: ', geoms[-1])
        else:
            print(mov_path ,' File NOT EXISTS.')  
    return movies,geoms
#end input check

def geoms_check(mov,lines_per_mol):   # fast number_of_lines_ reader exploiting limited buffer                 
    lines = 0
    geoms = 0
    buf_size = 1024 * 1024
    with open(mov,'r') as f:
      read_f = f.read
      buf = read_f(buf_size)
      while buf:
         lines += buf.count('\n')  
         #\n) is left at the end of the string, and is only omitted on the last line of the file
         buf = read_f(buf_size)
      f.close()
    if not (lines % lines_per_mol): geoms = lines / lines_per_mol  # 0 FALSE
    else:
       print('Nlines is not divisible by l_p_m: ',lines, lines_per_mol, mov)
       print('Check if there is empty line at the end of the file \n')
       sys.exit(1)
    return lines,geoms

# MAIN ROUTINE TO GO THROUGH EACH MOVIE - READ XYZ, CALCULATE DISTANCE, ANALYZE GEOMETRY
#@profile
def process_movies(movies,geoms):
    """
     Expecting .xyz file 
     first line = natoms
     second line = comment + time/timestep information, might require change in timestep assingment split index []
    """
    analyzed_geoms = np.array([[0, 0]],dtype=np.float64)                                          # main array with time and reaction channel for each geometry
    for m,mov in enumerate(movies):                                              # iterate over movies
         print("Processing ",m+1,"movie: ",mov)
         with open(mov,'r') as f:
                                       
          for g in range(1,int(geoms[m])+1):                                     # iterate over geoms in each mov file, first index is inclusive, last exclusive!
              #natoms = int(f.readline())
              atoms = f.readline()                                               # atoms
              timestep = f.readline().split()[2]                                 # comment + time 
              if g == 1:
                 xyz = np.zeros(shape=(natoms,3),dtype=np.float64)
                 if os.path.isfile(os.path.join(os.getcwd(),'dist_mat.dat')): os.remove('dist_mat.dat')     
              #print('geometry: ',g)

              for at in range(0,natoms):                                         # iterate over atoms in each geometry
                  line = f.readline().split()
                  xyz[at]=[float(line[1]),float(line[2]),float(line[3])]
              #print(time,"\n",xyz)
              
              dist_mat = distance_matrix(xyz)                                    #cal dist matrix
              
              ##### MODIFIE HERE FOR SPECIFIC MOLECULE: #####
              channel  = analyze_novec(dist_mat)        # analyze geometry
                                            
              analyzed_geoms = np.append(analyzed_geoms, [[int(timestep),int(channel)]], axis = 0)  # save analyzed data for statistics
          f.close()
          
    return(analyzed_geoms)
     
# DISTANCE MATRIX
#@profile
def distance_matrix(xyz):    
# all combinations of pairs: list(itertools.combinations(range(natoms),2)) - yet still need to loop over two-indices to call dist func brute force number of combinations len(<-)
# with open('dist_mat.dat','a') as file_dist_save:  # save dist_mat in file for check if needed
# np.savetxt(file_dist_save, dist_mat, newline='\n', fmt='%.8e',footer =" ")
    dist_mat = np.zeros(shape=(natoms-1,natoms),dtype=np.float64) # create empty dist matrix - matrix is not stored for future
    for k in range(0,natoms):
          for l in range(k+1,natoms):
                v1, v2 = np.array(xyz[k],dtype=np.float64), np.array(xyz[l],dtype=np.float64)
                dist = [(a - b)**2 for a, b in zip(v1, v2)]
                dist_mat[k][l] = math.sqrt(sum(dist))
                # print(v1,v2,l,k) # combination check

    return dist_mat

###############################################
#ConstantS - CRITERIA FOR GEOMETRY ANALYSIS 
C_CF3_diss_dist  =  4.30  # distance of CF3 group from center
C_F_diss_dist = 4.30
C_CN_diss_dist  =  4.30
###############################################

# GEOMETRY ANALYSIS
#@profile 
def analyze_novec(dist_mat):

  channel = 0
  cf3_diss = 0 # how many cf3 groups diss 1,2
  cf_bonds = 0 # bonds in cf3 group = 3
  cf2_diss = 0
  f_diss = 0
  cn_diss = 0
  # 1    C CF3       C-CF3  0,1
  if dist_mat[0][1] > C_CF3_diss_dist:  
      for cf3_atoms in range(8,11):
          if dist_mat[1][cf3_atoms] < C_F_diss_dist:  
            cf_bonds += 1
      if cf_bonds == 3:
         cf3_diss += 1
      if cf_bonds == 2:
         cf2_diss += 1 
  cf_bonds = 0 
  if dist_mat[0][4] > C_CF3_diss_dist:  
      for cf3_atoms in range(5,8):
          if dist_mat[4][cf3_atoms] < C_F_diss_dist:  
            cf_bonds += 1
      if cf_bonds == 3:
         cf3_diss += 1 
      if cf_bonds == 2:
         cf2_diss += 1

  if dist_mat[0][2] > C_CN_diss_dist:
      cn_diss += 1  

  for cf3_atoms in range(5,11):
      if (dist_mat[0][cf3_atoms] > C_F_diss_dist 
          and dist_mat[1][cf3_atoms] > C_F_diss_dist 
          and dist_mat[2][cf3_atoms] > C_F_diss_dist 
          and dist_mat[4][cf3_atoms] > C_F_diss_dist 
          and dist_mat[cf3_atoms][11] > C_F_diss_dist):  
                f_diss += 1 
                            
  if   cf3_diss == 1 and f_diss == 0: channel = 1
  elif cf3_diss == 2 and f_diss == 0: channel = 2
  elif cf3_diss == 0 and f_diss == 1: channel = 3
  elif cf3_diss == 0 and f_diss == 2: channel = 4
  elif cf2_diss == 1 and f_diss == 1: channel = 5 
  elif cf3_diss == 0 and f_diss == 0 and cn_diss == 1: channel = 6 
  """
  atom order:
  0    C CENTRAL   
  1    C CF3       C-CF3  0,1
  2    C CN        C-CN   0,2
  3    F CF        C-F    0,3
  4    C CF3       C-CF3  0,4
  5-7  F in CF3 with C4
  8-10 F in CF3 with C1
  11   N in CN 2,1 C-N    2,11       
 
  Channels:
  0 nothing happened
  1 ONE CF3 DISS
  2 two CF3 DISS
  3 one F diss
  4 two F diss
  5 ONE CF2 and F DISS 
  """

  return channel
     
           
def channel_statistics(analyze_geoms):
    """
    MODIFIE PARAMETERS FOR EACH TYPE OF MOLECULE (n_channels)
    nstep, timestep depends on simulation number of steps (e.g. nsteps in input.in)
    """
    channel_pop = np.zeros(shape=(n_steps,n_channels),dtype=np.float64)   # 2D array, 0 column time, rest {1,n_channel} are channels
    procentual = 1         # 0 - 1 or 0-100
    
    print("Collectiong data. Total number of geoms: ",len(analyze_geoms)-1)
    for rec in range(0,len(analyze_geoms)-1):              # first row is 0,0 entry from array init
      channel = int(analyze_geoms[rec][1])
      step    = int(analyze_geoms[rec][0])-1 #timerow from 1, not from 0
      channel_pop[step][channel] = channel_pop[step][channel] + 1  
    np.savetxt("analyze_geoms", analyze_geoms, fmt='%7.3f',delimiter=' ')
    nsteps = n_steps
    nstates = n_channels  
    timerow = np.arange(1,nsteps+1,1,dtype=np.float64).reshape(nsteps,1) * 0.024 * timestep
    stde = np.zeros(shape=(nsteps,nstates),dtype=np.float64)
    populations = np.copy(channel_pop)
    row_norms = populations.sum(axis=1) # row sum is axis=1 opposite to numpe row priority
    row_norms_mat = np.resize(row_norms,(nstates,nsteps)).transpose() # to broadcast the following division, need the same dimensions
    norm_pop = np.divide(populations,row_norms_mat)  # Normalization populations
    
    #Calc std binomidal distribution 95% reliability
    for step in range(0,nsteps):     
        for st in range(nstates):
            if (step == 100) or (not step % (200)):
                stde[step][st]=math.sqrt(norm_pop[step][st]*(1-norm_pop[step][st])/row_norms.reshape(nsteps,1)[step])*1.96  # binobidal normal dist    
   
    alive_mat=np.concatenate((timerow,row_norms.reshape(nsteps,1)),axis=1)	
    with open ("alive.dat", "w") as ad:
        np.savetxt(ad, alive_mat, fmt='%7.3f',delimiter=' ') 
   
    for st in range(0,nstates):  # 0 is time axis
        # take only individual state due to xmgrace import XYDY: time pop_st std_st
        outfile="channel_{}.dat".format(st+1)
        norm_pop_st = norm_pop[:,st]
        norm_pop_st_smooth = scipy.signal.savgol_filter(norm_pop_st, 33, 7).reshape(nsteps,1)
        st_stde = stde[:,st].reshape(nsteps,1) 
        st_pop_stde = np.concatenate((timerow, norm_pop_st_smooth, st_stde),axis=1)	
        with open (str(outfile), "w") as pd:
            np.savetxt(pd, st_pop_stde, fmt='%7.3f', delimiter=' ') 
   
        # matrix with all channels
        if st == 0: 
             time_norm_pop = np.concatenate((timerow, norm_pop_st_smooth),axis=1)
        else:
             time_norm_pop = np.concatenate((time_norm_pop, norm_pop_st_smooth),axis=1)
             
    with open (results_file, "w") as ac:
        np.savetxt(ac, time_norm_pop, fmt='%7.3f',delimiter=' ')   #matrix time pop_channel1 pop_channel_2 ..... pop_channel_n 
    
    return()          
##############################################
     ##########  MAIN   ##########
##############################################
if __name__ == "__main__":
      np.set_printoptions(linewidth  = 150)  # avoid text wrapping in console when printing np.array for checks
      natoms = 12
      lines_per_mol = 14
      AU_TO_FS   = 0.024189
      n_channels = 7
      n_steps    = 2100 # number of simulation steps, +1 since upper limit index is exluded
      timestep   = 10        # 
      movies,geoms=input_check()
      print("Geoms: ",geoms)
      print("#######################\n")
      
      analyze_geoms  = process_movies(movies,geoms)   # np.array returning time, channel over all geoms
      #print(analyze_geoms)
      statistic      = channel_statistics(analyze_geoms)
      
      
      sys.exit(0)
      #distance_matrix(movies)
