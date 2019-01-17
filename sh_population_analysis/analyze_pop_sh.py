# usage python analyze_pop.py 8 100 2100
import math
import sys
import numpy as np
import random
import time
import os
import string
import scipy
from scipy import signal
import numpy as np
##############################################
##############################################
##############################################

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Error: not enought parameters.\nUsage: python ",sys.argv[0]," nstate ntrajs nsteps ")
        sys.exit(1)
    nstates = int(sys.argv[1])  
    ntrajs = int(sys.argv[2])
    nsteps = int(sys.argv[3])
    dt = 10 #timestep 10 au
    populations = np.zeros(shape=(nsteps,nstates),dtype=np.float64)
    timerow = np.arange(1,nsteps+1,1,dtype=np.float64).reshape(nsteps,1) * 0.024 * dt
    stde = np.zeros(shape=(nsteps,nstates),dtype=np.float64)
    print(timerow)
    pop_files = []
    for ntr in range(1,ntrajs+1):
        pop_file = "TRAJ.{}/pop.dat".format(ntr)
        cwd = os.getcwd()
        pop_path = os.path.join(cwd,pop_file)
        if os.path.isfile(pop_path):
            pop_files.append(pop_path)
            print(pop_path)
        else:
            print(pop_path ,' File NOT EXISTS.')    
    
    for popfile in pop_files:                                              # iterate over movies
        step = 0
        print("Processing {}".format(popfile))
        with open(popfile,'r') as pf:
            for line in pf:
                if "CurrentState" in line or "#" in line:
                    print("Found {}. in {}. Skipping".format(line, popfile))
                else:
                    state = int(line.split()[1])-1
                    #print(state)
                    populations[step][state] += 1
                    #alive[step] += 1
                    step += 1
    print("##################\nReading files done \n ####################################")
    
       
    row_norms = populations.sum(axis=1) # row sum is axis=1 opposite to numpe row priority
    row_norms_mat = np.resize(row_norms,(nstates,nsteps)).transpose() # to broadcast the following division, need the smae dimension
    norm_pop = np.divide(populations,row_norms_mat)  # Normalize population
        
    for step in range(0,nsteps):     
        for st in range(nstates):
            if (step == 100) or (not step % (200)):
                stde[step][st]=math.sqrt(norm_pop[step][st]*(1-norm_pop[step][st])/row_norms.reshape(nsteps,1)[step])*1.96  # binobidal normal dist    
   
    alive_mat=np.concatenate((timerow,row_norms.reshape(nsteps,1)),axis=1)	
    with open ("alive.dat", "w") as ad:
        np.savetxt(ad, alive_mat, fmt='%7.3f',delimiter=' ') 
           
    # take only individual state due to xmgrace import XYDY: time pop_st std_st
    for st in range(0,nstates):  # 0 is time axis
        outfile="pop{}.dat".format(st+1)
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
             
    with open ("all_pops.dat", "w") as ac:
        np.savetxt(ac, time_norm_pop, fmt='%7.3f',delimiter=' ')   #matrix time pop_channel1 pop_channel_2 ..... pop_channel_n 
