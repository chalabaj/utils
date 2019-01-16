# usage python analyze_pop.py 8 100 2100
import math
import sys
import numpy as np
import random
import time
import os
import string
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
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
    populations = np.zeros(shape=(nsteps,nstates))
    timerow = np.arange(1,nsteps+1,1).reshape(nsteps,1) * 0.024 * dt
    stde = np.zeros(shape=(nsteps,nstates))
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
    
    for step in range(0,nsteps):    
        row_norms = populations.sum(axis=1) # row sum is axis=1 opposite to numpe row priority
        row_norms_mat = np.resize(row_norms,(nstates,nsteps)).transpose() # to broadcast the following division, need the smae dimension
        norm_pop = np.divide(populations,row_norms_mat)
        #print(row_norms.reshape(nsteps,1)[step])
        if (not step % 200):
            for st in range(nstates):
                stde[step][st]=math.sqrt(norm_pop[step][st]*(1-norm_pop[step][st])/row_norms.reshape(nsteps,1)[step])*1.96  # binobidal normal dist    
    print(norm_pop.shape,  timerow.shape, row_norms.reshape(nsteps,1).shape)
    alive_mat=np.concatenate((timerow,row_norms.reshape(nsteps,1)),axis=1)	
    final_mat=np.concatenate((timerow,norm_pop),axis=1)	   # add norm
    
    #line = ( str('%.4f ' %time) + ("  ".join("%.3f" %n for n in channel_pop[step])) + "\n")
    for st in range(1,nstates+1):  # 0 is time axis
        outfile="pop_st{}.dat".format(st)
        popst = final_mat[:,[0,st]]
        stdest = stde[:,1].reshape(nsteps,1) # take only individual state due to xmgrace load
        pop_stde = np.concatenate((popst,stdest),axis=1)	
        with open (str(outfile), "w") as pd:
            np.savetxt(pd, pop_stde, fmt='%7.3f', delimiter=' ')   # X is an array
    with open ("alive.dat", "w") as ad:
        np.savetxt(ad, alive_mat, fmt='%7.3f',delimiter=' ')   # X is an array
        
    #WRITE EACH LINE
    
    #sqrt(pop[j,i]*(1-pop[j,i])/anorm[i])*1.96       
    #plt.errorbar(timerow, norm_pop, yerr=stde, fmt='.k');


