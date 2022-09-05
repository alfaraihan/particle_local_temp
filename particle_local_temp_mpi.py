from ovito import data
from ovito.io import import_file
from ovito.data import CutoffNeighborFinder
import numpy as np
from numpy import (save,load)
from mpi4py import MPI
import time
import datetime

def particle_mass(ptype, mass):
    for i in range(len(mass)):
        if ptype==i+1: m=mass[i]
        else: continue
    return m

def units_conv(units):
    if units=="metal":
        v_conv=1e02
        m_conv=1.66054e-27
    if units=="real":
        v_conv=1e05
        m_conv=1.66054e-27
    return v_conv, m_conv
    

def local_par_temp(data, saveto, cutoff=5.1, units="metal", mass=[28.0855,15.999,1.00794,15.999], npy=True):
    #var
    cutoff=cutoff
    v_conv, m_conv=units_conv(units)

    # MPI
    comm=MPI.COMM_WORLD
    rank=comm.rank
    size=comm.size

    tstart=time.perf_counter()

    # file
    openfile=data
    savefile=saveto
    
    if rank==0:print("Local Partilce Temperature Finder \nUsing ", size, "processors")

    pipeline=import_file(openfile)

    for frame_index in range(pipeline.source.num_frames):
        if rank==0:print("\nFrame "+str(frame_index+1)+"/"+str(pipeline.source.num_frames))
        tstart2=time.perf_counter()

        data=pipeline.compute(frame_index)

        npart=data.particles.count
        pos=data.particles.position
        ptypes=data.particles.particle_type
        v=data.particles.velocities

        ntask=int(np.ceil(npart/size))

        if rank==0:print("read from : ", openfile,", npart = ", npart)

        finder = CutoffNeighborFinder(cutoff, data)
        vpmfinder = CutoffNeighborFinder(10, data)

        Ti=[]

        for rank_i in range(size):
                if rank==rank_i:
                    for index in range(rank_i*ntask,rank_i*ntask+ntask):
                        if index < npart :
                            if rank==0:
                                #progress
                                progress=index/(ntask-1)*100
                                if round(progress%2,4)==0:
                                    bar="⊂(◉‿◉)つ  Progress :" + str(int(progress)) +" %"
                                    print(bar)


                            vmx, vmy, vmz = 0, 0, 0
                            vpm_x, vpm_y, vpm_z= 0, 0, 0
                            mtot=0

                            ## Center of mass ## 
                            # center of mass velocity of a particle is based on average velocity of surrounding particles
                            # divided by total mass of surrounding particles

                            # Iterate over the neighbors of the current particle:
                            for neigh in finder.find(index):
                                m = particle_mass(ptypes[neigh.index], mass)

                                vi= v[neigh.index]*v_conv
                                vmx+=vi[0]*m 
                                vmy+=vi[1]*m
                                vmz+=vi[2]*m
                                mtot+=m

                            vpm_x=vmx/mtot if mtot!=0 else 0
                            vpm_y=vmy/mtot if mtot!=0 else 0
                            vpm_z=vmz/mtot if mtot!=0 else 0

                            Nn=0 
                            mv=0
                            kb=1.38064852e-23 # boltzman's constant

                            ## Temperature ##

                            for neigh in finder.find(index):

                                m = particle_mass(ptypes[neigh.index], mass)*m_conv

                                vi= v[neigh.index]*v_conv

                                mv+=m*np.power(vi[0]-vpm_x,2) #1/2 m v^2 
                                mv+=m*np.power(vi[1]-vpm_y,2)
                                mv+=m*np.power(vi[2]-vpm_z,2)

                                Nn+=1

                            Ti.append(mv/(3*Nn*kb) if Nn!=0 else 0)

        gTi = comm.gather(Ti, root=0)

        if rank==0:
            if npy:
                tTi=np.array([])
                for index_write in range(len(gTi)):
                    tTi=np.append(tTi,gTi[index_write]) 

                save(savefile+str(frame_index)+".npy", tTi)  
                print("saved to : "+savefile+str(frame_index)+".npy")

            file=open(savefile+str(frame_index)+".xyz", "w")
            file.write("%d \n"%(npart))
            file.write("Temp data xyz: ptype x y z temp\n")
            for file_index in range(npart):
                file.write("%d %f %f %f %f\n"%(ptypes[file_index],pos[file_index][0],pos[file_index][1],pos[file_index][2],tTi[file_index]))
            print("saved to : " + savefile+str(frame_index)+".xyz   ")
            tstop2=time.perf_counter()
            print(str(datetime.timedelta(seconds=int(tstop2-tstart2))))


    if rank==0:
        tstop=time.perf_counter()  
        print("Total time : ", datetime.timedelta(seconds=int(tstop-tstart)), "(ಠ益ಠ)")
    
    return gTi
