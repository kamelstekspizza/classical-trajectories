import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from trajectory import Trajectory,Laser
import concurrent.futures #Used for multiprocessing
from multiprocessing import Pool
import os
import time

class Simulation:
    def __init__(self,N_traj,E_0,omega,u,tau,pulse_type,t_0_limits,p_t_limits,**kwargs):
        self.Laser = Laser(E_0,omega,u,tau,pulse_type)
        self.t_0 = np.linspace(t_0_limits[0],t_0_limits[1],N_traj)
        self.p_t = np.linspace(p_t_limits[0],p_t_limits[1],N_traj)
        
        #Keyword arguments
        self.coulomb = kwargs.get('coulomb',True)
        self.N_proc = kwargs.get('N_proc',1)

        self.input = [N_traj,E_0,omega,u,tau,pulse_type,t_0_limits,p_t_limits,self.coulomb,self.N_proc]

    def energy_map(self):
        #fig,ax = plt.subplots()
        #colormap = plt.cm.jet #or any other colormap
        #normalize = matplotlib.colors.Normalize(vmin=-0.5, vmax=0)

        T_0,P_t = np.meshgrid(self.t_0,self.p_t)
        parameters = [(t_0,p_t) for t_0,p_t in np.nditer((T_0,P_t))]
        self.t_0_list = []
        self.p_t_list = []
        self.energy_list = []

        print()
        print('Simulating trajectories... ')

        time_start = time.time()
        with Pool(processes = self.N_proc) as pool:
            results = pool.starmap(self.simulate,parameters)
        
        for result in results:
            if result[2] < 0:
                self.t_0_list.append(result[0])
                self.p_t_list.append(result[1])
                self.energy_list.append(result[2])

        time_end = time.time()
        print()
        print(f'Wall time: {time_end-time_start}')

        #time_start = time.time()
        #with concurrent.futures.ProcessPoolExecutor(max_workers = self.N_proc) as executor:
            #for t_0 in self.t_0:
            #parameters = [[t_0,p_t] for t_0,p_t in np.nditer((T_0,P_t))]
            #for p in parameters:
            #    print(p)
            #results = [executor.submit(self.simulate,*parameter) for parameter in parameters]
            #for f in concurrent.futures.as_completed(results):
                #print('A calculation has finished')

             #   if f.result()[2] < 0:
             #       self.t_0_list.append(f.result()[0]*self.Laser.omega/np.pi*180)
             #       self.p_t_list.append(f.result()[1])
             #       self.energy_list.append(f.result()[2])
                
        """             for p_t in self.p_t:
                trajectory = Trajectory(t_0,self.Laser.tau/2,p_t,self.Laser,coulomb = self.coulomb)
                trajectory.solve()
                if trajectory.final_energy < 0 and trajectory.final_energy > -0.5:
                    self.t_0_list.append(t_0)
                    self.p_t_list.append(p_t)
                    self.energy_list.append(trajectory.final_energy)
                    print(f'(t_0,p_t,E) = ({t_0},{p_t},{trajectory.final_energy})')
                    #ax.scatter(t_0_list,p_t_list,c = energy_list,cmap = colormap,norm= normalize) """

        #time_end = time.time()
        #print(f'Wall time: {time_end-time_start}')

        #sc = ax.scatter(self.t_0_list,self.p_t_list,c = self.energy_list,cmap = colormap,norm= normalize,s = 1) #Make marker color depend on energy!
        #cb = fig.colorbar(sc)

        self.save_output()

        print('Done!')
        print('')

        #trajectory = Trajectory(t_0_list[0],self.Laser.tau/2,p_t_list[0],self.Laser,coulomb = self.coulomb,save_traj = True)
        #trajectory.solve()
        #trajectory.plot()

        #trajectory = Trajectory(t_0_list[-1],self.Laser.tau/2,p_t_list[-1],self.Laser,coulomb = self.coulomb,save_traj = True)
        #trajectory.solve()
        #trajectory.plot()

        return

    def simulate(self,t_0,p_t):
        if self.Laser.pulse_type == 'gaus':
            trajectory = Trajectory(t_0,self.Laser.tau*6,p_t,self.Laser,coulomb = self.coulomb)
        else:
            trajectory = Trajectory(t_0,self.Laser.tau/2,p_t,self.Laser,coulomb = self.coulomb)
        trajectory.solve()
        if trajectory.final_energy < 0:# and trajectory.final_energy > -0.5:
            #print(f'(t_0,p_t,E) = ({t_0},{p_t},{trajectory.final_energy})')
            return (t_0,p_t,trajectory.final_energy)
            
        else:
            return (t_0,p_t,10)
        

    def plot_E(self,time):
        fig_E,ax_E = plt.subplots()
        ax_E.plot(time,self.Laser.E(time))

        return
    
    def make_outputfolder(self):
        #Make a numbered directory XXXX in output directory to store calculation results
        for i in range(1,10000):
            this_run = '{:d}'.format(i).zfill(4)
            output_folder = f'./OUTPUT/{this_run}'
            if not(os.path.isdir((output_folder))):
                os.makedirs(output_folder)
                self.output_folder = output_folder
                break
        return

    def save_output(self):
        #Save the output of calculations, and copy input file to directory where output is stored
        self.make_outputfolder()
        print('')
        print(f'Saving output in {self.output_folder}')
        print('')

        self.DATA = np.transpose(np.array([self.t_0_list,self.p_t_list,self.energy_list]))
        np.savetxt(f'{self.output_folder}/output.dat',self.DATA,header = 't_0 [a.u.], p_t [a.u.], Energy [a.u.]')

        with open(f'{self.output_folder}/input.dat','w') as file:
            for line in self.input:
                file.write(str(line)+'\n')

        return
