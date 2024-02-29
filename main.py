import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('./src/')
from simulation import Simulation
from trajectory import Trajectory

N_traj = 10
N_proc = 7
E_0 = 6.5372045046E-02 #1.5e14 W/cm^2
#E_0 = 0.076
#E_0 = 0.01
omega = 0.057
period = 2*np.pi/omega
u = 0
tau = 4*period
pulse_type = 'cos_2'
t_0_limits = [-0.25*period,0.25*period]
t_0_limits = [-30,30]
p_t_limits = [0,0.3]



sim = Simulation(N_traj,E_0,omega,u,tau,pulse_type,t_0_limits,p_t_limits,N_proc = N_proc)
sim.energy_map()

sim = Simulation(N_traj,E_0,omega,u,tau,pulse_type,t_0_limits,p_t_limits,coulomb = False,N_proc = N_proc)
sim.energy_map()

time = np.linspace(-tau/2,tau/2,200)
sim.plot_E(time)

traj = Trajectory(-2.3,tau/2,0.05,sim.Laser,save_traj = True)
traj.solve()
traj.plot()
print(traj.energy(traj.y[:,0]))
print(traj.final_energy)

traj = Trajectory(-2.3,tau/2,0.05,sim.Laser,save_traj = True, coulomb = False)
traj.solve()
traj.plot()
print(traj.energy(traj.y[:,0]))
print(traj.final_energy)

traj = Trajectory(2.3,tau/2,0.05,sim.Laser,save_traj = True)
traj.solve()
traj.plot()
print(traj.energy(traj.y[:,0]))
print(traj.final_energy)

traj = Trajectory(2.3,tau/2,0.05,sim.Laser,save_traj = True,coulomb = False)
traj.solve()
traj.plot()
print(traj.energy(traj.y[:,0]))
print(traj.final_energy)

plt.show()