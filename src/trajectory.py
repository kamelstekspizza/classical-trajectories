import numpy as np
import scipy.integrate as si
import matplotlib.pyplot as plt


class Trajectory:
    def __init__(self,t_0,t_end,p_t,Laser_in,**kwargs):
        #Required arguments  
        self.t_0 = t_0
        self.t_end = t_end
        self.p_t = p_t
        self.Laser = Laser_in

        #keyword arguments
        self.coulomb = kwargs.get('coulomb',True)
        self.save_traj = kwargs.get('save_traj',False)


        if self.coulomb:
            self.RHS = self.RHS_coulomb
        else:
            self.RHS = self.RHS_no_coulomb

        return
        

    def RHS_coulomb(self,t,y):
        r_3 = np.sqrt(y[0]**2+y[1]**2)**3
        RHS = np.zeros(4)
        RHS[0] = y[2]
        RHS[1] = y[3]
        RHS[2] = -y[0]/r_3 -self.Laser.E(t)
        RHS[3] = -y[1]/r_3

        return RHS
    
    def RHS_no_coulomb(self,t,y):
        RHS = np.zeros(4)
        RHS[0] = y[2]
        RHS[1] = y[3]
        RHS[2] = -self.Laser.E(t)
        return RHS
    
    def energy(self,y):
        r = np.sqrt(y[0]**2+y[1]**2)
        potential = -1.0/r
        kinetic = 0.5*(y[2]**2+y[3]**2)

        return potential + kinetic

    def solve(self):
        if np.abs(self.Laser.E(self.t_0)) < 1e-9:
            #Electric field is too small for tunneling!
            self.final_energy  = 10
            return
        
        limits = (self.t_0,self.t_end)
        y_0 = np.zeros(4)
        y_0[0] = -(0.5 + 0.5*self.p_t**2)/self.Laser.E(self.t_0) #Effective tunneling barrier, 0.5 is the I_p of hydrogen
        #y_0[0] = -(0.5 + np.sqrt(0.5**2-4*self.Laser.E(self.t_0)))/(2*self.Laser.E(self.t_0)) #Effective tunneling barrier, 0.5 is the I_p of hydrogen, Eichmann papers
        #y_0[0] = 0.5/self.Laser.E(self.t_0) #Effective tunneling barrier, 0.5 is the I_p of hydrogen
        #y_0[2] = -self.Laser.A(self.t_0)
        y_0[3] = self.p_t
        solution = si.solve_ivp(self.RHS,limits,y_0,method = 'RK45',max_step = 0.25*self.Laser.period)

        

        self.final_energy = self.energy(solution.y[:,-1])
        
        if self.save_traj:
            self.time = solution.t
            self.y = solution.y

    def coulomb_orbit(self):
        E = self.energy(self.y)
        ang_mom = self.y[0]*self.y[3]
        e = np.sqrt(1+ 2*E*ang_mom**2)

    def plot(self):
        fig_pos,ax_pos = plt.subplots()
        fig_mom,ax_mom = plt.subplots()
        fig_energy,ax_energy = plt.subplots()

        ax_pos.plot(self.y[0],self.y[1])
        ax_mom.plot(self.time,self.y[2])
        ax_energy.plot(self.time,self.energy(self.y))

        return
    
    
    
    


class Laser:
    def __init__(self,E_0,omega,u,tau,pulse_type):
        self.E_0 = E_0
        self.omega = omega
        self.u = u
        self.tau = tau
        self.pulse_type = pulse_type

        self.A_0 = self.E_0/self.omega
        self.period = 2*np.pi/self.omega

        if self.pulse_type == 'mono':
            self.A = self.A_mono
            self.E = self.E_mono
        elif self.pulse_type == 'cos_2':
            self.A = self.A_cos_2
            self.E = self.E_cos_2        
        elif self.pulse_type == 'cos_6':
            self.A = self.A_cos_6
            self.E = self.E_cos_6

        elif self.pulse_type == 'gaus':
            self.A = self.A_gaus
            self.E = self.E_gaus

    def A_mono(self,t):
       return -self.E_0/self.omega*np.sin(self.omega*t)
         
    def A_cos_2(self,t):
        return self.E_0/(self.omega)*np.sin(self.omega*t+self.u)*(np.cos(np.pi*t/self.tau))**2

    def A_cos_6(self,t):
        return self.E_0/(self.omega)*np.sin(self.omega*t+self.u)*(np.cos(np.pi*t/self.tau))**6

    def A_gaus(self,t):
        return self.A_0*np.sin(self.omega*t+self.u)*np.exp(-2.0*np.log(2)*(t/self.tau)**2)

    def E_mono(self,t):
        return self.E_0*np.cos(self.omega*t)
    
    def E_cos_2(self,t):
        return self.E_0/(self.omega)*(-self.omega*np.cos(self.omega*t+self.u)*(np.cos(np.pi*t/self.tau))**2+2*np.pi/self.tau*np.sin(self.omega*t +self.u)*np.sin(np.pi*t/self.tau)*np.cos(np.pi*t/self.tau))

    def E_cos_6(self,t):
        return self.E_0/(self.omega)*(-self.omega*np.cos(self.omega*t+self.u)*(np.cos(np.pi*t/self.tau))**2+6*np.pi/self.tau*np.sin(self.omega*t +self.u)*np.sin(np.pi*t/self.tau)*np.cos(np.pi*t/self.tau)**5)

    def E_gaus(self,t):
        return self.A_0*np.exp(-2.0*np.log(2)*(t/self.tau)**2)*(-self.omega*np.cos(self.omega*t+self.u) + 4.0*np.log(2)*t/self.tau**2*np.sin(self.omega*t+self.u))

