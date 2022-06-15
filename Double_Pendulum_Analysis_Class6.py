#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 21:23:38 2022

@author: jabejar
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from threading import Thread
import time
import os
class Double_Pendulum():
    
    def __init__(self,d = 'dMass',dt =0.1, maxTime = 5,maxMass = 10, minMass = 0.1,massFraction = 1 \
                 , maxDegree =90, minDeg = 0.1,degFraction = 1, maxCrossings = 10 \
                     ,M1 =1,M2 =1,L1 =1,L2 = 1,g = 9.8,Tht1 = 1,Tht2 = 1,Anim = False):
        self.Anim = Anim # bool
        self.ani = None # intended to hold the animation function constant 
        self.Data = None # Data array if data storage is desired
        
        self.Lock = True
        self.d = d
        self.dt = dt
        self.tmax = maxTime
        self.numtimes=int(self.tmax/self.dt)
        self.times = np.linspace(0, self.tmax, self.numtimes)


        self.dof=4 #Number of variables needed
        self.maxCrossings = maxCrossings




        self.Mass1 = np.linspace(minMass,maxMass,int((maxMass- minMass)/massFraction)) if d == 'dMass' else M1 if d == 'Single' else 1
        
        if d == 'dTheta1' or d == 'dTheta2':
            self.ThetaVals = np.linspace(minDeg, maxDegree,int((maxDegree - minDeg)/degFraction)) *(np.pi/180)


        try:
            self.Instances = self.Mass1.shape[0]
        except:
            try:
                self.Instances = self.ThetaVals.shape[0]
            except:
                self.Instances = 1
       
        self.newVariables = np.zeros([self.Instances,self.dof], dtype = float)
        self.oldVariables = np.zeros([self.Instances,self.dof], dtype = float)
        self.energies = np.zeros([self.numtimes,self.Instances],dtype=float)  
        self.AngularVel= np.zeros([self.Instances,self.maxCrossings],dtype = float )
        self.poincare = np.zeros((self.Instances,self.maxCrossings,2),dtype = float)



        self.M1 = self.Mass1

    
        if d == 'dTheta1' or d == 'dTheta2':
            self.theta1 = 1 if d == 'dTheta1'  else 0
            self.theta2 = 1 if d == 'dTheta2'  else 0
            self.omega1 = 0.0
            self.omega2 = 0.0

        else:
            self.theta1 = Tht1*(np.pi/180)
            self.theta2 = Tht2*(np.pi/180)
            self.omega1 = 0.0
            self.omega2 = 0.0



        #Parameters
        self.M2 = M2 
        self.L1 = L1
        self.L2 = L2
        self.M = self.M1 + self.M2
        self.g = g

       



        self.params = [self.L1,self.L2,self.M1,self.M2,self.M,self.g]
      
        self.oldVariables=np.array([self.theta1,self.theta2,self.omega1,self.omega2]) * np.ones((self.Instances,self.dof))#The first values of theta and its derivative will be the initial conditions

        if not d == 'dMass' :
            if d == 'dTheta1' or d == 'dTheta2':
                self.oldVariables[:,0 if d == 'dTheta1'  else 1] = self.ThetaVals  
                self.energies[0,:] = (self.M)*self.g*self.L1*(1-np.cos(self.theta1 if d == 'dTheta2'  else self.ThetaVals)) + self.M2*self.g*self.L2*(1-np.cos(self.theta2 if d == 'dTheta1'  else self.ThetaVals))
            else:
                self.energies[0,:] = (self.M)*self.g*self.L1*(1-np.cos(self.theta1)) + self.M2*self.g*self.L2*(1-np.cos(self.theta2))
        else:
            self.energies[0,:] = (self.M)*self.g*self.L1*(1-np.cos(self.theta1)) + self.M2*self.g*self.L2*(1-np.cos(self.theta2))
            



        
        self.step = 0
        self.E =1




        self.counts = np.zeros(self.Instances,dtype=np.int16)

        self.prevCounts = np.zeros(self.Instances,dtype=np.int16)

        self.ValidCounts = np.ones(self.Instances,dtype=np.int16)
        
        
        
    def RK4(self,numPoints = None , saveDat = False):
        if saveDat and numPoints is not None:
            self.Data = np.zeros([self.Instances,int(numPoints),self.dof],dtype=float) 
        
        while sum(self.ValidCounts) != 0 or (numPoints is not None and self.step < numPoints ) :
            
            self.Lock = True
            self.ValidCounts = np.less(self.counts,self.maxCrossings)
            #The next four lines calculate derivatives at different times, and multiply the derivatives by the time step, to get the changes in the variables
            k1 = self.dt*self.derivatives(self.oldVariables)
            k2 = self.dt*self.derivatives(self.oldVariables+k1/2)
            k3 = self.dt*self.derivatives(self.oldVariables+k2/2)
            k4 = self.dt*self.derivatives(self.oldVariables+k3)
            #Now that we've done the RK method and solved the ODE with the old time and data at the previous step, we up the time and step and store for the new time
            
            self.newVariables = self.oldVariables+k1/6+k2/3+k3/3+k4/6 #This is the Runge-Kutta estimate of the new values of theta and d(theta)/dt
            self.newVariables[:,0:2] = self.newVariables[:,0:2] - 2*np.pi*np.greater(self.newVariables[:,0:2], np.pi)
            self.newVariables[:,0:2] = self.newVariables[:,0:2] + 2*np.pi* np.less(self.newVariables[:,0:2], -np.pi)
            
            if not self.Anim:
                self.PoincareFunc(MassVar=1)
                if self.E < self.numtimes-1:
                    self.energies[self.E,:] = self.energy(self.newVariables) #Another place where science shows up
                    self.E += 1
                
                
                
                
            else:
                
                self.oldVariables = self.newVariables
                self.Lock = False
                
                time.sleep(0.1)
                
                
            if numPoints is not None and self.step < numPoints:
                self.Data[:,self.step,:] = self.newVariables
                self.step +=1
                if self.step%1000 == 0:
                    print(self.step)
                
               
                
               
        if saveDat:
            
            #This block genreates a new directory for the dataset
            i = 0
            DoneMD = False
            while not DoneMD:
                try:
                    os.mkdir('/home/jabejar/Documents/Datasets/ds_%003i'%(i))
                    DoneMD = True
                except:
                    i += 1
                    
                    
            
            np.save('/home/jabejar/Documents/Datasets/ds_%003i'%(i) + '/Parameters',np.array(self.params,dtype = object))
            
            if self.d == 'dTheta1' or self.d == 'dTheta2':
                np.save('/home/jabejar/Documents/Datasets/ds_%003i'%(i) + '/Theta1_arr' if self.d == 'dTheta1' else 'Theta2_arr',self.ThetaVals)
                
            np.save('/home/jabejar/Documents/Datasets/ds_%003i'%(i) + '/db_dataset',self.Data)
            np.save('/home/jabejar/Documents/Datasets/ds_%003i'%(i) + '/Time_arr',np.linspace(0,numPoints*self.dt,int(numPoints)))
       
        
       
    def PoincareFunc(self,MassVar = 1):
        
        # Keeping the angles between pi and -pi
        # self.newVariables[:,0:2] = self.newVariables[:,0:2] - 2*np.pi*np.greater(self.newVariables[:,0:2], np.pi)
        # self.newVariables[:,0:2] = self.newVariables[:,0:2] + 2*np.pi* np.less(self.newVariables[:,0:2], -np.pi)
        
        
        self.counts += np.less(self.oldVariables[:,MassVar],0) * np.greater(self.newVariables[:,MassVar],0) * np.less(abs(self.oldVariables[:,MassVar]),2.0) * self.ValidCounts
        self.cross = np.where( self.counts -self.prevCounts == 1)[0]
        if self.cross.size > 0:
            print(sum(self.ValidCounts))
            
            self.poincare[self.cross,self.counts[self.cross]-1,:] = self.newVariables[self.cross,0 if MassVar else 1][0],self.newVariables[self.cross,2 if MassVar else 3][0]
                
            self.AngularVel[self.cross,self.counts[self.cross]-1] = self.newVariables[self.cross,3]
        
        
        self.oldVariables = self.newVariables
        self.prevCounts = self.counts*1
        
        
    def derivatives(self,variables,params = None,t = None):
        
        

        
        #Unpack variables
        theta1 = variables[:,0]
        theta2 = variables[:,1]   
        omega_1 = variables[:,2]
        omega_2 = variables[:,3]

        
        
        
        omega1 = omega_1
        omega2 = omega_2    
       
        omega1dot = (-self.M2*self.L1*(omega1**2)*np.sin(theta1-theta2)*np.cos(theta1-theta2) + self.g*self.M2*np.sin(theta2)*np.cos(theta1-theta2) - self.M2*self.L2*(omega2**2)*np.sin(theta1-theta2) - (self.M)*self.g*np.sin(theta1))/(self.L1*(self.M) - self.M2*self.L1*(np.cos(theta1-theta2)**2))

        omega2dot = (self.M2*self.L2*(omega2**2)*np.sin(theta1-theta2)*np.cos(theta1-theta2) + self.g*(self.M)*np.sin(theta1)*np.cos(theta1-theta2) + (self.M)*self.L1*(omega1**2)*np.sin(theta1-theta2) - (self.M)*self.g*np.sin(theta2))/(self.L2*(self.M) - self.M2*self.L2*(np.cos(theta1-theta2)**2))

        return np.array([omega1,omega2,omega1dot,omega2dot]).T


    def energy(self,variables,params = None, t = None):
       
        #Unpack variables   
        theta1 = variables[:,0]
        theta2 = variables[:,1]   
        omega1 = variables[:,2]
        omega2 = variables[:,3]


        
        KE = 0.5*(self.M)*(self.L1**2)*(omega1**2) + self.M2*self.L1*self.L2*omega1*omega2*np.cos(theta1-theta2) + 0.5*self.M2*(self.L2**2)*(omega2**2) 
        U = (self.M)*self.g*self.L1*(1-np.cos(theta1)) + self.M2*self.g*self.L2*(1-np.cos(theta2))
        
        return KE+U 
    
    
    def PlotBifurcation(self,):
        
        if self.d != 'Single':
            plt.figure()
            plt.title('Bifurcation Diagram')
            plt.xlabel('Mass1 [g]' if self.d == 'dMass' else r'$\theta$1 [rad]' if self.d == 'dTheta1' else r'$\theta$2 [rad]')
            plt.ylabel(r'$\omega$2 [rad/s]')
            plt.plot(self.Mass1 if self.d == 'dMass' else self.ThetaVals ,self.AngularVel,'.b',markersize=0.8)
            
            
    def PlotPoincare(self, n = 1):
        plt.figure()
        plt.title('Poincare Section for when Mass 2 goes through '+''r'$\theta$1 = 0')
        plt.xlabel(r'$\theta$1 [rad]')
        plt.ylabel(r'$\omega$1 [rad/s]')
        for i in range(0,self.Instances,n):
            plt.plot(self.poincare[i,:,0],self.poincare[i,:,1],'.',markersize = 1)
            
            
    
    def AnimatePendulum(self,):
        def init():
            self.Xp1.append(0)
            self.Yp1.append(0)
            
            self.Xp2.append(0)
            self.Yp2.append(0)
            
            maxL = self.L1 +self.L2
            self.ax.set_xlim(-maxL-1, maxL +1)
            self.ax.set_ylim(-maxL-1, maxL +1)
            return self.ln,self.ln2,self.ln3,self.ln4,
        
        def update(frame):
            if not self.Lock:
                # print(self.newVariables[0,0])
                self.Xp1[0] = -self.L1*np.sin(self.newVariables[0,0])
                self.Yp1[0] = -self.L1*np.cos(self.newVariables[0,0])
                
                
                self.Xp2[0] = -self.L2*np.sin(self.newVariables[0,1]) +self.Xp1[0] 
                self.Yp2[0] = -self.L2*np.cos(self.newVariables[0,1]) + self.Yp1[0] 
                
            # print(self.Xp1)
            self.ln.set_data(self.Xp1[0], self.Yp1[0])
            self.ln2.set_data(self.Xp2[0], self.Yp2[0])
            
            self.ln3.set_data([0,self.Xp1[0]],[ 0, self.Yp1[0]])
            self.ln4.set_data([self.Xp1[0], self.Xp2[0]],[self.Yp1[0], self.Yp2[0]])
            return self.ln,self.ln2,self.ln3,self.ln4,
        if self.Anim:
            
            
            
            
            self.figAnim, self.ax = plt.subplots()
            self.Xp1, self.Yp1 , self.Xp2, self.Yp2 = [], [] , [], []
            
            self.ln, = plt.plot([], [], 'ro',label = 'Mass1 = {}g'.format(self.M1))
            self.ln2, = plt.plot([], [],'bo',label = 'Mass2 = {}g'.format(self.M2))
            plt.title('Double Pendulum Simulation')
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.legend()          
            self.ln3, = plt.plot([], [],'k-')
            self.ln4, = plt.plot([], [],'k-')
            self.ani = FuncAnimation(self.figAnim, update, init_func=init,fargs =(),frames=None, blit=True,interval = 10)
        else:
            print("Animation Parameter Disabled")
        

        
if __name__ == "__main__":
    print('Select:\nTest1:1\nTest2:2\nTest3:3\nTest4:4\nTest5:5\n\n')
    try:
        input_TestVal = int(input('Your Selection:'))
    except:
        input_TestVal = 10
    if input_TestVal == 1:
        """ Tests a single pendulum case and plots an animated simulation"""
        
        pendulumSet1 = Double_Pendulum(dt = 0.1, d = 'Single',Tht1=180, Tht2 = 168,M1 = 10,M2 = 0.3,L1 =2,L2=5,Anim = True)
        thrd1 = Thread(target = pendulumSet1.RK4, args = (10,True,))
        thrd1.start()
        pendulumSet1.AnimatePendulum()
       
    elif input_TestVal == 2:
        """ Tests a dMass case where we plot the poincare sections for N mass values as well as the resulting Bifurcation diagram  """
        
        pendulumSet1 = Double_Pendulum(dt = 0.01, d = 'dMass',maxMass = 40, minMass = 10,massFraction = 0.01,Tht1=10, Tht2 = 125,maxCrossings=100)
        pendulumSet1.RK4()
        pendulumSet1.PlotBifurcation()
        pendulumSet1.PlotPoincare()
        
    elif input_TestVal == 3:
        """ Tests a dTheta1 case where we plot the poincare sections for N Theta1 values as well as the resulting Bifurcation diagram  """
        
        pendulumSet1 = Double_Pendulum(dt = 0.1, d = 'dTheta1',maxDegree =90, minDeg = 0.1,degFraction = 1, maxCrossings = 100)
        pendulumSet1.RK4()
        pendulumSet1.PlotBifurcation()
        pendulumSet1.PlotPoincare()
        
    elif input_TestVal == 4:
        """ Tests a dTheta2 case where we plot the poincare sections for N Theta2 values as well as the resulting Bifurcation diagram  """
        
        pendulumSet1 = Double_Pendulum(dt = 0.1, d = 'dTheta2',maxDegree =90, minDeg = 0.1,degFraction = 1, maxCrossings = 100,M1 = 1,M2 = 1)
        pendulumSet1.RK4()
        pendulumSet1.PlotBifurcation()
        pendulumSet1.PlotPoincare()
    
    elif input_TestVal == 5:
        """ Tests the save data feature where we input tyhe numPoints we want saved into the RK4 function.
            Note that the data can be accessed through the object after the RK4 algorithm terminates"""
            
        pendulumSet1 = Double_Pendulum(dt = 0.01, d = 'dMass',maxMass = 20, minMass = 0.1,massFraction = 0.01,Tht1=10, Tht2 = 125,maxCrossings=1)
        pendulumSet1.RK4(numPoints=100_000,saveDat=True)
        pendulumSet1.PlotBifurcation()
        pendulumSet1.PlotPoincare()
    
    
    else:
        print('Invalid Selection')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    