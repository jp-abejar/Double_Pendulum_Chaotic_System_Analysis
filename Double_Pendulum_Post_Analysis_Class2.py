    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 21:46:47 2022

@author: jabejar
"""
import numpy as np
import matplotlib.pyplot as plt


class DB_Post_Analysis():
        def __init__ (self,dataSetFile):
            self.dPendulum_data = np.load(dataSetFile + '/db_dataset.npy').swapaxes(1,2)
            self.DatShape = self.dPendulum_data.shape
            self.Time = np.load(dataSetFile + '/Time_arr.npy')
            
        
            self.Params = np.load(dataSetFile + '/Parameters.npy',allow_pickle=True)
            
            
        def Poincare(self, dreplace = None,MassId = 1):
            '''Determine the Poincare sections of the temporal distribution array'''
            if dreplace is None:
                poincare_bool = (self.dPendulum_data[:,MassId,:-1] < 0) & (self.dPendulum_data[:,MassId,1:] > 0) & (self.dPendulum_data[:,MassId,:-1] < 2)
            else:
                poincare_bool = (dreplace[:,MassId,:-1] < 0) & (dreplace[:,MassId,1:] > 0) & (dreplace[:,MassId,:-1] < 2)
        
            poincare_index = np.where(poincare_bool == True) # returns the index where the poincare points occur
            
            poincare_index = np.array([poincare_index[0],poincare_index[1]]) # convert to array for easier analysis
            print(poincare_index)
            poincare_sum = np.sum(poincare_bool,axis = -1) # get the sum of poncare sections per row
            pcMin = min(poincare_sum) # find the minimum crosses for all rows
            pcDiff = (poincare_sum - pcMin).tolist() # get the excess amount of values per row
            pcDiff.insert(0,0) # insert 0 to beginning
            poincare_index_adj = [] 
            mpci = max(poincare_index[0])+1
            
            
            #Removing excess poincare values to make the array symmetric
            for i in range(mpci):
                poincare_index_adj.append(poincare_index[:,i*pcMin + sum(pcDiff[:i+1]):pcMin*(i+1)+ sum(pcDiff[:i+1])])
                
                
            #Concatenate the data and reshape to index the main data array
            poincare_index_adj = np.concatenate(poincare_index_adj,axis = -1).reshape(2,mpci*pcMin)
            print(poincare_index_adj)
            
            
            # Slice the main data array to grab poincare sections
            poincare = self.dPendulum_data[poincare_index_adj[0],0,poincare_index_adj[1]+1].reshape(mpci,pcMin)
            poincare2 = self.dPendulum_data[poincare_index_adj[0],1,poincare_index_adj[1]+1].reshape(mpci,pcMin)
            poincare3 = self.dPendulum_data[poincare_index_adj[0],2,poincare_index_adj[1]+1].reshape(mpci,pcMin)
            poincare4 = self.dPendulum_data[poincare_index_adj[0],3,poincare_index_adj[1]+1].reshape(mpci,pcMin)
            
            return np.array([poincare,poincare2,poincare3,poincare4])
        
        def RPCA(self,):
            ''' For future use'''
            from RPCA3 import rpca
            poincareData = dataObj.Poincare()
            org_shape = poincareData.shape
            rpcaInitData = poincareData.reshape(1,org_shape[0],org_shape[1]*org_shape[2])[0]
            
            rpca_dat = rpca(M = rpcaInitData.T,max_iter=100,lam = 0.01)
            
            return rpca_dat,org_shape
        
        
        def Fourier_analysis(self,):
            fftv = []
            
            for i in range(4):
                if len(self.dPendulum_data.shape) ==3:
                    fftv.append(np.fft.fftshift(np.fft.fft(self.dPendulum_data[:,i,:],axis = -1)))
                else:
                    fftv.append(np.fft.fftshift(np.fft.fft(self.dPendulum_data[i,:],axis = -1)))
            
            
            freq = np.fft.fftshift(np.fft.fftfreq(self.Time.shape[0],d = self.Time[1]-self.Time[0]))
            
            return fftv,freq
    

if __name__ == "__main__":
        
        dataObj = DB_Post_Analysis('/home/jabejar/Documents/GitHub/Double_Pendulum_Simulation_Analysis/Datasets/ds_000')
        a = dataObj.dPendulum_data
        b = dataObj.Time
        c = dataObj.Params
        
        fft,freqs = dataObj.Fourier_analysis()
        
        for i in range(a.shape[0]):
            plt.figure(5)
            plt.plot(b,a[i,0,:],'.',markersize = 0.3)
        
            
            
            for j in range(4):
                plt.figure(j)
                plt.plot(freqs,abs(fft[j][i]))
                plt.xlim(0,5)
        
        
        
        
        dataObj = DB_Post_Analysis('/home/jabejar/Documents/GitHub/Double_Pendulum_Simulation_Analysis/Datasets/ds_005')
        a = dataObj.dPendulum_data
        b = dataObj.Time
        c = dataObj.Params
        
        poincareData = dataObj.Poincare()
        
        plt.figure()
        
        for i in range(len(c[2])):
            plt.subplot(12,12,(i%144)+1)
            plt.plot(poincareData[0][i].T,poincareData[2][i].T,'.',markersize = 1)
        
        plt.plot(poincareData[0][::1].T,poincareData[2][::1].T,'.',markersize = 1)
        
        plt.figure()
        plt.plot(c[2],poincareData[2],'.b',markersize = 0.4)
        
        plt.figure()
        plt.plot(c[2],np.std(poincareData[2],axis = -1),'.b',markersize = 0.4)
        
       
        plt.figure()
        plt.plot(b,a[0,0,:])
        
   
      
        (L,S,(U,s,V),iter_,lam,err),org_shape = dataObj.RPCA()
        
        
        
        S = L.T.reshape(org_shape)
        
        
        plt.figure()
        plt.plot(S[0][::200],S[2][::200],'.',markersize = 0.4)
        
        
            
            
            
            
            
 