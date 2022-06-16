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
        
        self.Time = np.load(dataSetFile + '/Time_arr.npy')
        
        # params = [ L1 , L2 , M1 , M2 , M1+M2 , g ]
        self.Params = np.load(dataSetFile + '/Parameters.npy',allow_pickle=True)
    def Poincare(self,):
        ''' To be Implemented . Determine the Poincare sections of the temporal distribution array'''
        poincare = (self.dPendulum_data[:,1,:-1] < 0) & (self.dPendulum_data[:,1,1:] >0)
        return poincare
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
    '''
    dataObj = DB_Post_Analysis('/home/jabejar/Documents/Datasets/ds_000')
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
    '''
    

    dataObj = DB_Post_Analysis('/home/jabejar/Documents/Datasets/ds_002')
    a = dataObj.dPendulum_data
    b = dataObj.Time
    c = dataObj.Params
    
    poincare_bool = (a[:,1,:-1] < 0) & (a[:,1,1:] >0)
    
    poincare_index = np.where(poincare_bool == True)
    
    poincare_index = np.array([poincare_index[0],poincare_index[1]])
    poincare_sum = np.sum(poincare_bool,axis = -1)
    pcMin = min(poincare_sum)
    pcDiff = (poincare_sum - pcMin).tolist()
    pcDiff.insert(0,0)
    poincare_index_adj = []
    for i in range(max(poincare_index[0])):
        print(sum(pcDiff[:i+1]))
        poincare_index_adj.append(poincare_index[:,i*pcMin + sum(pcDiff[:i+1]):pcMin*(i+1)+ sum(pcDiff[:i+1])])
        
    poincare_index_adj = np.concatenate(poincare_index_adj,axis = -1).reshape(2,max(poincare_index[0])*pcMin)
    
    
    poincare = a[poincare_index_adj[0],0,poincare_index_adj[1]].reshape(max(poincare_index[0]),pcMin)
    poincare2 = a[poincare_index_adj[0],2,poincare_index_adj[1]].reshape(max(poincare_index[0]),pcMin)
    
    plt.figure()
    plt.plot(poincare,poincare2,'.',markersize = 1)
    
    
    
    
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    