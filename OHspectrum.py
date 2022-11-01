# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 22:01:44 2020

@author: Kenneth
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob as glob
from natsort import natsorted
from scipy import interpolate


Jnum=15  #Number of lines derived

# #OH(3,1) vibrational band
# lambda0=15046.6 #nm

# Bnu=np.array([17.807,16.414])
# Dnu=np.array([0.0018,0.0018])
# Ynu=np.array([-7.876,-8.568])
# Anu=np.array([-140.2,-140.6])
# qnu=np.array([0.0399,0.0351])
# Gnu=np.array([5417.41,12061.61])
# DeltaG=np.array([3569.63,10213.83])  #Change in G from nu=0 state


#OH(2,0) Band

Bnu=np.array([18.515,17.108])
Dnu=np.array([0.0018,0.0018])
Ynu=np.array([-7.547,-8.214])
Anu=np.array([-139.7,-140.5])
qnu=np.array([0.0417,0.0377])
Gnu=np.array([1847.78,8821.36])
DeltaG=np.array([0.0,6973.58])  #Change in G from nu=0 state


F1lower=np.zeros((Jnum))    #3/2 states
F2lower=np.zeros((Jnum))    #1/2 states

F1upper=np.zeros((Jnum))    #3/2 states
F2upper=np.zeros((Jnum))    #1/2 states

Erel=np.zeros((2,Jnum))   #Energy difference for each K as functions of ONLY K equations
Intensity=np.zeros((2,Jnum))
Intdouble=np.zeros(Jnum)
Pint=np.zeros((2,2*(Jnum-1)))
Qint=np.zeros((2,2*(Jnum-1)))
Rint=np.zeros((2,2*(Jnum-1)))
WP=np.zeros((2,2*(Jnum-1)))
WQ=np.zeros((2,2*(Jnum-1)))
WR=np.zeros((2,2*(Jnum-1)))

WG=np.zeros((2,2*(Jnum-1)))
Gint=np.zeros((2,2*(Jnum-1)))
Gb=np.zeros((2,2*(Jnum-1)))

Pb=np.zeros((2,2*Jnum-2))
Rb=np.zeros((2,2*Jnum-2))
Qb=np.zeros((2,2*Jnum-2))

Pcross=np.zeros((2,2*Jnum-2))
Wcross=np.zeros((2,2*Jnum-2))
Cint=np.zeros((2,2*Jnum-2))

Q=np.zeros((2,Jnum))
Part=np.zeros(2)


Fdouble=np.zeros((2,Jnum))

Pij=np.array([1.62,2.9,4.04,5.14,6.21,7.26,8.3,9.32,10.34,11.36,12.37,13.38,14.39,15.4,16.41,17.42,18.42,19.43,20.43,21.44])

Temp=190.0 #K

for J in range(0,Jnum):
    J1=(J*2+2)/2.0
#    F1lower[0,J]=Bnu[0]*((J1+0.5)**2-1-0.5*np.sqrt(4*(J1+0.5)**2+Ynu[0]*(Ynu[0]-4.0)))-Dnu[0]*J1**4
    J2=(J*2+1)/2.0
#    F2lower[0,J]=Bnu[0]*((J2+0.5)**2-1+0.5*np.sqrt(4*(J2+0.5)**2+Ynu[0]*(Ynu[0]-4.0)))-Dnu[0]*J2**4
    K1=J1-0.5
    F1lower[J]=Bnu[0]*((J+1.0)**2-1-0.5*np.sqrt(4*(J+1)**2+Ynu[0]*(Ynu[0]-4.0)))-Dnu[0]*J**2*(J+1.0)**2
    K2=J2+0.5
    F2lower[J]=Bnu[0]*((K2**2-1.0)-1+0.5*np.sqrt(4*(K2)**2+Ynu[0]*(Ynu[0]-4.0)))-Dnu[0]*K2**2*(K2+1.0)**2

for J in range(0,Jnum):
    J1=(J*2+2)/2.0
#    F1upper[0,J]=Bnu[1]*((J1+0.5)**2-1-0.5*np.sqrt(4*(J1+0.5)**2+Ynu[1]*(Ynu[1]-4.0)))-Dnu[1]*J1**4
    J2=(J*2+1)/2.0
#    F2upper[0,J]=Bnu[1]*((J2+0.5)**2-1+0.5*np.sqrt(4*(J2+0.5)**2+Ynu[1]*(Ynu[1]-4.0)))-Dnu[1]*J2**4
    K1=J1-0.5
    F1upper[J]=Bnu[1]*((J+1.0)**2-1-0.5*np.sqrt(4*(J+1)**2+Ynu[1]*(Ynu[1]-4.0)))-Dnu[1]*J**2*(J+1.0)**2
    K2=J2+0.5
    F2upper[J]=Bnu[1]*((K2**2-1.0)-1+0.5*np.sqrt(4*(K2)**2+Ynu[1]*(Ynu[1]-4.0)))-Dnu[1]*K2**2*(K2+1.0)**2
    
    Fdouble[0,J]=qnu[0]*K1*(K1+1)   #3/2 state splitting
    Fdouble[1,J]=qnu[1]*K1*(K1+1)   #1/2 state splitting

DG=Gnu[1]-Gnu[0]


Erel[0,:]=F1upper[:]-F1lower[:]+DG
Erel[1,:]=F2upper[:]-F2lower[:]+DG

for i in range(0,2):
    for J in range(0,Jnum-1):
        
        Pb[0,J+(Jnum-1)*i]=F1upper[J]-F1lower[J+1]+DG+(-1)**i*Fdouble[0,J]
        Rb[0,J+(Jnum-1)*i]=F1upper[J+1]-F1lower[J]+DG+(-1)**i*Fdouble[0,J]
        
        Pb[1,J+(Jnum-1)*i]=F2upper[J]-F2lower[J+1]+DG+(-1)**i*Fdouble[0,J]
        Rb[1,J+(Jnum-1)*i]=F2upper[J+1]-F2lower[J]+DG+(-1)**i*Fdouble[0,J]
#    for J in range(0,Jnum):
        Qb[1,J+(Jnum-1)*i]=F1upper[J]-F1lower[J]+DG+(-1)**i*Fdouble[0,J]
        Qb[0,J+(Jnum-1)*i]=F2upper[J]-F2lower[J]+DG+(-1)**i*Fdouble[0,J]
        
        
        Pcross[0,J+(Jnum-1)*i]=F1upper[J]-F2lower[J+1]+DG+(-1)**i*Fdouble[0,J]
        
    
    
    
Wdouble=Fdouble+Erel[0,:]
WP[0,:]=1/Pb[0,:]*10**7
WR[0,:]=1/Rb[0,:]*10**7

WQ[0,:]=1/Qb[0,:]*10**7
WQ[1,:]=1/Qb[1,:]*10**7
WP[1,:]=1/Pb[1,:]*10**7
WR[1,:]=1/Rb[1,:]*10**7


Wcross[0,:]=1/Pcross[0,:]*10**7


Wavelength=1/Erel*10**7
Wdouble=1/Wdouble*10**7

for J in range(0,Jnum):
    J1=(J*2+1)/2.0+1.0
    J2=(J*2+1)/2.0
    Intensity[0,J]=(2.0*J1+1)*np.exp(-1.44*F1upper[J]/Temp)
    Intensity[1,J]=(2.0*J2+1)*np.exp(-1.44*F2upper[J]/Temp)
#    Intdouble[J]=(2.0*J1+1)*np.exp(-1.44*Fdouble[J]/Temp)
    Q[0,J]=(2.0*J1+1)*np.exp((-1.44*(F1upper[J])/Temp))
    Q[1,J]=(2.0*J1+1)*np.exp((-1.44*(F2upper[J])/Temp))


Part[0]=np.sum(Q[0,:])  #Partition function of excited states
Part[1]=np.sum(Q[1,:])  #Partition function of excited states

for i in range(0,2):
    for J in range(0,Jnum-1):
        J1=(J*2+1)/2.0+1.0
        J2=(J*2+1)/2.0
        
        Pint[0,J+(Jnum-1)*i]=(2.0*J1+1)*(np.exp((-1.44*(F1upper[J])/Temp)))/Part[0] #Pij[J]*
        Rint[0,J+(Jnum-1)*i]=(2.0*J1+2)*np.exp((-1.44*F1upper[J+1])/Temp)/Part[0]
        
        Pint[1,J+(Jnum-1)*i]=(2.0*J2+1)*np.exp((-1.44*F2upper[J])/Temp)/Part[1]
        Rint[1,J+(Jnum-1)*i]=(2.0*J2+2)*np.exp((-1.44*F2upper[J+1])/Temp)/Part[1]
        
#    for J in range(0,Jnum):
        Qint[0,J+(Jnum-1)*i]=(2.0*J1+1)*np.exp((-1.44*F1upper[J])/Temp)/Part[0]
        Qint[1,J+(Jnum-1)*i]=(2.0*J2+1)*np.exp((-1.44*F2upper[J])/Temp)/Part[1]
        

        Cint[0,J+(Jnum-1)*i]=(2.0*J1+1)*np.exp((-1.44*F1upper[J]+(-1)**i*Fdouble[0,J])/Temp)
    
inttot=1.0#np.sum(Pint[:,:])+np.sum(Rint[:,:])+np.sum(Qint[:,:])+np.sum(Cint[0,:])

plt.figure()
#plt.scatter(Wavelength[0,:],Intensity[0,:]/inttot)
#plt.scatter(Wavelength[1,:],Intensity[1,:]/inttot)
#plt.scatter(Wdouble[:],Intdouble[:]/inttot)
plt.stem(WP[0,:],Pint[0,:]/inttot,label='P(3/2)',linefmt='-') 
plt.stem(WR[0,:],Rint[0,:]/inttot,label='R(3/2)',linefmt='--')
plt.stem(WP[1,:],Pint[1,:]/inttot,label='P(1/2)',linefmt='-')
plt.stem(WR[1,:],Rint[1,:]/inttot,label='R(1/2)',linefmt='--')
plt.stem(WQ[1,:],Qint[1,:]/inttot,label='Q(3/2)',linefmt=':')
plt.stem(WQ[1,:],Qint[1,:]/inttot,label='Q(1/2)',linefmt=':')

#plt.stem(Wcross[0,:],Cint[0,:]/inttot,label='Pcross',linefmt='-.')

#plt.xlim(1500,1560)

plt.legend()













