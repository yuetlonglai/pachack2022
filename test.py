import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
from scipy import signal
from scipy import optimize

filename='cacr21.txt'

def cubic(x,a,b,c):
    return a*(x-b)**2 +c

wl,R,ep=np.loadtxt(filename,skiprows=1,unpack=True)

wlm=[]
Rm=[]
for i in range(len(wl)):
    if wl[i]<925 and wl[i]>435:
        wlm.append(wl[i])
        Rm.append(R[i])
albedo=np.mean(Rm)
albedoerr=np.std(Rm)

peaks,_=signal.find_peaks(Rm)

pwl = wl[np.argmax(Rm)]
#print('Location of peak = %.3e nm ' %(pwl))
#popt,pcov=optimize.curve_fit(cubic,wlm,Rm,p0=[-1e-4,700,0])
#print(popt)
popt1,pcov1=np.polyfit(wlm,Rm,2,cov=True)
f1=np.poly1d(popt1)

plt.figure()
plt.grid()
plt.title(filename)
plt.ylabel("Normalised Reflectance")
plt.xlabel("Wavelength (nm)")
plt.plot(wl[peaks],R[peaks],'.',color='red')
#plt.errorbar(wl,R,yerr=ep,fmt='-',color='black')
plt.plot(wl,R,color='blue')
plt.plot(wlm,f1(wlm),color='cyan')
plt.xlim(300,1000)
plt.show()


plt.figure()
plt.title('')
plt.plot(wl,np.gradient(R))
plt.show()




#print("Albedo = %.3e +/- %.3e" %(albedo,albedoerr))
