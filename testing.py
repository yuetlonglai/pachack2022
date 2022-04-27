import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import signal
from amcclass import asteroid

#stype1=asteroid('c1mt117a.tab')
stype2=asteroid('c1dd01.txt') #S
stype3=asteroid('c1dd09.txt') #S
stype4=asteroid('c1dd109.txt') #S
stype5=asteroid('c1tb56.tab') #Xm
stype6=asteroid('c1dd65.txt') #Xp
stype7=asteroid('c1dd66.txt') #Xp
stype8=asteroid('c1ag11.txt') #Xe
stype9=asteroid('c1dd85.txt') #Xe
stype10=asteroid('cacr21.txt') #C
stype11=asteroid('cacr22.txt') #C
stype12=asteroid('c1er34.txt') #Xp

#stype3.classify()

plt.figure()
plt.grid()
plt.title('Reflectance Spectrum')
plt.ylabel("Normalised Reflectance")
plt.xlabel("Wavelength (nm)")
stype4.showplot().quadfit().linearfit().albedo().classify()
#stype7.albedo()
#plt.xlim(300,1000)
plt.legend()
plt.show()

#stypelist=[stype2,stype3,stype4,stype5,stype6,stype7,stype8,stype9,stype10,stype11,stype12]
#for i in stypelist:
#    print(i.classify())