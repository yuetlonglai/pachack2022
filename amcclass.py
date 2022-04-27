import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import signal
from scipy import stats

'''Note that this programme runs on the expectation that the reflectance of the data is normalised'''
'''This is because we discovered that reflectance should be given as a ratio => it should already be a normalised quantity'''

class asteroid:
    '''This is a class used to identify an asteroid's group by looking at the reflectance spectrum'''
    def __init__(self,filename):
        self._filename = filename
        self._wl,self._R,self._ep=np.loadtxt(filename,skiprows=1,unpack=True)

        self._wlm=[]
        self._Rm=[]
        self._epm=[]
        for i in range(len(self._wl)):
            if self._wl[i]<925 and self._wl[i]>435:
                self._wlm.append(self._wl[i])
                self._Rm.append(self._R[i])
                self._epm.append(self._ep[i])
        self._epm = np.array(self._epm)

        #albedo
        self._al = np.mean(self._Rm)
        self._alerr = np.sqrt(sum(self._epm**2))

        #quadratic fit
        self._popt2,self._pcov2=np.polyfit(self._wlm,self._Rm,2,cov=True)
        self._f2=np.poly1d(self._popt2)

        #linear fit
        self._popt1,self._pcov1=np.polyfit(self._wlm,self._Rm,1,cov=True)
        self._f1=np.poly1d(self._popt1)

    def albedo(self):
        '''print albedo'''
        print("Albedo = %.3e +/- %.3e" %(self._al,self._alerr))
        return self

    def showplot(self):
        '''plot reflectance spectrum'''
        self._peaks,_=signal.find_peaks(self._Rm)

        self._pwl = self._wl[np.argmax(self._Rm)]
        #print('Location of peak = %.3e nm ' %(self._pwl))
                
        plt.plot(self._wl[self._peaks],self._R[self._peaks],'.',color='red')
        #plt.errorbar(wl,R,yerr=ep,fmt='-',color='black')
        plt.plot(self._wl,self._R,label=self._filename)
        return self
    
    def quadfit(self):
        '''plot the quadratic fit within the visible region to see if it is an S group or not'''
        peak=(-self._popt2[1])/(2*self._popt2[0])
        plt.plot(self._wlm,self._f2(self._wlm),label=str(self._filename) + ' quadratic fit')
        print('Location of peak = %.3e nm ' %(peak))
        return self
        
    
    def linearfit(self):
        '''plot the linear fit within the visible region to see if it is an S group or not'''
        plt.plot(self._wlm,self._f1(self._wlm),label=str(self._filename) + ' linear fit')
        return self

    def normal_dist(x , mean , sd):
        prob_density = (np.pi*sd) * np.exp(-0.5*((x-mean)/sd)**2)
        return prob_density

    def classify(self):
        '''method used to make a classification verdict'''
        corr1=stats.chisquare(self._Rm,self._f1(self._wlm))
        corr2=stats.chisquare(self._Rm,self._f2(self._wlm))
        #print(np.gradient(self._f2(self._wlm))[-1])
        print('Chi-squared value for linear fit = ' +str(corr1[0]))
        print('Chi-squared value for quadratic fit = ' + str(corr2[0]))
        #print(corr1[1],corr2[1])

        secondbound=0.30
        firstbound=0.10
        sbound=-0.001 #gradient of the quadratic fit beyond which is negative enough to say it is having a reddish decline
        
        self._clas=''
        #if corr2<corr1 and np.gradient(self._f2(self._wlm))[-1] < sbound and self._al <=secondbound:
        if corr2[0]/corr1[0] <= 0.2 and self._al >= firstbound and self._al <= secondbound: #see if the chi squared value differ a lot which makes the S group classification verdict certain
            self._clas='S'
            print('Classification = ' + str(self._clas))
        elif corr2[0]<corr1[0] and np.gradient(self._f2(self._wlm))[-1] >= sbound or corr1[0]<corr2[0]:
            if self._al <= secondbound and self._al >= firstbound:
                self._clas='Xm'
                print(np.gradient(self._f2(self._wlm))[-1])
                print('Classification = ' + str(self._clas))
            elif self._al < firstbound:
                #print('Classification = C or Xp')
                if corr1[0]<0.01:
                    self._clas='Xp'
                    print('Classification = ' + str(self._clas))
                else:
                    self._clas='C'
                    print('Classification = ' + str(self._clas))
            elif self._al > secondbound:
                self._clas='Xe'
                print('Classification = ' + str(self._clas))
        else:
            self._clas='Uncertain'
            print('Classification = ' + str(self._clas))

        self._Prob=0

        if self._clas == 'S':
            if corr2[0]<0.016:
                prob1=0.90
            else:
                prob1=0.10

            cdf_ul=stats.norm(loc=0.26,scale=0.09).cdf(self._al+self._alerr)
            cdf_ll=stats.norm(loc=0.26,scale=0.09).cdf(self._al-self._alerr)
            prob2=cdf_ul-cdf_ll

            self._Prob=prob1*prob2

        if self._clas == 'Xm':
            if corr1[0]<0.016:
                prob1=0.90
            else:
                prob1=0.1

            cdf_ul=stats.norm(loc=0.14,scale=0.05).cdf(self._al+self._alerr)
            cdf_ll=stats.norm(loc=0.14,scale=0.05).cdf(self._al-self._alerr)
            prob2=cdf_ul-cdf_ll

            self._Prob=prob1*prob2

        if self._clas == 'Xp':
            if corr1[0]<0.016:
                prob1=0.9
            else:
                prob1=0.1

            cdf_ul=stats.norm(loc=0.05,scale=0.01).cdf(self._al+self._alerr)
            cdf_ll=stats.norm(loc=0.05,scale=0.01).cdf(self._al-self._alerr)
            prob2=cdf_ul-cdf_ll

            self._Prob=prob1*prob2
        if self._clas == 'C':
            if corr1[0]<0.016:
                prob1=0.9
            else:
                prob1=0.1
            
            cdf_ul=stats.norm(loc=0.08,scale=0.08).cdf(self._al+self._alerr)
            cdf_ll=stats.norm(loc=0.08,scale=0.08).cdf(self._al-self._alerr)
            prob2=cdf_ul-cdf_ll

            self._Prob=prob1*prob2
        if self._clas == 'Xe':
            if corr1[0]<0.016:
                prob1=0.9
            else:
                prob1=0.1
             
            cdf_ul=stats.norm(loc=0.54,scale=0.25).cdf(self._al+self._alerr)
            cdf_ll=stats.norm(loc=0.54,scale=0.25).cdf(self._al-self._alerr)
            print(self._al+self._alerr, self._al-self._alerr)
            prob2=cdf_ul-cdf_ll

            self._Prob=prob1*prob2
        
        print('Classification = ' + str(self._clas) + ', with certainty of %.3f  ' %(self._Prob*100))
        print(self._clas,self._Prob,prob1,prob2)

        return self


        