#K.Grover
# Plots the frequecny evolution of a SMBHB for a hard coded masses and an observed Earth frequency

#!/usr/bin/env python

from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)
plt.rcParams['font.size'] = 20

#converts solar mass to seconds
def solarMassToSec(m): return m*4.925*pow(10.,-6.)

#converts kpc to light seconds
def kpcToSec(d): return d*1.02927214e11

#calculates eta for two masses
def calcEta(m_1, m_2): return (m_1*m_2)/((m_1+m_2)*(m_1+m_2))

#calculates chirp mass for two masses
def calcChirpMass(m_1,m_2): return pow((m_1*m_2)/((m_1+m_2)*(m_1+m_2)),3./5.)*(m_1+m_2)
	
#calculates the time to coalesence if the time at the earth is 0
def calctc(fE, mchirp): return 5*(8*pi*fE)**(-8./3.)*mchirp**(-5./3.)

#calculate the frequency at a given time from observed earth frequency
def calcFreq(t, mchirp, fE): 
	return (-(256./5.)*pi**(8./3.)*mchirp**(5./3.)*t+fE**(-8./3.))**(-3./8.)	

#plot the difference in frequency between earth and some time t
def plotFreqDiff(t_array, freq_diff_array, t_max):

#Renormalise to shift tc to 0  so array is seconds from coalesence, then change to years
	renorm_t_array=np.empty_like(t_array)
	for k in range(t_array.size):
		renorm_t_array[k]=-t_array[k]+t_max
		renorm_t_array[k]=renorm_t_array[k]*3.16887646e-8 
	
	plt.figure(1)
	plt.semilogx(renorm_t_array, freq_diff_array)
	plt.grid(True, which='both')
	plt.xlabel("Time from coalesence /years")
	plt.ylabel(r"$f_E-f$  /Hz")
	plt.axvline(x=t_max*3.16887646e-8 , c='r',linestyle="--", label="earth")
	plt.axvline(x=((kpcToSec(1)+t_max)*3.16887646e-8), c='c' ,linestyle="--",label="earth-1kpc")
	plt.axvline(x=(t_max*3.16887646e-8-10), c='g',linestyle="--", label="earth+10yrs")
	plt.axhline(y=(3.171e-9), c='y', linestyle="-", label="fE+fmin") 
	plt.axhline(y=(-3.171e-9), c='y', linestyle="-", label="fE-fmin") 
	#plt.axhline(y=(1.653e-6), c='g', label="1/1week")
	plt.legend(loc="lower right")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))	
	plt.xlim(1,1e4)
	plt.savefig("freq_evolution.png")
	
#plot zoomed in version around the earth and + 10 years
	plt.figure(2)
	plt.semilogx(renorm_t_array, freq_diff_array)
	plt.grid(True, which='both')
	plt.xlabel("Time from coalesence /years")
	plt.ylabel(r"$f_E-f$  /Hz")
	plt.xlim((t_max*3.16887646e-8-10)-1.,t_max*3.16887646e-8+1.)
	plt.ylim(-1e-8,1e-8)
	plt.axvline(x=t_max*3.16887646e-8 , c='r', linestyle="--",label="earth")
	#plt.axvline(x=((kpcToSec(1)+t_max)*3.16887646e-8), c='c' ,label="earth-1kpc")
	plt.axvline(x=(t_max*3.16887646e-8-10), c='g', linestyle="--",label="earth+10yrs")
	plt.axhline(y=(3.171e-9), c='y', linestyle="-", label="fE+fmin") 
	plt.axhline(y=(-3.171e-9), c='y',  linestyle="-", label="fE-fmin") 
	#plt.axhline(y=(1.653e-6), c='g', label="1/1week")
	plt.legend(loc="lower right")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))	
	plt.savefig("freq_evolution_zoom.png")


if __name__=="__main__":
	
	m1=1e9			# in solar mass
	m2=1e9
	lookback=5.		# in kpc
	fE=1e-7			# in Hz
	
	#convert parameters to seconds
	m1=solarMassToSec(m1)
	m2=solarMassToSec(m2)
	lookback = kpcToSec(lookback)
	
	mchirp=calcChirpMass(m1,m2)
	
	steps=100000
	t_array=np.empty(steps)
	t_min=-lookback
	t_max=calctc(fE,mchirp)

	freq_diff_array=np.empty(steps)
	
	for k in range(steps):
		t_array[k]=t_min+((t_max-t_min)/(steps))*k
		freq_diff_array[k] = fE - calcFreq(t_array[k], mchirp, fE) 
	
	plotFreqDiff(t_array,freq_diff_array, t_max)
	
