#K Grover
# Calculates differentce between frequency at earth adn at 1kpc in past for newtonain order and a 
# range of masses and earth frequencies. Produces a contour plot.

#!/usr/bin/env python

from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, ticker
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

#plot the results
def plotDeltaf(log_m1_solar_array, log_fE_array, diff_array,q, fignum):
	V=([3.17e-9])
	fmtdict={3.17e-9:"fmin"}
	plt.figure(fignum)
	plt.contourf(log_m1_solar_array,log_fE_array,diff_array, locator=ticker.LogLocator()) 
	plt.colorbar()
	plt.grid(True, which='both')
	CS = plt.contour(log_m1_solar_array,log_fE_array,diff_array, V,
			colors='k', 
		   )
	plt.clabel(CS, fmt=fmtdict)
	plt.xlabel(r"$log(m_1/M_\odot)$")
	plt.ylabel(r"$log(f_E/Hz)$")
	plt.savefig("deltaf_contour_q"+str(q)+".png")

#calc difference and plot for a given mass ratio
def createPlot(q, fignum):
	lookback= 1  #kpc
	lookback = -kpcToSec(lookback)

	steps=500
	fE_min=1e-9
	fE_max=1e-7
	m1_min=1e8
	m1_max=1e9
	
	fE_array=np.empty(steps)
	log_fE_array=np.empty(steps)
	m1_solar_array=np.empty(steps)
	log_m1_solar_array=np.empty(steps)
	m1_array=np.empty(steps)
	mchirp_array=np.empty(steps)
	
	for k in range(steps):
		fE_array[k]=fE_min +k*((fE_max-fE_min)/(steps-1.))
		log_fE_array[k]=log(fE_array[k],10)
		m1_solar_array[k]=m1_min +k*((m1_max-m1_min)/(steps-1.))
		log_m1_solar_array[k]=log(m1_solar_array[k],10)
		m1_array[k]=solarMassToSec(m1_solar_array[k])
		mchirp_array[k]=calcChirpMass(m1_array[k], solarMassToSec(m1_solar_array[k])/q)
		
	deltaf_array=np.empty([steps,steps])
	
	for k in range(steps):
		for j in range(steps):
			deltaf_array[k][j]=fE_array[j]-calcFreq(lookback, mchirp_array[k], fE_array[j])
			
	plotDeltaf(log_m1_solar_array, log_fE_array, deltaf_array,q, fignum)

#To calculate specific frequency differences 
def specificfreqdiff(m1, fE, lockback):
	mchirp=calcChirpMass(solarMassToSec(m1),solarMassToSec(m1))
	return fE-calcFreq(lookback, mchirp, fE)

if __name__=="__main__":
	createPlot(1,1)
	createPlot(3,2)

#Calculate the difference between fE and fpulsar at 9 different mass freq combos
	#lookback=1  #kpc
	#lookback = -kpcToSec(lookback)
	#print "M=1e9	fe=1e-7	", specificfreqdiff(1e9, 1e-7, lookback)
	#print "M=1e9	fe=5e-8	", specificfreqdiff(1e9, 5e-8, lookback)
	#print "M=5e8	fe=5e-8	", specificfreqdiff(5e8, 5e-8, lookback)

	

	