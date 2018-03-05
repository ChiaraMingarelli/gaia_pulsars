#KGrover
# Calculates differentce between frequency at earth adn at 1kpc in past at 2PN order and a 
# range of masses and earth frequencies. Produces 3 text files for the mass range, frequency range and change in frequency.
# (didn't plot in same file as takes a while to run so easier to use plot_PN.py seperatley) 

#!/usr/bin/env python

from math import *
import numpy as np
import scipy.integrate as integrate

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

def calcBeta(m1, m2, eta, Lchi1, Lchi2):
	return (1./12.)*((113.*(m1/(m1+m2))**2. + 75.*eta)*Lchi1 + (113.*(m2/(m1+m2))**2. + 75.*eta)*Lchi2 )
	
def calcSigma(eta, chi1chi2, Lchi1, Lchi2):
	return (eta/48.)*(-247.*chi1chi2+721.*Lchi1*Lchi2)

#calculate the frequency at a given time from observed earth frequency
def calcFreq(t, mchirp, fE): 
	return (-(256./5.)*pi**(8./3.)*mchirp**(5./3.)*t+fE**(-8./3.))**(-3./8.)

def dfdt_PN(f, t0, mchirp, eta,  beta, sigma):	
	b0=1.
	b1=0.
	b2=-((743./336.)+(11./4.)*eta)*(pi*mchirp*f)**(2./3.)*eta**(-2./5.)
	b3=(4*pi-beta)*(pi*mchirp*f)**(3./3.)*eta**(-3./5.)
	b4=((34103./18144.)+(13661./2016.)*eta+(59./18.)*eta*eta+sigma)*(pi*mchirp*f)**(4./3.)*eta**(-4./5.)
	return (96./5.)*pi**(8./3.)*mchirp**(5./3.)*f**(11./3.)*(b0+b1+b2+b3+b4)

def calcFreqPN(fE, t_array, mchirp, eta, beta, sigma):
	fp_PN_array=integrate.odeint(dfdt_PN,fE,t_array,args=(mchirp,eta,beta,sigma))
	return fp_PN_array
	
	
def createPlot(q, fignum):
	t_size=50
	lookback=1  #kpc
	lookback = -kpcToSec(lookback)

	steps=1000
	fE_min=1e-9
	fE_max=1e-7
	m1_min=1e8
	m1_max=1e9
	
	chi1=1.
	chi2=1.
	cosLchi1=1.
	cosLchi2=1.
	coschi1chi2=1.
	Lchi1=chi1*cosLchi1
	Lchi2=chi1*cosLchi2
	chi1chi2=chi1*chi2*coschi1chi2
	
	fE_array=np.empty(steps)
	log_fE_array=np.empty(steps)
	m1_solar_array=np.empty(steps)
	log_m1_solar_array=np.empty(steps)
	m1_array=np.empty(steps)
	mchirp_array=np.empty(steps)
	eta_array=np.empty(steps)
	beta_array=np.empty(steps)
	sigma_array=np.empty(steps)
	
	t_array=np.empty(t_size)
	
	for k in range(t_size):
		t_array[k]=k*(lookback)/(t_size-1.)
	
	for k in range(steps):
		fE_array[k]=fE_min +k*((fE_max-fE_min)/(steps-1.))
		log_fE_array[k]=log(fE_array[k],10)
		m1_solar_array[k]=m1_min +k*((m1_max-m1_min)/(steps-1.))
		log_m1_solar_array[k]=log(m1_solar_array[k],10)
		m1_array[k]=solarMassToSec(m1_solar_array[k])
		mchirp_array[k]=calcChirpMass(m1_array[k], solarMassToSec(m1_solar_array[k])/q)
		eta_array[k]=calcEta(m1_array[k], m1_array[k]/q)
		beta_array[k]=calcBeta(m1_array[k], m1_array[k]/q, eta_array[k], Lchi1, Lchi2)
		sigma_array[k]=calcSigma(eta_array[k], chi1chi2, Lchi1, Lchi2)
		
	deltaf_array=np.empty([steps,steps])
	
	outfile1=open("deltaf_array.txt", 'w')
	outfile2=open("log_m1_solar_array.txt", 'w')
	outfile3=open("log_fE_array.txt", 'w')
	for k in range(steps):
		outfile2.write(str(log_m1_solar_array[k])+",")
		outfile3.write(str(log_fE_array[k])+",")
		print k
		for j in range(steps):
			deltaf_array[k][j]=abs(calcFreq(lookback, mchirp_array[k], fE_array[j]) - calcFreqPN(fE_array[j], t_array, mchirp_array[k], eta_array[k], beta_array[k], sigma_array[k])[t_size-1])
			outfile1.write(str(deltaf_array[k][j])+",")
		outfile1.write("\n")
	
	outfile1.close()
	outfile2.close()
	outfile3.close()
	
if __name__=="__main__":
	createPlot(1.,1)