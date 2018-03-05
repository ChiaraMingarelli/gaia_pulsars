#Chiara's edited K Grover edited version of Chiara's code.
# Used to produce alpha plots for the paper
#also used to make S/L plots.

#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from matplotlib import rc, text, ticker
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
rc('text', usetex=True)
#plt.rcParams['font.size'] = 20

### CONSTANTS ### will use geometric units, all in seconds.
c=2.99792458*(10**8)
G=6.67428*(10**(-11))
s_mass=G*(1.98892*10**(30))/(c**3)

def mu(m1,m2): return s_mass*(m1*m2)/(m1+m2) 

def M(m1,m2): return s_mass*(m1+m2)

def mchirp(m1,m2): return ((mu(m1,m2))**(3./5))*((M(m1,m2))**(2./5))

def parsec(d): return d*3.08568025e16/299792458   #converts parsecs to secs

def time(t): return t*299792458/3.08568025e16     #converts light-seconds to parsecs

def years(t): return t*31556926.			#converts years to seconds

#Calculate time from coalesence at a given frequency 
def tf(m1,m2,f): return -5.*((8*math.pi*f)**(-8./3))*mchirp(m1,m2)**(-5./3)

#calculate frequency at a time from coalesence, using above (note time given should be negative back in time from coalesence)
def freq(t, m1,m2): return ((1./5.)*(-t)*mchirp(m1,m2)**(5./3.))**(-3./8.)*(1./(8.*math.pi)) 


# S/L calculation to see magnitude for simple precession calculation (for a given frequency)
def SoverL_freq(f,a,m1,m2):
	ans=(a/2.)*(8*math.pi*f)**(1./3)*mchirp(m1,m2)**(5./24)*(M(m1,m2)/mu(m1,m2))**(5./4)
	return ans

# S/L calculation to see magnitude for simple precession calculation (for a given time before coalescence, t_c=0 time is negative)
def SoverL_time(t,a,m1,m2):
	ans=a*((M(m1,m2)/mu(m1,m2))**(5./4))*((-t)**(-1./8))*(256./5)**(-1./8)
	#print ans
	return ans

def thom_pre(t,a,m1,m2):
	ans= 2*math.pi(1-math.cos(SoverL_time(t,a,m1,m2)))
	return ans
	

if __name__=="__main__":
#masses of the system
	m1=1e9
	m2=1e9
#frequency observed on earth
	f_earth=1e-7
	
#the time parameters in seconds all relative to tc=0
	t_earth= tf(m1,m2,f_earth)
	t_1kpc=t_earth-parsec(1e3)
	t_2kpc= t_earth-parsec(2e3)
	t_max= t_earth-parsec(8e6)
	t_10yrs= t_earth+years(10)

#create arrays
	steps=9000000
	lookback_step=(-t_max)/(steps-1.)
#goes from 0 to t_max (0 to negative numbers)
	lookback_array=np.empty([steps])
#lookback array but positive and in years, i.e. years from coalesence
	lookback_yrs_array=np.empty([steps])
#S/L vector
	SoverL=np.empty([steps])
	thom_p=np.empty([steps])
#fill arrays
	for k in range(steps):
		lookback_array[k]= t_max+k*lookback_step-1
		#print lookback_array[k]
		lookback_yrs_array[k]= -lookback_array[k]*3.171e-8
		#print lookback_yrs_array[k]
		SoverL[k]=SoverL_time(lookback_array[k],1.0,m1,m2)

	change=tf(m1,m2,1e-9)-tf(m1,m2,1e-7)
	#print change*3.171e-8

	#print "parsecs in 26million years", time(years(2.6e7))

	plt.figure(1)
	plt.loglog(lookback_yrs_array, SoverL,c="m", label=r'$S/L(t)$, $a=1.0$')
	plt.loglog(lookback_yrs_array, thom_p,c='c', label="Thom.prec. equiv.")
	plt.ylabel(r"$S/L(t)$")
	plt.xlabel(r"Time from coalesence (years)")
	plt.axvline(x=-t_earth/31556926., c="b",linestyle='-',label=r"$t_E,10^{-7}$ Hz")
	plt.axvline(x=-t_1kpc/31556926.,c="r",linestyle='-',label="$1$ kpc")
	plt.axvline(x=-t_10yrs/31556926., c="g",linestyle='-', label="$10$ yrs")
	plt.legend(loc="upper right",fancybox="True")
	plt.ylim(0,0.35)
	plt.grid(True)
	plt.hold(True)
	plt.xlim(10,1e7)
	plt.savefig("SoverL_8Mpc_thomas.png")

	

