#K Grover edited version of Chiara's code.
# Used to produce alpha plots for the paper


#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from matplotlib import rc, text, ticker
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
rc('text', usetex=True)
plt.rcParams['font.size'] = 20

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
    
print "freq for 10^5yrs, 10^8Msunbinary,\n", freq(years(-10**5), 10**8,10**8)
print "freq for 10^5yrs, 10^10Msunbinary,\n", freq(years(-10**5), 10**10,10**10)
print "time to coal for 10**8M binary at 30 nHz,\n", tf(10**8,10**8, 3e-8)/31556926
print "time to coal for 10**10M binary at 2 nHz,\n", tf(10**10,10**10, 2e-9)/31556926

#an approximate alpha for L>>S
def alpha(f,lookback,m1,m2):
    t2=0.
    t1=lookback
    alpha=(-8./3)*((256./5)**(-5./8))*(2.+((3.*m2)/(2*m1)))*(mu(m1,m2)**(3./8))*(M(m1,m2)**(-3./4))*(((-t2)**(3./8))-((-t1)**(3./8)))
    return alpha

#A full expression for d_alpha/dt not making any assumptions on L or S.
def d_alpha(t,a,m1,m2,Q): 
    c1=mu(m1,m2)**(3./8)*M(m1,m2)**(-3./4)*((256./5)**(-5./8))*(2.+1.5*m2/m1)
    c2=2.*Q*a*mu(m1,m2)**(-9./8)*M(m1,m2)**(5./4)*((256./5)**(-1./8))
    c3=a**2*mu(m1,m2)**(-9./4)*M(m1,m2)**(5./2)*((256./5)**(-1./4))
    d_alpha=c1*(-t)**(-5./8)*math.sqrt(1+c2*(-t)**(-1./8)+c3*(-t)**(-1./4))
    return d_alpha

# the integrated d_alpha between two times.
def big_alpha(f,lookback,a,m1,m2,Q):
    t2=-1.
    t1=lookback
    alpha=quad(d_alpha,t1,t2,args=(a,m1,m2,Q))
    return alpha[0]

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
	t_max= t_earth-parsec(3e3)
	t_10yrs= t_earth+years(10)

#create arrays
	steps=2000
	lookback_step=(-t_max)/(steps-1.)
#goes from 0 to t_max (0 to negative numbers)
	lookback_array=np.empty([steps])
#lookback array but positive and in years, i.e. years from coalesence
	lookback_yrs_array=np.empty([steps])
#Different alphas as function of time
	alpha_approx=np.empty([steps])
	alpha_exact_large_spin=np.empty([steps])
	alpha_exact_small_spin=np.empty([steps])
	
#fill arrays
	for k in range(steps):
		lookback_array[k]= t_max+k*lookback_step
		lookback_yrs_array[k]= -lookback_array[k]/31556926.
		alpha_approx[k]=alpha(f_earth, lookback_array[k], m1, m2)
		alpha_exact_large_spin[k]=big_alpha(f_earth,lookback_array[k],0.98,m1,m2,1.)
		alpha_exact_small_spin[k]=big_alpha(f_earth,lookback_array[k],0.1,m1,m2,1.)

#renormalise alpha so Dalpha = 0 at earth
	for k in range(steps):
		alpha_approx[k]=alpha_approx[k]- alpha(f_earth, t_earth, m1, m2)
		alpha_exact_large_spin[k]=alpha_exact_large_spin[k]-big_alpha(f_earth,t_earth,0.98,m1,m2,1.)
		alpha_exact_small_spin[k]=alpha_exact_small_spin[k]-big_alpha(f_earth,t_earth,0.1,m1,m2,1.)
	

	#plot alpha(t) on a semilog plot
	fig=plt.figure(1)
	ax=fig.add_subplot(111)
	majorFormatter = FormatStrFormatter('%d')
	plt.semilogx(lookback_yrs_array, alpha_exact_small_spin, 'k--', label=r'$a=0.10$')	
	plt.semilogx(lookback_yrs_array, alpha_exact_large_spin,'k', label=r'$a=0.98$')
	plt.ylabel(r"$\alpha-\alpha_{Earth}$ (radians)")
	plt.xlabel(r"Time from coalesence (years)")
	plt.axvline(x=-t_earth/31556926., c="b",linestyle='-')
	plt.axvline(x=-t_1kpc/31556926.,c="r",linestyle='-')
	plt.axvline(x=-t_10yrs/31556926., c="g",linestyle='-')
	plt.legend(loc="upper left",fancybox="True")
#	plt.grid(True, which="both")
	plt.hold(True)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.ylim(-1e2,4e2)
	plt.xlim(10,1e4)
#	plt.text(-t_earth/31556926. ,410,r"$t_E$", rotation=0., size=15)
#	plt.text(-t_1kpc/31556926.-10,420,r"$t_E- 1kpc$", rotation=80., size=10)
#	plt.text(-t_10yrs/31556926.-10,420,r"$t_E +10yrs$", rotation=80., size=10)
#	ax.xaxis.set_major_formatter(majorFormatter)
	ax.yaxis.set_major_formatter(majorFormatter)
	plt.savefig("semilog_alpha_cm.png")

	print "Earth is at",- t_earth/31556926.
	print "+10 years at",-t_10yrs/ 31556926.

	print "For a = 0.1"
	print "delta alpha at Earth + 10 years alpha = ", big_alpha(f_earth, t_10yrs,0.1,m1,m2,1.)-big_alpha(f_earth,t_earth,0.1,m1,m2,1.)
	print "delta alpha at Eath - 1kpc alpha =" , big_alpha(f_earth,t_earth-parsec(1e3),0.1,m1,m2,1.)-big_alpha(f_earth,t_earth,0.1,m1,m2,1.)

	print "For a = 0.98"
	print "delta alpha at Earth + 10 years alpha = ", big_alpha(f_earth,t_10yrs,0.98,m1,m2,1.)-big_alpha(f_earth,t_earth,0.98,m1,m2,1.)
	print "delta alpha at Eath - 1kpc alpha =" , big_alpha(f_earth,t_earth-parsec(1e3),0.98,m1,m2,1.)-big_alpha(f_earth,t_earth,0.98,m1,m2,1.)
