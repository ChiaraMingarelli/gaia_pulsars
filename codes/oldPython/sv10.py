#!/usr/bin/env python

#version 1.1 alpha calc

import math
import matplotlib.pyplot as plt
from pylab import figure, axes, plot, xlabel, ylabel, title, grid, savefig, show
import numpy as np
from scipy.integrate import quad
from matplotlib import rc
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
rc('text', usetex=True)

### CONSTANTS ### will use geometric units, all in seconds.
c=2.99792458*(10**8)
G=6.67428*(10**(-11))
s_mass=G*(1.98892*10**(30))/(c**3)


def mu(m1,m2): return s_mass*(m1*m2)/(m1+m2) 

def M(m1,m2): return s_mass*(m1+m2)

def mchirp(m1,m2): return ((mu(m1,m2))**(3./5))*((M(m1,m2))**(2./5))

def parsec(d): return d*3.08568025e16/299792458   #converts parsecs to secs

def time(t): return t*299792458/3.08568025e16     #converts light-seconds to parsecs

def tf(m1,m2,f): return -5.*((8*math.pi*f)**(-8./3))*mchirp(m1,m2)**(-5./3)
print tf(1e9,1e9,1e-7)/3.1e7

def freq(tc,t, m1,m2): return ((1./5.)*(tc-t)*mchirp(m1,m2)**(5./3.))**(-3./8.)*(1./(8.*math.pi)) #gives frequency given a time in seconds  and mass in solar masses(note time given should be negative as assume tc=0)

def alpha(f,lookback,m1,m2):
    #t2=tf(m1,m2,f)
    t2=0.
    #print "t2 L>Sis %e", t2
    #t1=-lookback+t2
    t1=lookback
   # print "t1 is %e" %t1
    alpha=(-8./3)*((256./5)**(-5./8))*(2.+((3.*m2)/(2*m1)))*(mu(m1,m2)**(3./8))*(M(m1,m2)**(-3./4))*(((tc-t2)**(3./8))-((tc-t1)**(3./8)))
    return alpha

def alpha_S(f,lookback,m1,m2,a=1.):
    #t2=tf(m1,m2,f)
    t2=0.
    #print "t2, S<L %e", t2
    t1=lookback
    alpha=-4*a*((256./5)**(-3./4))*(2.+((3.*m2)/(2*m1)))*(mu(m1,m2)**(-3./4))*(M(m1,m2)**(1./4))*(((tc-t2)**(1./4))-((tc-t1)**(1./4)))
    return alpha

def al_alpha1(m1,m2,t1,t2,f): #The alpha(t) as calculated in Sesana Vecchio 2010.
    part1=2*(math.pi**(5./3))*(1.+((3.*m2)/(4.*m1)))*(mu(m1,m2))*(M(m1,m2)**(-1./3))*(f**(5./3))*(t2-t1)
    return (part1)

def al_alpha2(m1,m2,t1,t2,f): #The alpha(t) as calculated in Sesana Vecchio 2010 plus a 3rd order term 
    part1=2*(math.pi**(5./3))*(1.+((3.*m2)/(4.*m1)))*(mu(m1,m2))*(M(m1,m2)**(-1./3))*(f**(5./3))*(t2-t1)
    part2=32*(math.pi**(13./3))*(1.+((3.*m2)/(4.*m1)))*(mu(m1,m2)**(2))*(M(m1,m2)**(1./3))*(f**(13./3))*(t2-t1)*(t1+t2)
    return (part1+part2)
tc=0.
def d_alpha(t,a,m1,m2,Q): #A full expression for d_alpha/dt not making any assumptions on L or S.
    c1=mu(m1,m2)**(3./8)*M(m1,m2)**(-3./4)*((256./5)**(-5./8))*(2.+1.5*m2/m1)
    c2=2.*Q*a*mu(m1,m2)**(-9./8)*M(m1,m2)**(5./4)*((256./5)**(-1./8))
    c3=a**2*mu(m1,m2)**(-9./4)*M(m1,m2)**(5./2)*((256./5)**(-1./4))
    d_alpha=c1*(tc-t)**(-5./8)*math.sqrt(1+c2*(tc-t)**(-1./8)+c3*(tc-t)**(-1./4))
    return d_alpha

def big_alpha(f,lookback,a,m1,m2,Q):# the integrated d_alpha between two times.
    #t2=tf(m1,m2,f)
    t2=-1.
    #print "t2, BIG %e", t2
    #t1=-lookback+t2
    t1=lookback
    alpha=quad(d_alpha,t1,t2,args=(a,m1,m2,Q))
    return alpha[0]

def gamma(a,lookback,m1,m2):
    gamma=a*(M(m1,m2)**(5./4))*(mu(m1,m2)**(-9./8))*((256./5)**(-1./8))*(tc-lookback)**(-1./8)
    return gamma


#print "SV10 with 2 leading order terms says",  al_alpha1(1e9,1e9,0,3.1e8,1e-7)
#print "SV10 with 3 leading order terms says",  al_alpha2(1e9,1e9,0,3.1e8,1e-7)
#print "we say ", alpha(1e-7,-3.1e8,1e9,1e9)
#print "Ours-SV10 orignial:  %e" %(abs(alpha(1e-7,3.1e8,1e9,1e9)) - al_alpha1(1e9,1e9,0,3.1e8,1e-7))
#print "Ours-SV10 +3rd term:  %e" %(abs(alpha(1e-7,3.1e8,1e9,1e9)) - al_alpha2(1e9,1e9,0,3.1e8,1e-7))

#masses of the system
m1=1e9
m2=1e9
#frequency observed on earth
f_earth=1e-7

#the time parameters in seconds all relative to tc=0
tc=0.
t_earth=tc+tf(m1,m2,f_earth)
t_1kpc=tc+t_earth-parsec(1e3)
t_2kpc=tc+t_earth-parsec(2e3)
t_max=tc+t_earth-parsec(2.3e4)  #chosen so that the freq axis spans one order of magnitude
t_10yrs=tc+t_earth+10.*31556926.

#the frequency parameters in Hz (to draw vlines with on graph)
f_1kpc=freq(tc,t_1kpc,m1,m2)
f_2kpc=freq(tc,t_2kpc,m1,m2)
f_10yrs=freq(tc,t_10yrs, m1, m2)

#create arrays
steps=1000
lookback_step=(tc-t_max)/(steps-1.)

lookback_array=np.empty([steps])
lookback_yrs_array=np.empty([steps])
freq_array=np.empty([steps])
alpha_array=np.empty([steps])
alphaS_array=np.empty([steps])
LandS=np.empty([steps])
LandS_2=np.empty([steps])
gamma_array=np.empty([steps])
gamma_array2=np.empty([steps])

#test=M(m1,m2)**(5./4)*mu(m1,m2)**(-9./8)*(256./5)**(-1./8)
#print test


#fill arrays
for k in range(steps):
	lookback_array[k]= t_max+k*lookback_step
	lookback_yrs_array[k]= -lookback_array[k]*3.171e-8
	freq_array[k]=freq(tc,lookback_array[k],m1,m2)
	alpha_array[k]=alpha(f_earth, lookback_array[k], m1, m2)      # L>>S approximation
	#alphaS_array[k]=alpha_S(f_earth, lookback_array[k], m1, m2)
        LandS[k]=big_alpha(f_earth,lookback_array[k],0.01,m1,m2,1.)   #full solution with a=0.01
        LandS_2[k]=big_alpha(f_earth,lookback_array[k],0.99,m1,m2,1.) #full solution with a=0.99
        gamma_array[k]=gamma(0.99,lookback_array[k],m1,m2)
        gamma_array2[k]=gamma(0.01,lookback_array[k],m1,m2)
        
new=LandS_2-(LandS_2[0]- LandS[0])
#### TEST STATEMENTS #####

#print "zero value for gamma", gamma_array[0]
#print "gamma is", gamma(0.99,-3.12e12,1e9,1e9)
#print "alpha is", alpha(f_earth,-3.12e12,1e9,1e9)
#print "time is", (3.12e12/3.12e8)
#print "alpha for 1e9 at 1e-7 at 1kpc", alpha(f_earth,-parsec(1e3),1e9,1e9)
#print "Big alpha for 1e9 at 1e-7 at 1kpc", big_alpha(f_earth,-parsec(1e3),0.99,1e9,1e9,1.)
#print "new's last value",(LandS_2[0]- LandS[0])


#plot alpha(t) on a semilog/loglog plot


plt.figure(1)
#ax = self.figure.add_subplot(111)
ax1=plt.subplot(111)
plt.loglog(lookback_yrs_array, alpha_array,'b', label=r'$L\gg S$')
#plt.plot(lookback_yrs_array, LandS, 'r--', label=r'$L+S,a=0.01$')
plt.loglog(lookback_yrs_array, new, 'c', label=r'$L+S,a=0.99$')
#ax1.fill_between(lookback_yrs_array,new,alpha_array,where=alpha_array>=new,facecolor='0.89')
plt.xlabel(r"Time from coalesence (years)")
plt.xlim(lookback_yrs_array[steps-2]-1000,lookback_yrs_array[0])
plt.ylim(5,LandS[0])  #starts at 5 when using loglog plotting, can use next one when doing regular plots
#plt.ylim(LandS[steps-2],LandS_2[0])
plt.ylabel(r"$\alpha$ (radians) ")
#plt.title(r"The change in $\alpha$ from $t_c=0$ to 1 Mpc for $m_1=m_2=10^9M_\odot$")
#leg = ax1.legend()
#leg.draggable()
plt.legend(loc="center right",fancybox="True")
plt.grid(True)
plt.hold(True)

#ax2=plt.subplot(112)
ax2=plt.twiny()
plt.grid(True,which='both')
#ax2.xaxis.set_minor_locator(minorLocator)
plt.loglog(LandS_2,freq_array)
plt.axvline(x=f_earth, c="red",linestyle='-.',label= r"$Earth, 10^{-7} Hz$")
plt.axvline(x=f_1kpc,c="k",linestyle='--',label=r"$Earth - 1kpc$")
plt.axvline(x=f_2kpc,c="m",linestyle='--',label=r"$Earth - 2kpc$")
plt.axvline(x=f_10yrs, c="green",linestyle='-.',label=r"$Earth + 10yrs$")
plt.xlabel("frequency (Hz)")
plt.xlim(freq_array[steps-2],freq_array[0])
plt.ylim(5,LandS[0])
plt.grid(True)
plt.hold(True)
#plt.ylim(new[-2],LandS[0])
#plt.ylim(LandS[steps-2],LandS_2[0])
#plt.ylim(alphaS_array[steps-2],LandS[0])
plt.legend(loc="lower right",fancybox="True")
#plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc="upper center",ncol=4, mode="expand", borderaxespad=0.)

ax3=ax1.twinx()
#plt.grid(True)
plt.semilogx(lookback_yrs_array,gamma_array,c='yellow',label=r'$\Gamma,a=0.99$')
#plt.loglog(lookback_yrs_array,gamma_array2,label=r'$\Gamma,a=0.01$')
plt.xlim(lookback_yrs_array[steps-2]-1000,lookback_yrs_array[0])
plt.ylabel(r"$\Gamma=S/L(t)$")
plt.xlabel("lookback time (years)")
plt.grid(True)
plt.legend(loc='upper right',fancybox="True")
plt.savefig("loglogSV10gamma.png")
plt.show()
