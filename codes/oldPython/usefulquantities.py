#
#  usefulquantities.py
#  
#
#  Created by alberto vecchio on 14/09/2011.
#  Copyright (c) 2011 university of birmingham. All rights reserved.
#

# Import modules
from math import *
import scipy as sp
import scipy.integrate as ig
import numpy as np
from scipy.linalg import inv
import matplotlib.pyplot as plt     # If you wish to plot a graph of the integrands

# Converts solar mass to seconds
def MassSeconds(m):
    m_secs=m*4.925*pow(10.,-6.)
    return m_secs

# The main program

if __name__=="__main__":

	# masses in solar masses
	m1 = 1.0e9
	m2 = 1.0e9

	# other mass parameters
	mtot = m1 + m2
	mred = m1*m2/mtot
	mchirp = pow(mtot,(2.0/5.0))*pow(mred,(3.0/5.0))
	eta = mred/mtot
		
	# GW frequency
	f = 1.0e-7
		
	# observation time
	Tobs = 3.1e7*10.0
	L1kpc = 3.08e18*1000.0/3.e10
	
	# orbital velocity (v/c) 
	vorb = pow((pi*MassSeconds(mtot)*f),(1.0/3.0))
	print "orbital velocity (v/c)",vorb
	
	# PN contributions
	dfdt_nw = (96.0/5.0)*pow(pi,(8.0/3.0))*pow(MassSeconds(mchirp),(5.0/3.0))*pow(f,(11.0/3.0))
	print "fdot and deltaF at NW",dfdt_nw, dfdt_nw*pow(Tobs,2), dfdt_nw*L1kpc*Tobs	
	
	dfdt_1pn = dfdt_nw*((743.0/336.0) + (11.0*eta/4.0))*pow((pi*MassSeconds(mtot)*f),(2.0/3.0))
	print "fdot and deltaF at 1PN",dfdt_1pn,dfdt_1pn*pow(Tobs,2), dfdt_1pn*L1kpc*Tobs
	
	dfdt_15pn = dfdt_nw*4*pi*(pi*MassSeconds(mtot)*f)
	print "fdot and deltaF at 1.5 (no spin) PN",dfdt_15pn, dfdt_15pn*pow(Tobs,2), dfdt_15pn*L1kpc*Tobs
	print "fdot and deltaF at 1.5 SO (in unit of beta) PN", -dfdt_15pn/(4.0*pi), -dfdt_15pn*pow(Tobs,2)/(4.0*pi), -dfdt_15pn*L1kpc*Tobs/(4.0*pi)
	
	dfdt_2pn =  dfdt_nw*((34103.0/18144.0) + (13661.0/2016.0)*eta + (59.0/18.0)*pow(eta,2))*pow((pi*MassSeconds(mtot)*f),(4.0/3.0))
	dfdt_2pn_ss =  dfdt_nw*pow((pi*MassSeconds(mtot)*f),(4.0/3.0))
	print "fdot and deltaF at 2PN SS (in unit of sigma) PN", dfdt_2pn_ss, dfdt_2pn_ss*pow(Tobs,2), dfdt_2pn_ss*L1kpc*Tobs
	