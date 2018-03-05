#K Grover
# Plots the files created in defltaf_PN_contour.py

#!/usr/bin/env python

from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, ticker, rc


rc('text', usetex=True)
plt.rcParams['font.size'] = 20


def plotDeltaf(log_m1_solar_array, log_fE_array, diff_array, fignum):
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
	plt.xlabel(r"$log(\frac{m_1}{M_\odot})$")
	plt.ylabel(r"$log(\frac{f_E}{Hz})$")
	plt.savefig("deltaf_PN.png")
	
if __name__=="__main__":
	file1=open("deltaf_array.txt", 'r')
	file2=open("log_m1_solar_array.txt", 'r')
	file3=open("log_fE_array.txt", 'r')
	
	deltaf_array=np.loadtxt(file1, dtype=float, delimiter=',')
	log_m1_solar_array=np.loadtxt(file2, dtype=float, delimiter=',') 
	log_fE_array=np.loadtxt(file3, dtype=float, delimiter=',')
	
#uncomment to find largest difference
	#maximum=np.amax(deltaf_array)
	#print "\n The maximum number in the deltaf_array is %e" %maximum
	#a=np.argmax(deltaf_array)
	#(i,j)=np.unravel_index(a,deltaf_array.shape)
	#print "\n This occurs at mass = %e" %(10**(log_m1_solar_array[i]))
	#print "and freq = %e" %(10**(log_fE_array[j]))
	#print "\n It corresponds to an fmin of %e s" %(1./maximum)
	#print "This is %e years" %((1./maximum)/31556926)
	
	plotDeltaf(log_m1_solar_array, log_fE_array, deltaf_array, 1)

	