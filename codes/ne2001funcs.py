#!/bin/env python

from subprocess import Popen, PIPE
import numpy as np

bindir='/Users/laura/Research/mygit/Transients/NE2001/bin.NE2001/'

def list_parameters():
    """
    List possible parameters to return from NE2001
    """

    print "DIST      (kpc)                    ModelDistance"
    print "DM        (pc-cm^{-3})             DispersionMeasure"
    print "DMz       (pc-cm^{-3})             DM_Zcomponent"
    print "SM        (kpc-m^{-20/3})          ScatteringMeasure"
    print "SMtau     (kpc-m^{-20/3})          SM_PulseBroadening"
    print "SMtheta   (kpc-m^{-20/3})          SM_GalAngularBroadening"
    print "SMiso     (kpc-m^{-20/3})          SM_IsoplanaticAngle"
    print "EM        (pc-cm^{-6})             EmissionMeasure_from_SM"
    print "TAU       (ms)                     PulseBroadening @1GHz"
    print "SBW       (MHz)                    ScintBW @1GHz"
    print "SCINTIME  (s)                      ScintTime @1GHz @100 km/s"
    print "THETA_G   (mas)                    AngBroadeningGal @1GHz"
    print "THETA_X   (mas)                    AngBroadeningXgal @1GHz"
    print "NU_T      (GHz)                    TransitionFrequency"

    return 

def get_param(l, b, DM, param, ndir=1, verbose=False):
  """
  Runs NE2001 for a set of galactic coordinates and DM (or distance) and returns a single
  parameter. The parameter is choosen using the parameter name given as the second column
  of the NE2001 output. 

  Input: 
  l = galactic longitude
  b = galactic latitude
  DM = dispersion measure (or distance if ndir=0)
  param = string that defines which parameter to return
  ndir = The DM vs D direction (default=1, i.e. DM->D)
  verbose= Print the entire output of NE2001

  Output:
  Parameter requested

  Based on test_ne2001 by J. Cordes
  """

  process=Popen("cd %s; ./NE2001 %s %s %s %d" %(bindir, l, b, DM, ndir), shell=True, stdout=PIPE)
  output = process.communicate()[0]		# output = a long string

  if verbose: print output

  x=output.split()				# x = a list
  val=float(x[x.index(param)-1])

  return val

def get_param_map(ddeg, D=8.3, param='DM'):
    '''Note: originally just for DM; that's way returned value is in "DM" variable'''
    glong=np.array([])
    glat=np.array([])
    DM=np.array([])
    longlist=np.arange(-180, 180+ddeg, ddeg)
    latlist=np.arange(-90, 90+ddeg, ddeg)
    for long in longlist:
        for lat in latlist:
            glong=np.append(glong, long)
            glat=np.append(glat, lat)
            DM=np.append(DM, get_param(long, lat, 50, param, ndir=-1))

    #Reshape
    glong.shape=(len(longlist), len(latlist))
    glat.shape=(len(longlist), len(latlist))
    DM.shape=(len(longlist), len(latlist))

    return glong,glat,DM

