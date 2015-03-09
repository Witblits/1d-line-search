#!/usr/bin/python
#This script is based on the Particle Swarming algorithm which is used to approximate the discrete global optimum.

from __future__ import division				#ensures division always returns the correct answer (eg. 1/2 = 0.5 instead of 1/2 = 0)
import argparse
import pylab as pl
import os
##########################################################################################
total_runs = 0
pop = 0                  		#population size
n_dim = 0                 		#number of dimensions
##########################################################################################