from math import sqrt
import numpy
import math
import random
from numpy import mean
import scipy
import pylab as pl
from utilities import *

#simulation details

Dia = 3e-6 #checked
r = Dia/2
eta_oil = 920e-6 #wikipedia, checked
eta_water = 1.002e-3 #wikipedia, checked
theta = math.radians(130) #contact angle,  Raynaert, 2006
sigma = 7.5e-2 #charge density, checked
alpha_oil = 3.3e-4 #degree of dissosiation, checked, Aveyard 2002
epsilon_zero = 8.8542e-12 #permittivity of free space, checked
epsilon_oil = 2.0 #relative dieelectric constant
gamma_ow = 50e-3 #surface tension oil water, checked
H = 10e-9 #undulation amplitude, Danov 2005, checked
del_phi = math.radians(0)
N = 144 # number of Particles
L = float(18*Dia) # box side length
M = float(((4*math.pi*r**3)/3)*1.055) # mass
T0 = 293 # temperature (Kelvin) 



R = numpy.loadtxt('Equilibrated final position2.txt', delimiter=',')

#x5 = R[:,0]
#y5 = R[:,1]
#
#
#pl.scatter(x5, y5,marker='.', color='k', s=20, alpha=1)
#pl.xlim( (-L/2, L/2) )
#pl.ylim( (-L/2, L/2) )
#pl.show()

plot_adsorbed_circles(R[:,0]+L/2, R[:,1]+L/2, r, L, R)