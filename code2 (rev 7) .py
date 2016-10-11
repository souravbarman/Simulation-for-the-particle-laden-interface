# 2D Molecular dynamics simulations  of 

from math import sqrt
import numpy
import math
import random
from numpy import mean
import scipy
import pylab as pl

#simulation details

Dia = 3e-6
r = Dia/2
eta_oil = 920e-6
eta_water = 1.002e-3
theta = math.radians(130)
sigma = 7.5e-2 #charge density
alpha_oil = 3.3e-4 #degree of dissosiation
epsilon_zero = 8.8542e-12
epsilon_oil = 2.0
gamma_ow = 50e-3 #surface tension oil water
H = 10e-9 #undulation amplitude
del_phi = math.radians(0)
N = 144 # number of Particles
L = float(18*Dia) # box side length
M = float(4*math.pi*r**3/3*1.055) # mass
T0 = 293 # temperature (Kelvin) 

eta_eff = float(eta_oil*(1-math.cos(theta)) + eta_water*(1+math.cos(theta)))

b = float(6*math.pi*r*eta_eff)

zeta = float(r*(3+math.cos(theta))/2) 

q_oil = float(2*math.pi*(r**2)*sigma*(1-math.cos(theta))*alpha_oil)

constant1 = float(3*(q_oil**2)*(zeta**2)/(8*math.pi*epsilon_zero*epsilon_oil))

constant2 = -float(48*math.pi*gamma_ow*(H**2)*math.cos(del_phi)*(r**4))

constant3 = float(q_oil**2/(4*math.pi*epsilon_zero*epsilon_oil))

surface_coverage = float((math.pi*r**2*N)/L**2)

L_offset = float(L*(1-1/sqrt(N))/2)

#def InitPositionSquare(N,L):
#
#  position = numpy.zeros((N,2)) + 0.0
#  added = 0
#  for x in range(0, int(sqrt(N))):
#    for y in range(0, int(sqrt(N))):
#        if(added < N):
#          position[added, 0] = x * float(L)/(float(sqrt(N))-0.999) - float(L)/2
#          position[added, 1] = y * float(L)/(float(sqrt(N))-1) - float(L)/2
#           
#          added += 1
#  
#         
#  x = position[:,0]*1e6
#  y = position[:,1]*1e6
#  
#  pl.subplot(3,1,1)
#  pl.scatter(x,y)
#  pl.show()
#  
#  return position

def InitPositionSquare(N,L):
  position = numpy.zeros((N,2)) + 0.0
  Nsquare = 1
  while(N > (Nsquare*Nsquare)):
    Nsquare += 1
  if(Nsquare**2 != N):
    print "Warning: Your particle number",N, \
          "is not a perfect square; this may result " \
          "in a lousy initialization"
  rs = float(L)/Nsquare
  roffset = float(L)/2 - rs/2
  added = 0
  for x in range(0, Nsquare):
    for y in range(0, Nsquare):
        if(added < N):
          position[added, 0] = rs*x - roffset 
          position[added, 1] = rs*y - roffset 
           
          added += 1
          
  x = position[:,0]*1e6
  y = position[:,1]*1e6
  
  
  pl.subplot(3,2,1)
  pl.scatter(x,y)
  pl.show()
  
  return position
#
#def InitPositionSquare(N,L):
#    
#  position = numpy.zeros((N,2)) + 0.0
#  position = -L*numpy.random.rand(N,2)+L/2
#  
#  x = position[:,0]*1e6
#  y = position[:,1]*1e6
#  
#  pl.subplot(2,1,1)
#  pl.scatter(x,y)
#  pl.show()
#  
#  return position

#def InitVelocity(N,T0,mass=1.):
#  initNDIM = 2
#  velocity = numpy.zeros((N,2)) + 0.0
#  random.seed(1)
#  netP = numpy.zeros((2,)) + 0.
#  netE = 0.0
#  for n in range(0, N):
#    for x in range(0, initNDIM):
#      newP = random.random()-0.5
#      netP[x] += newP
#      netE += newP*newP
#      velocity[n, x] = newP
#  netP *= 1.0/N
#  vscale = math.sqrt(3*N*T0/(mass*netE))
#  for n in range(0, N):
#    for x in range(0, initNDIM):
#      velocity[n, x] = (velocity[n, x] - netP[x]) * vscale
#  return velocity
  


#def PutInBox(Ri):
#    
#        while Ri[0]>L/2:
#           Ri[0]=Ri[0]-L 
#           if Ri[0]<0:
#            break
#
#        while Ri[0]<-L/2:
#           Ri[0]=Ri[0]+L 
#           if Ri[0]>0:
#            break
#        
#            
#        return Ri

def PutInBox(Ri):
    
        while Ri[0]>L/2:
           Ri[0]=Ri[0]-L 
           if Ri[0]<0:
            break

        while Ri[0]<-L/2:
           Ri[0]=Ri[0]+L 
           if Ri[0]>0:
            break
        
        while Ri[1]>(L/2-Dia):
           Ri[1]=Ri[1]-L+2*Dia 
           if Ri[1]<0:
            break

        while Ri[1]<-(L/2-Dia):
           Ri[1]=Ri[1]+L-2*Dia 
           if Ri[1]>0:
            break
            
        return Ri
        
        

def new_particle_velocity(i, R, V,delta_t):     
    
    Vx = 0.0
    Vy = 0.0
    
    average_V = [0.0,0.0]
    
    FE = 0.0
    FC = 0.0  
    #print V
    for j in range (0,len(R)):
      if j!=i:       
       average_V = (V.mean(axis=0) * N - V[i])/(N-1)
#       b = float(5.2853e-8)
#       constant1 = float(2.78433e-32)
#       constant2 = float(3.81704e-39)
       
       D = R[i]-R[j]
       D = PutInBox (D) ##  interpartcile displacement 
       d = sqrt((D[0]**2)+(D[1]**2))  #interparticle distance
       
#       print float((zeta/d)**2)
       
       if float((zeta/d)**2) < 1e-6:
           FE = FE + constant1 * (1/float(d)**4) + float((-1)*12*math.pi*gamma_ow*H**2*math.cos(del_phi)*r**4/(d**4))
       else:
           FE = FE + constant3 * (1/float(d)**2 - float(d)/float((4*zeta**2 + d**2)**1.5)) + float((-1)*12*math.pi*gamma_ow*H**2*math.cos(del_phi)*r**4/(d**4))
           
       FC =  FC + constant2 * (1/float(d)**5)

    Vx =  average_V[0] + 1/float(b) * (FE *D[0]/float(d)+ FC*D[0]/float(d))
    Vy = average_V[1]+ 1/float(b) * (FE *D[1]/float(d)+ FC*D[1]/float(d))
      
    #print x_new, y_new
    return Vx, Vy

#def ComputeEnergy(R,V):
#
#    KE=0.0
#    PE=0.0
#    TE=0.0
#
#    for i in range (0,len(R)):
#       for j in range (i+1, len(R)):
#         d= Distance(R[i],R[j])
#         PE = PE + 4*(1/(d**12)-1/(d**6))
#                      
#
#    for i in range (0,len(V)):
#        for j in range (0,len(V[i])):
#           KE=KE+0.5*M*V[i][j]**2
#           
#    TE = PE+KE        
#    
#    return TE

def ComputeEnergy(R,V):

    KE=0.0
    PE=0.0
    TE=0.0

    for i in range (0,len(R)):
       for j in range (i+1, len(R)):
         D = R[i]-R[j] # Distance
         D = PutInBox (D) # Periodic BC
         d = sqrt((D[0]**2)+(D[1]**2))
         if float((zeta/d)**2) < 1e-6:
             PE = PE + constant1 * (-1/(3*float(d)**3)) + float((-1)*12*math.pi*gamma_ow*H**2*math.cos(del_phi)*r**4/(d**4))
         else:
             PE = PE + constant3 * (1/float(d) - float(d)/float((4*zeta**2 + d**2)**0.5)) + float((-1)*12*math.pi*gamma_ow*H**2*math.cos(del_phi)*r**4/(d**4))                  

    for i in range (0,len(V)):
        for j in range (0,len(V[i])):
            KE = KE + 0.5*M*V[i][j]**2    
            
           
    TE = PE + KE        
    
    return TE

#MAIN PROGRAM 

R = numpy.zeros((N,2)) + 0.0


R = InitPositionSquare(N, L)    # non-zero position initilization 
R_list = R.tolist() 
R_list_for_bulk = R.tolist() 

Upper_R = numpy.asarray([v for v in R_list if v[1] == max(R[:,1])])
Lower_R = numpy.asarray([v for v in R_list if v[1] == min(R[:,1])])
 
nR = numpy.zeros((N,2)) + 0.0
nV = numpy.zeros((N,2)) + 0.0

Upper_V = numpy.zeros((sqrt(N),2)) + 0.0
Lower_V = numpy.zeros((sqrt(N),2)) + 0.0

Upper_nV = numpy.zeros((sqrt(N),2)) + 0.0
Lower_nV = numpy.zeros((sqrt(N),2)) + 0.0


Upper_nR = numpy.zeros((sqrt(N),2)) + 0.0
Upper_nR[:,1] = Upper_nR[:,1]+max(R[:,1])

Lower_nR = numpy.zeros((sqrt(N),2)) + 0.0
Lower_nR[:,1] = Lower_nR[:,1]+min(R[:,1])

Bulk_nR = numpy.zeros(((N-2*sqrt(N)),2)) + 0.0
Bulk_nV = numpy.zeros(((N-2*sqrt(N)),2)) + 0.0
#V = InitVelocity(N, T0, M)

############################        
delta_t = float(1e-6)
steps = 1000
shear_rate = 0.98

############################


##### Wall equilibration ########

for i in range(0, len(Upper_R)):

    R_list_for_bulk.remove(Upper_R[i].tolist())
    R_list_for_bulk.remove(Lower_R[i].tolist())

Bulk_R = numpy.asarray(R_list_for_bulk)
    

for t_equilibration_walls in range(0,steps): 
    for i in range(0,len(Upper_R)):     
        Upper_nV[i] = new_particle_velocity(i, Upper_R, Upper_V, delta_t)
        Upper_nR[i][0] = Upper_R[i][0] + nV[i][0] * delta_t   
        Lower_nV[i] = new_particle_velocity(i, Lower_R, Lower_V, delta_t)
        Lower_nR[i][0] = Lower_R[i][0] + nV[i][0] * delta_t  
        Upper_V[i] = Upper_nV[i].copy()
        Lower_V[i] = Upper_nV[i].copy()
        PutInBox(Upper_nR[i])
        PutInBox(Lower_nR[i])  

    Upper_R = Upper_nR.copy() #instead of R=nR, nR.copy() made it lot faster
    Upper_V = Upper_nV.copy()
    Lower_R = Lower_nR.copy() #instead of R=nR, nR.copy() made it lot faster
    Lower_V = Lower_nV.copy()


Upper_R_for_nonequil = Upper_R.copy()
Upper_V_for_nonequil = Upper_V.copy()
Lower_R_for_nonequil = Lower_R.copy()
Lower_V_for_nonequil = Lower_V.copy()



##### Bulk equilibration #######

V = numpy.zeros((N,2)) + 0.0
 
R = numpy.asarray(Upper_R.tolist() + Bulk_R.tolist() + Lower_R.tolist())
Init_R = R
Bulk_R = (R_list_for_bulk)

R_List1 = R.tolist()

#Upper_V = numpy.zeros((sqrt(N),2)) + 0.0
#Lower_V = numpy.zeros((sqrt(N),2)) + 0.0
TE = numpy.zeros((steps,1)) + 0.0
#Energy = []
for t_equilibration_bulk in range(0,steps): 
    
    #print t_equilibration_bulk, ComputeEnergy(R,V)
        
    TE[t_equilibration_bulk] = ComputeEnergy(R,V).copy()  
    
    pl.subplot(3,2,3)
    pl.scatter(t_equilibration_bulk, TE[t_equilibration_bulk]*1e13)
    pl.show()

    for I in range(0,len(Bulk_R)): 
#        p = R_List1.index(Bulk_R[I].tolist())
        p = numpy.shape(Upper_R)[0] + I 
#        print p        
        Bulk_nV[I] = new_particle_velocity(p, R, V, delta_t)
        V[p] = Bulk_nV[I].copy()
        Bulk_nR[I] = Bulk_R[I] + Bulk_nV[I] * delta_t          
        PutInBox(Bulk_nR[I])
        

    Bulk_R = Bulk_nR.copy() #instead of R=nR, nR.copy() made it lot faster
    Bulk_V = Bulk_nV.copy()
    Bulk_R_non = Bulk_nR.copy() 
    Bulk_V_non = Bulk_nV.copy()
    R = numpy.asarray(Upper_R.tolist() + Bulk_R.tolist() + Lower_R.tolist())
    V = numpy.asarray(Upper_V.tolist() + Bulk_V.tolist() + Lower_V.tolist())
    
    
R = numpy.asarray(Upper_R.tolist() + Bulk_R.tolist() + Lower_R.tolist())
V = numpy.asarray(Upper_V.tolist() + Bulk_V.tolist() + Lower_V.tolist())

Bulk_R_for_nonequil = Bulk_R_non.copy() 
Bulk_V_for_nonequil = Bulk_V_non.copy()

R_equil = R.copy()
V_equil = V.copy()
Energy = TE.copy()

x5 = R[:,0]*1e6
y5 = R[:,1]*1e6

pl.subplot(3,2,4)
pl.scatter(x5, y5)
pl.show()



numpy.savetxt('Equilibrated Upper_R.txt', Upper_R_for_nonequil, delimiter=',')
numpy.savetxt('Equilibrated Lower_R.txt', Lower_R_for_nonequil, delimiter=',')
numpy.savetxt('Equilibrated Upper_V.txt', Upper_V_for_nonequil, delimiter=',')
numpy.savetxt('Equilibrated Lower_V.txt', Lower_V_for_nonequil, delimiter=',')

numpy.savetxt('Equilibrated Bulk_R.txt', Bulk_R_for_nonequil, delimiter=',')
numpy.savetxt('Equilibrated Bulk_V.txt', Bulk_V_for_nonequil, delimiter=',')

numpy.savetxt('Equilibrated Energy state.txt', Energy)

##### Applying Shear ####

Upper_R = Upper_R_for_nonequil
Lower_R = Lower_R_for_nonequil
Upper_V = Upper_V_for_nonequil
Lower_V = Lower_V_for_nonequil

Bulk_R = Bulk_R_for_nonequil
Bulk_V = Bulk_V_for_nonequil

R = numpy.asarray(Upper_R.tolist() + Bulk_R.tolist() + Lower_R.tolist())
V = numpy.asarray(Upper_V.tolist() + Bulk_V.tolist() + Lower_V.tolist())

numpy.savetxt('Equilibrated final position.txt', R, delimiter=',')
numpy.savetxt('Equilibrated final velocity.txt', V, delimiter=',')

shear_displacement = numpy.zeros((sqrt(N),2)) + 0.0
TE_non = numpy.zeros((steps,1)) + 0.0
 
#for t_non_equilibration_shear in range(0,steps):
#    
#    shear_displacement[:,0] = shear_rate*L*delta_t*0.5 # assuming couette cell like flow and shear rate = (relative velocity)/L
#    Upper_R = Upper_R + shear_displacement
#    Lower_R = Lower_R - shear_displacement
#    Upper_V[:,0] = Upper_V[:,0] + shear_rate*L*0.5
#    Lower_V[:,0] = Lower_V[:,0] - shear_rate*L*0.5
#    TE_non[t_non_equilibration_shear] = ComputeEnergy(R,V).copy()    
#  
#    
#    pl.subplot(3,2,5)
#    pl.scatter(t_non_equilibration_shear, TE_non[t_non_equilibration_shear]*1e13)
#    pl.show()
#    
#           
#    for i in range(0,len(Upper_R)):     
#        
#        PutInBox(Upper_R[i])
#        PutInBox(Lower_R[i])  
#
#    Upper_R = Upper_R.copy() #instead of R=nR, nR.copy() made it lot faster
#    Lower_R = Lower_R.copy() #instead of R=nR, nR.copy() made it lot faster
#    R = numpy.asarray(Upper_R.tolist() + Bulk_R.tolist() + Lower_R.tolist())
#    V = numpy.asarray(Upper_V.tolist() + Bulk_V.tolist() + Lower_V.tolist())
#    
#    
#    for I in range(0,len(Bulk_R)): 
##        p = R_List1.index(Bulk_R[I].tolist())
#        p = numpy.shape(Upper_R)[0] + I 
##        print p        
#        Bulk_nV[I] = new_particle_velocity(p, R, V, delta_t)
#        V[p] = Bulk_nV[I].copy()
#        Bulk_nR[I] = Bulk_R[I] + Bulk_nV[I] * delta_t          
#        PutInBox(Bulk_nR[I])
#        
#
#    Bulk_R = Bulk_nR.copy() #instead of R=nR, nR.copy() made it lot faster
#    Bulk_V = Bulk_nV.copy()
##    Upper_R = PutInBox(Upper_R)    
##    Lower_R = PutInBox(Lower_R)
#    R = numpy.asarray(Upper_R.tolist() + Bulk_R.tolist() + Lower_R.tolist())
#    V = numpy.asarray(Upper_V.tolist() + Bulk_V.tolist() + Lower_V.tolist())
#
#R = numpy.asarray(Upper_R.tolist() + Bulk_R.tolist() + Lower_R.tolist())
#V = numpy.asarray(Upper_V.tolist() + Bulk_V.tolist() + Lower_V.tolist()) 
  

x1 = R[:,0]*1e6
y1 = R[:,1]*1e6


#x2 = Upper_R[:,0]*1e6
#y2 = Upper_R[:,1]*1e6
#x3 = Lower_R[:,0]*1e6
#y3 = Lower_R[:,1]*1e6
pl.subplot(3,2,6)
pl.scatter(x1, y1)
#pl.scatter(x2, y2)
#pl.scatter(x3, y3)
pl.show()

