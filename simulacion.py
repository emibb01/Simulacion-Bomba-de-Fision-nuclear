#Importy the necesary libreries

import numpy as np
from matplotlib import pyplot as plt
import scipy


#Variables for the simulation

Fission_CS=1.235e-28   #Fission cross section units meters
Scattering_CS=4.566e-28    #Scattering cross section units meters
CS_total= Fission_CS+Scattering_CS #Total cross section area (each nuclei presents to neutron)
N_PerFission=5 #Neutrons per fission
IN_distence=160e-15    #Inter nuclear distance (d)(fm)
N_speed= 19.56e-15     #Speed of neutrons
Core_R=42000e-15       #Core Radius (fm)
Ti=1               #initial time for the simulation
Tmax= 25000    #Max time the simulation can run(zs)
N_initial=100      #Initial amount of neutrons
No_of_nuclei_D= (2*Core_R/IN_distence)+1       #number of nuclei in the core diameter (max 600)
NCounter=1; #Neutron counter for arrays
Nuclear_R=np.sqrt(CS_total/np.pi)    #Finds the nuclear radius in the core in which events can happen
N_Live_new=0 #Number of live neutrons in the core
i=1 #index counter for arrays
counter=1; #counter for arrays
Nuclei_in_CR=(4/3)*np.pi*(No_of_nuclei_D/2)**3 # number of nuclei in the core radius
deltat=0.5  # time step.


    # Generating random radial angles between 0 and 2π
phi = 2 * np.pi * np.random.rand(N_initial, 1)

    # Generating theta between 0 and π
theta = np.pi * np.random.rand(N_initial, 1)

    # Generating random radii between 0 and Core_R
rad = Core_R * np.random.rand(N_initial, 1)

    # Converting spherical coordinates to Cartesian coordinates
x = rad * np.sin(theta) * np.cos(phi)
y = rad * np.sin(theta) * np.sin(phi)
z = rad * np.cos(theta)
    
    #contains the Cartesian coordinate   
cordenates=[x, y, z] #contains the Cartesian coordinates
    
    # Creating arrays filled with a specific value
Vx = np.full(N_initial, N_speed)  # Velocity of x
Vy = np.full(N_initial, N_speed)  # Velocity of y
Vz = np.full(N_initial, N_speed)  # Velocity of z

#Start of simulation

for t in np.arange(Ti, Tmax - deltat, deltat): # Loop through the time steps
    while N_Live_new < 1e5: # Stop when the number of live
        # neutron are bigger than 1x10^{5}
        
        
        
        break  # Remove this line after implementing actual operations
