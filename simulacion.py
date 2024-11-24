#Importy the necesary libreries

import numpy as np
from matplotlib import pyplot as plt
import scipy
import math



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
counter=0 #counter for arrays
Nuclei_in_CR=(4/3)*np.pi*(No_of_nuclei_D/2)**3 # number of nuclei in the core radius
deltat=0.5  # time step.
n=0

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
TimeSteps=int((Tmax - Ti) / deltat)
LiveNeutronsOverTime=np.zeros(TimeSteps)
#Start of simulation

for t in np.arange(Ti, Tmax, deltat): # Loop through the time steps
    while N_Live_new < 1e5: # Stop when the number of live
        # neutron are bigger than 1x10^{5}
        Neutron_x_d2=x/IN_distence # X neutron position in terms of inter nucleaer distance (d).
        Neutron_y_d2=y/IN_distence
        Neutron_z_d2=z/IN_distence

        Max_Neutron_x_d = np.ceil(Neutron_x_d2) #Max integer number for x values (nuclei poition)
        Min_Neutron_x_d= np.floor(Neutron_x_d2) #Min integer number for x values (nuclei poition)

        Max_Neutron_y_d = np.ceil(Neutron_x_d2)
        Min_Neutron_y_d= np.floor(Neutron_x_d2)

        Max_Neutron_z_d = np.ceil(Neutron_x_d2)
        Min_Neutron_z_d= np.floor(Neutron_x_d2)

        DistanceA=np.sqrt((x-(Max_Neutron_x_d*IN_distence))**2+(y-(Max_Neutron_y_d*IN_distence))**2+(z-(Max_Neutron_z_d*IN_distence))**2)

        DistanceB=np.sqrt((x-(Max_Neutron_x_d*IN_distence))**2+(y-(Max_Neutron_y_d*IN_distence))**2+(z-(Min_Neutron_z_d*IN_distence))**2)
        
        DistanceC=np.sqrt((x-(Max_Neutron_x_d*IN_distence))**2+(y-(Min_Neutron_y_d*IN_distence))**2+(z-(Min_Neutron_z_d*IN_distence))**2)
        
        DistanceD=np.sqrt((x-(Min_Neutron_x_d*IN_distence))**2+(y-(Min_Neutron_y_d*IN_distence))**2+(z-(Min_Neutron_z_d*IN_distence))**2)
        
        DistanceE=np.sqrt((x-(Min_Neutron_x_d*IN_distence))**2+(y-(Max_Neutron_y_d*IN_distence))**2+(z-(Min_Neutron_z_d*IN_distence))**2)
        
        DistanceF=np.sqrt((x-(Min_Neutron_x_d*IN_distence))**2+(y-(Max_Neutron_y_d*IN_distence))**2+(z-(Max_Neutron_z_d*IN_distence))**2)
        
        DistanceG=np.sqrt((x-(Min_Neutron_x_d*IN_distence))**2+(y-(Min_Neutron_y_d*IN_distence))**2+(z-(Max_Neutron_z_d*IN_distence))**2)
        
        DistanceH=np.sqrt((x-(Max_Neutron_x_d*IN_distence))**2+(y-(Min_Neutron_y_d*IN_distence))**2+(z-(Max_Neutron_z_d*IN_distence))**2)

        
        if (DistanceA<Nuclear_R).any():
            n=np.random.rand()
        elif (DistanceC<Nuclear_R).any():
            n=np.random.rand()
        elif (DistanceD<Nuclear_R).any():
            n=np.random.rand()
        elif (DistanceE<Nuclear_R).any():
            n=np.random.rand()
        elif (DistanceF<Nuclear_R).any():
            n=np.random.rand()
        elif (DistanceG<Nuclear_R).any():
            n=np.random.rand()
        elif (DistanceH<Nuclear_R).any():
            n=np.random.rand()
        
        if n<Fission_CS/CS_total: #if true the collision is declared a fission event (appox 0.217) 
            N_Live_new=N_initial+N_PerFission
            np.append(x,3354e-15)
            np.append(y,2502e-15)
            np.append(z,1564e-15)
            
            np.append(x,N_speed*np.random.rand())
            np.append(y,N_speed*np.random.rand())
            np.append(z,N_speed*np.random.rand())
        
        
            LiveNeutronsOverTime[counter]=N_Live_new
            N_initial=N_Live_new
            counter+=1
        
        else:
            continue #Agregar el proceso de scattering
    break  # Remove this line after implementing actual operations
