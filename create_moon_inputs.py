import os
import pandas as pd
import numpy as np

# Define variables in each contrived history
sim_name = "moon-test-05"
moon_index = 1
final_moon_distance = 2.9

# Indicies for each moonlet, used in distance data and NB data
moonlet_indicies = np.arange(1,101,1)

# Time in NB_data column
time = np.linspace(0, 100, 100)

# Cumulative mass of the growing moon following each impact
mass = np.linspace(0.1, 1.00, 100)
final_moon_mass = mass[-1] # moon_mass is defined as the last item in here

# Impactor mass, simplified to even homogenous impact events
impactor_mass = (1.0 - 0.1)/100

distances = np.linspace(2.9, 10.0, 100)




# ---------------------Create lists to store each file's data-----------------------------
'''
FINAL PLANETS
'''

final_planets = [sim_name, 
                 1, 
                 " 1 " + str(moon_index) + " " + str(final_moon_mass) + " " + str(final_moon_distance) + " 0.0376  2.1195"]


'''
DISTANCE DATA
'''
distance_data = [sim_name,
                 0.0000]
# Get info for each line in distance data
# each row contains index of moonlet and distance in Earth radii 
for i in range(len(distances)):
    distance_data.append(f"{moonlet_indicies[i]}\t{distances[i]:.3e}")
    
'''
NB_DATA
'''
NB_data = [sim_name,
           str(moon_index) + "\t (planet #)",
           "1.0 \t (mass adjust)",
           "0.005 \t (Maximum planetesimal mass)"]

for i in range(len(distances)-1):
    NB_data_entry = (
    f"{time[i]:.2e}\t"
    f"{moon_index}\t"
    f"{mass[i+1]:.3e}\t"
    f"{moonlet_indicies[i+1]}\t"
    f"{impactor_mass:.3e}\t"
    "2.6774\t1.006\t79.59"
)
    NB_data.append(NB_data_entry)
 # Get info for each line in NB_data
 

# --------------WRITE FILES-------------------- 
with open("final_planets.txt", "w") as f:
    for item in final_planets:
        f.write(f"{item}\n")
        
with open("distance_data.txt", "w") as f:
    for item in distance_data:
        f.write(f"{item}\n")
        
        
with open("NB_data.txt", "w") as f:
    for item in NB_data:
        f.write(f"{item}\n")

        
        

