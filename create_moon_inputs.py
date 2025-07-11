import os
import pandas as pd
import numpy as np

# Define variables in each contrived history
sim_name = "moon-test-05"
moon_index = 1
final_moon_distance = 2.9 #  Earth radii
final_moon_mass = 1.0 # Moon masses
time_step = 100e3 # years
N = 100

# Mass increment
mass_increment = final_moon_mass / N

# Mass of starting body
start_mass = mass_increment * 1.01

# cumulatively add the masses of the bodies
cumulative_mass = np.arange(start_mass, start_mass+mass_increment * N, mass_increment)

# Array of times
times = np.arange(0, time_step * N, time_step)

# Distance array
distances = np.linspace(2.9, 10.0, N+1)

mass_adjust = 1/ cumulative_mass[-1]


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
    distance_data.append(f"{i+1}\t{distances[i]:.3e}")
    
'''
NB_DATA
'''
NB_data = [sim_name,
           str(moon_index) + "\t (planet #)",
           str(mass_adjust)+ "\t (mass adjust)",
           "0.005 \t (Maximum planetesimal mass)"]

for i in range(len(cumulative_mass)):
    NB_data_entry = (
    f"{times[i]:.2e}\t"
    f"{moon_index}\t"
    f"{cumulative_mass[i]:.3e}\t"
    f"{i+2}\t"
    f"{mass_increment:.3e}\t"
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

        
        

