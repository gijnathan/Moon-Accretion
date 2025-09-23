Using NB-Differentiate to model lunar accretion matching observed lunar mantle from Warren 2005 and Delano 1986
21 September 2025: Finding out why ALL versions of code are segfaulting during the minimzation procedure. Able to narrow the issue down to an issue with func(x) (the function that calculates the chi sq inside of the code). One potential reason might be the new_starting_molar_proportions function, since it was the only thing we really edited. Not sure what is causing this segfault though as calculations are still being done.  

3 July 2025: Make sure future INPUT_Moon_v2 files have dist_3 set in Earth radii. Before, dist_3 was seet to 4.5 (AU), which caused bodies beyond 4.5 to gain water too fast and oxidize the mantle. 

2 July 2025: Created new python script called create_moon_inputs.py which generates a contrived accretion history for the moon. From 100 homogenous moonlets, create moon that is at 2.9 Earth radii, with mass of 1.0 Moon masses. Script generates three files: distance_data, NB_data, and final_planets. In final_planets, there is only one body, which is the moon. Three final columns for NB_data were kept constant since tthey currently are unknown variables. 

moon-test-05 is the first version to use this contrived accretion history with 100 impacts. Checked and works fine with code. 

27 June 2025: Modified N\_body\_differentiate\_Moon to version 2. Now INPUT is INPUT\_Moon\_v2.txt. For all future INPUT files, use v2.

