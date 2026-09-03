#!/bin/bash --login

python /mnt/research/planets/hung/Moon-Accretion/create_moon_inputs.py;
/mnt/research/planets/hung/Moon-Accretion/gen_param_MSE_moon_2;
cp /mnt/research/planets/hung/Moon-Accretion/INPUT_Moon_v5.txt .
