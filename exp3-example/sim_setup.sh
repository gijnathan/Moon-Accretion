#!/bin/bash
set -euo pipefail

N=100 # number of simplexes
firstbody_metal_min=0.1
firstbody_metal_max=0.9
dX_TC_1=0.05

ROOT_DIR="$(pwd)"
GEN_EXEC="$ROOT_DIR/gen_param_MSE_moon_2"

rm -rf MSE_fit/*
mkdir -p MSE_fit

# Loop over firstbody_metal values
for i in $(awk -v min="$firstbody_metal_min" -v max="$firstbody_metal_max" -v step="$dX_TC_1" \
    'BEGIN{
        n = int((max-min)/step + 0.5)  # number of steps, rounded to nearest int
        for(k=0; k<=n; k++) printf "%.2f\n", min + k*step
    }'); do
    echo "firstbody_metal $i"
    mkdir -p "MSE_fit/firstbody_metal_$i"

    # Create simplexes, copy data and generate PARAM_MSE for specific firstbody metal value
    for ((n=1; n<=N; n++)); do
        SIMDIR="$ROOT_DIR/MSE_fit/firstbody_metal_$i/Simplex_$n"
        mkdir -p "$SIMDIR"
        cp -r data/* "$SIMDIR"
        sed -i "\$s/.*/${i}\tfirstbody_metal_frac/" "$SIMDIR/INPUT_Moon_v5_firstbody_metal.txt"
        (
            cd "$SIMDIR"
            "$GEN_EXEC" 0.2
        )
        (( n % 10 == 0 )) && echo "$n"
    done
    cp submit.sb "$ROOT_DIR/MSE_fit/firstbody_metal_$i"
    # Generate parallel_commands.txt
    find "$ROOT_DIR/MSE_fit/firstbody_metal_$i" -type d -name "Simplex_*" | sort | \
        awk '{print "cd " $0 " && N_body_differentiate_Moon_v5_firstbody_metal > log.txt"}' \
        > $ROOT_DIR/MSE_fit/firstbody_metal_$i/parallel_commands.txt    
done