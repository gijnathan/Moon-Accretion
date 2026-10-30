#!/bin/bash
set -euo pipefail

N=100 # number of simplexes
X_TM_min=0.1
X_TM_max=0.9
dX_TM=0.05

ROOT_DIR="$(pwd)"
GEN_EXEC="$ROOT_DIR/gen_param_MSE_moon_2"

rm -rf MSE_fit/*
mkdir -p MSE_fit

# Loop over X_TM values
for i in $(awk -v min="$X_TM_min" -v max="$X_TM_max" -v step="$dX_TM" \
    'BEGIN{
        n = int((max-min)/step + 0.5)  # number of steps, rounded to nearest int
        for(k=0; k<=n; k++) printf "%.2f\n", min + k*step
    }'); do
    echo "X_TM $i"
    mkdir -p "MSE_fit/XTM_$i"

    # Create simplexes, copy data and generate PARAM_MSE for specific X_TM value
    for ((n=1; n<=N; n++)); do
        SIMDIR="$ROOT_DIR/MSE_fit/XTM_$i/Simplex_$n"
        mkdir -p "$SIMDIR"
        cp -r data/* "$SIMDIR"
        (
            cd "$SIMDIR"
            "$GEN_EXEC" "$i"
        )
        (( n % 10 == 0 )) && echo "$n"
    done
    cp submit.sb "$ROOT_DIR/MSE_fit/XTM_$i"
    # Generate parallel_commands.txt
    find "$ROOT_DIR/MSE_fit/XTM_$i" -type d -name "Simplex_*" | sort | \
        awk '{print "cd " $0 " && N_body_differentiate_Moon_v5 > log.txt"}' \
        > $ROOT_DIR/MSE_fit/XTM_$i/parallel_commands.txt    
done