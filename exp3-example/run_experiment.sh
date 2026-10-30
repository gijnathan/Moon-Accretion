#!/bin/bash

ROOT_DIR="$(pwd)"

for d in "$ROOT_DIR"/MSE_fit/firstbody_metal*; do
    # find the number of simplexes in each simulation to pass ntasks automatically
    n_simplex=$(find "$d" -maxdepth 1 -type d -name "Simplex_*" | wc -l)
    if [ "$n_simplex" -eq 0 ]; then
        echo "Warning: no simplex directories found in $d, skipping"
        continue
    fi
    (cd "$d" && sbatch --ntasks-per-node="$n_simplex" submit.sb)
done