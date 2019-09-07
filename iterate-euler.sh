#!/usr/bin/env bash

export SHELL=$(type -p bash)

function iterate {
    ncells=$1
    mid=$(( (ncells / 2) + 1))
    # Run algorithm - store time and midpoint in var
    var=$( TIMEFORMAT="%U"; { time ./1d_euler -n "${ncells}" -t 0.15 -d 1 -D 1 -v -2 -V 2 -p 0.4 -P 0.4 -f | cut -f 3 | sed -n "$mid p"; } 2>&1 )
    velocity=$(echo "${var}" | head -n1)
    time=$(echo "${var}" | tail -n1)
    echo -e ${ncells}"\t"${velocity}"\t"${time}
}

export -f iterate

binrange=( $(seq 1001 1000 10000) $(seq 10001 10000 100000) $(seq 100001 300000 1000001))

# Set -j to increase thread count for parrallel execution
parallel -j 1 iterate ::: "${binrange[@]}"
