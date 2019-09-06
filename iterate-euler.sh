#!/usr/bin/env bash

export SHELL=$(type -p bash)

function iterate {
    ncells=$1
    mid=$(( (ncells / 2) + 1))
    # Run algorithm - store time and midpoint in var
    var=$( TIMEFORMAT="%R"; { time ./euler_equations_commandline -n "${ncells}" -t 0.15 -d 1 -D 1 -v -2 -V 2 -p 0.4 -P 0.4 | cut -f 2 | sed -n "$mid p"; } 2>&1 )
  
    velocity=$(echo "${var}" | head -n1)
    time=$(echo "${var}" | tail -n1)
    echo -e ${ncells}"\t"${velocity}"\t"${time}
}

export -f iterate

#binrange=( $(seq 101 100 1000) $(seq 1001 1000 10000) $(seq 10001 10000 110000) $(seq 100001 100000 1100000) )
binrange=( 101 501 1001 5001 10001 50001 100001 500001 1000001 )
parallel -j 1 iterate ::: "${binrange[@]}"
