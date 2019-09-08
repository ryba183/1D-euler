#!/usr/bin/env bash

export SHELL=$(type -p bash)

function iterate {
    ncells=$1
    mid=$(( (ncells / 2) + 1))
    # Run algorithm - store time and midpoint in var
    var=$( TIMEFORMAT="%U"; { time ./1d_euler -n "${ncells}" -f | cut -f 3 | sed -n "$mid p"; } 2>&1 )
    velocity=$(echo "${var}" | head -n1)
    time=$(echo "${var}" | tail -n1)
    echo -e ${ncells}"\t"${velocity}"\t"${time}
}

export -f iterate

dat='velocity_test2/velocity_test2.tsv'
png="${dat%.*}".png

binrange=( $(seq 1001 1000 10000) $(seq 10001 10000 100000) $(seq 100001 300000 1000001))
# Set -j to increase thread count for parallel execution
parallel -j 1 iterate ::: "${binrange[@]}" >> "${dat}"

gnuplot << EOF
  set terminal png
  set output '${png}'
  set logscale xy 10
  set format xy "%e
  unset key
  set multiplot layout 2,1
    set xlabel 'cells'
    set ylabel 'velocity'
    plot '${dat}' u 1:(abs(\$2)) w l
    set ylabel 'time / seconds'
    set xlabel 'cells'
    plot '${dat}' u 1:3 w l
  unset multiplot
EOF


