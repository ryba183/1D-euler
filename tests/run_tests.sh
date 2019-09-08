#!/usr/bin/env bash

for boundary in transmissive periodic reflective; do

  dat=tests/test5-"${boundary}".dat
  gif="${dat%.*}".gif

  ./1d_euler --density_left 5.99924 \
             --velocity_left 19.5975 \
             --pressure_left 460.894 \
             --density_right 5.99242 \
             --velocity_right -6.19633 \
             --pressure_right 46.0950 \
             --time 0.15 \
             --boundary "${boundary}" \
    > "${dat}"

  gnuplot << EOF
  set terminal gif animate delay 3
  set output '${gif}'
  stats '${dat}' nooutput
  set xlabel 'position'
  set xrange [0:1]
  unset key
  do for [i=0:int(STATS_blocks)-2] {
    set multiplot layout 2,2
      set yrange [0:25]
      set ylabel 'density'
      plot '${dat}' u 1:2 index i w l
      set yrange [-20:25]
      set ylabel 'velocity'
      plot '${dat}' u 1:3 index i w l
      set yrange [0:2000]
      set ylabel 'pressure'
      plot '${dat}' u 1:4 index i w l
      set yrange [0:400]
      set ylabel 'internal energy'
      plot '${dat}' u 1:5 index i w l
    unset multiplot
  }
EOF

done


















