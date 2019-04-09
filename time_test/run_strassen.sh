#!/bin/bash

counter=1
max_size=10000000000000000
rm time_results_no_opt.csv
touch time_results_no_opt.csv
echo "n, time_strassen, time_blas" >> time_results_no_opt.csv    
while [ $counter -le $max_size ]
do    
    ./time_test $counter 0 >> time_results_no_opt.csv 
    counter=$(( $counter * 2 )) 
done

gnuplot plot_times
