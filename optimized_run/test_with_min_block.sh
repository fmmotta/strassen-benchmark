#!/bin/bash

counter=1
max_size=10000000000000000
rm opt_block_results.csv
touch opt_block_results.csv
echo "n, time_strassen, time_blas" >> opt_block_results.csv    
while [ $counter -le $max_size ]
do    
    ./time_test $counter 129 >> opt_block_results.csv 
    counter=$(( $counter * 2 )) 
done

gnuplot plot_times
