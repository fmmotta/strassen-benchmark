#!/bin/bash

counter=8192
max_size=2   
rm block_results.csv
touch block_results.csv
echo "n, time_strassen" >> block_results.csv    
while [ $counter -ge $max_size ]
do    
    ./find_opt_block 16384 $counter  >> block_results.csv 
    counter=$(( $counter / 2 )) 
done

gnuplot plot_min_block
