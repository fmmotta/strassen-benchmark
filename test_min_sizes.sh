#!/bin/bash

counter=16384
max_size=15000   
rm block_results.csv
touch block_results.csv
echo "n, time_strassen" >> block_results.csv    
while [ $counter -ge 2 ]
do    
    ./a.out 10000 $counter  >> block_results.csv 
    counter=$(( $counter / 2 )) 
done

gnuplot plot_min_block
