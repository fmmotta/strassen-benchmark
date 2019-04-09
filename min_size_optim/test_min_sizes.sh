#!/bin/bash

counter=2048
max_size=2   
rm block_results.csv
touch block_results.csv
echo "n, time_strassen" >> block_results.csv    
while [ $counter -ge $max_size ]
do    
    ./min_block_search 4096 $counter  >> block_results.csv 
    counter=$(( $counter / 2 )) 
done

gnuplot plot_min_block
