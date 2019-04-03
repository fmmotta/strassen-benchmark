#!/bin/bash

counter=1
max_size=4097   
rm block_results.csv
touch block_results.csv
echo "n, time_strassen" >> block_results.csv    
while [ $counter -le $max_size ]
do    
    ./a.out 2000 $counter  >> block_results.csv 
    counter=$(( $counter * 2 )) 
done

gnuplot plot_min_block
