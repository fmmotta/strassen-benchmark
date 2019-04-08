#!/bin/bash

counter=1
max_size=10000000000000000
rm time_results.csv
touch time_results.csv
echo "n, time_strassen, time_blas" >> time_results.csv    
while [ $counter -le $max_size ]
do    
    ./time_test $counter 64 >> time_results.csv 
    counter=$(( $counter * 2 )) 
done

gnuplot plot_times
