#!/bin/bash

counter=1
rm time_results.csv
touch time_results.csv
echo "n, time_strassen, time_blas" >> time_results.csv    
while [ $counter -le 1000 ]
do    
    ./a.out $counter 1 >> time_results.csv 
    counter=$(( $counter + 1 )) 
done

gnuplot plot_times
