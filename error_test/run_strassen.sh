#!/bin/bash

counter=1
max_size=10000000000000000
	rm residue.csv
	touch residue.csv
	echo "n, residue" >> residue.csv    
	while [ $counter -le $max_size ]
	do    
	    ./error_test $counter 64 >> residue.csv 
    counter=$(( $counter * 2 )) 
done

gnuplot plot_times
