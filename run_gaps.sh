#!/bin/bash	
export LC_NUMERIC="en_US.UTF-8"      #sistema le virgole dei decimali
                                     #i file vanno ottenuti compilando almeno una volta compiler.sh

gfortran gaps.f90 davidson.guess.lib2.o lapack_tebd.lib.f -O2 -o gaps.x   #compila linkando i file necessari

start=0.0
end=1.5
spacing=0.1

output_file="gaps_22.txt"                         #mettere L prima di ogni run (gaps_L.txt), segnarsi se obc, pbc
: > "$output_file"                                #pulisce il file
                          
for g in $(seq $start $spacing $end); do          #loop per i valori del campo g desiderati 
    echo "Sending input g= $g to simulation" 
    echo "$g" | ./gaps.x >> "$output_file"        #pipe
done
echo "Done simulating"
