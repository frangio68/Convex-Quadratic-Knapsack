#!/bin/bash
# requires the "generatore" random generator of Sensor Placement Problems
# that can bo found at http://www.di.unipi.it/optimize/Data/RDR.html

RANDOM=1         # Same seed for RANDOM...

for nitem in 10000 100000; do
for maxval in 100 1000; do
for seed in `seq 1 10`; do 

echo $nitem >testname.prm
echo $maxval >>testname.prm
echo $((RANDOM%1000+0)) >>testname.prm

./generatore testname

mv testname.rdr $nitem,$maxval,$seed.rdr

CQKnPSolve $nitem,$maxval,$seed.rdr

done
done
done


