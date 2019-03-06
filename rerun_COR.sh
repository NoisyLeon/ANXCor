#!/bin/bash


cppexe=/home/leon/code/ANXCor/ANXCor
#cd /lustre/janus_scratch/life9360/COR_TEST_forye_Juan_002
for i in `seq 1 100`;
        do
                echo $i >> log_anxcorr
		export OMP_NUM_THREADS=12; $cppexe parameters.txt <<- END
		Y
		END
        done
#
