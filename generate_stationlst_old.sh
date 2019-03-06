#!/bin/bash
seeddir=/work2/leon/COR_ALASKA_dir/seed_data
lstdir=/work2/leon/benchmark_old/station_old.lst
if [ -e station.lst ]; then
	        rm -f station.lst
fi
cd $seeddir
rm s1.lst
ls AKALL_*OCT.11.*seed >> seed_fname.lst
ls AKALL_*NOV.10.*seed >> seed_fname.lst
#ls AKALL_2017*seed >> seed_fname.lst
#ls AKALL_2018*seed >> seed_fname.lst
#ls AKALL_2019*seed >> seed_fname.lst
for sname in `cat seed_fname.lst`
do
	rdseed -S -f $sname
        cat rdseed.stations >> s1.lst
done
sort -k1,1 -u s1.lst > s2.lst
awk '{print $1,$4,$3}' s2.lst > station.lst
rm s1.lst s2.lst seed_fname.lst
cp station.lst $lstdir
