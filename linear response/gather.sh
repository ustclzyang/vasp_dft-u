#!/bin/bash
for v in +0.05 -0.05 +0.10 -0.10 +0.15 -0.15 +0.20 -0.20
do

printf "%s " "$v" >> mylog
awk -f gather.awk OUTCAR.V=$v.ICHARG=11 >> mylog
printf " " >> mylog
awk -f gather.awk OUTCAR.V=$v >> mylog
printf "\n" >> mylog

done
