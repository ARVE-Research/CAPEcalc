#!/usr/bin/env bash

# use this script to calculate CAPE and CIN from GCM output
# the script assumes that the GCM output is in the format of one file per month with data for all atmospheric levels
# the script generates an empty output file and then executes the CAPEcalc program to make the calculations

for ((yr=2001;yr<=2010;yr++))
do
  for ((m=1;m<=12;m++))
  do
    
    mon=`printf %02i $m`
    
    daydir=rawdata/$yr$mon

    output=CAPE/${daydir##*/}.nc
    
    echo $output

    ncgen -4 -o $output cape.cdl
    
    for infile in `ls $daydir/*$yr$mon*.nc4`
    do
      
      yyyymmdd=$(expr $infile : '.*\([0-9]\{8\}\)[^/]*/*$')

      ./CAPEcalc $infile $output $yyyymmdd
      
    done

    echo

  done
done
