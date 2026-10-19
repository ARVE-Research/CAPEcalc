#!/usr/bin/bash

for ((yr=2012;yr<=2023;yr++))
do
  for ((m=1;m<=12;m++))
  do
    
    mon=`printf %02i $m`
    
    daydir=rawdata/$yr$mon

    output=CAPE/${daydir##*/}.nc
    
    echo $output

    ncgen -4 -o $output cape.cdl
    
    for infile in `ls $daydir/*.nc4`
    do
      
      yyyymmdd=$(expr $infile : '.*\([0-9]\{8\}\)[^/]*/*$')

      ./merra2cape $infile $output $yyyymmdd
      
    done

    echo

  done
done
