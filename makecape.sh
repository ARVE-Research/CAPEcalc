#!/usr/bin/bash

daydir=${1}

output=CAPE/${daydir##*/}.nc

echo $output

ncgen -4 -o $output cape.cdl

for infile in `ls $daydir/*1990010*.nc4`
do
  
  yyyymmdd=$(expr $infile : '.*\([0-9]\{8\}\)[^/]*/*$')
    
  ./merra2cape $infile $output $yyyymmdd
  
done

echo