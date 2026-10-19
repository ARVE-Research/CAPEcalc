# CAPEcalc

The CAPEcalc program takes air temperature and specific humidity from an atmosphere model and estimates convective available potential energy (CAPE) and convective inhibition (CIN), which are two measures of the potential for thunderstorm activity. The program requires 4D fields (x,y,vertical level) of air temperature (K) and specific humidity (kg kg-1) and writes CAPE and CIN to a netCDF file that must be pre-generated before running the program.

The program is currently set up to work on a single GCM output file at a time that contains the temperature and humidity variables, it likewise writes a monthy file with CAPE and CIN. These files can be aggregated to a timeseries using tools such as `NCO` and `CDO`.

## Usage

1. Edit the netCDF template file for the calculation's output `cape.cdl` to reflect the source, temporal, and spatial resolution of the data.

2. Edit the wrapper script `CAPEcalc.sh` for the pathnames, filename pattern, and years of calculation for your input files.

3. Edit the `Makefile` with your preferred Fortran compiler and compiler options and the path to your netCDF-Fortran libraries

4. Run `CAPEcalc.sh`, which generates empty output files using `ncgen` and runs the Fortran executable `CAPEcalc` over the years/months specified in the script.
