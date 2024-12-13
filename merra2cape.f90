program gcm2cape

! compile with Makefile

! first version (gcm2cape.f90) by Alex Koch in early 2022
! updates by Jed Kaplan in late 2024 for MERRA-2 reanalysis

use parametersmod, only : sp,dp,i2
use netcdf
use capemod,       only : getcape
use calendarmod,   only : timestruct,ymdt2jd

implicit none

! parameters

real(sp), parameter :: sat_pressure_0c = 6.112 ! millibar
real(sp), parameter :: molweight_ratio = 0.622 ! molecular weight ratio water vapour to dry air; unitless

real(sp), parameter :: Tfreeze  =   273.15
real(sp), parameter :: rmissing = -9999.

real(sp), parameter :: minq = tiny(1._sp)

! input variables

real(dp), allocatable, dimension(:)     :: lon
real(dp), allocatable, dimension(:)     :: lat    
real(sp), allocatable, dimension(:)     :: lev    ! air pressure (hPa)
real(dp), allocatable, dimension(:)     :: time

real(sp), allocatable, dimension(:,:,:) :: Tk ! air temperature (K)
real(sp), allocatable, dimension(:,:,:) :: QV ! specific humidity (kg kg-1)

real(sp) :: qmissing

integer(i2) :: imissing

! intermediate variables

real(sp), allocatable, dimension(:) :: Tair              ! air temperature (degC)
real(sp), allocatable, dimension(:) :: Tdew              ! dewpoint temperature (degC)

! output variables

real(sp), allocatable, dimension(:,:) :: cape
real(sp), allocatable, dimension(:,:) :: cin

integer(i2), allocatable, dimension(:,:) :: ovar

! other local variables

character(200) :: infile
character(200) :: outfile

character(60)  :: status_line

real(sp) :: q     ! specific humidity (kg kg-1)
real(sp) :: p     ! pressure (hPa)

integer :: status
integer :: ifid
integer :: ofid
integer :: dimid
integer :: varid

integer :: xlen
integer :: ylen
integer :: tlen
integer :: nlev

integer :: x
integer :: y
integer :: l
integer :: t
integer :: l0
integer :: lm
integer :: ot
integer :: otb

integer :: tairid
integer :: humsid

integer :: capeid
integer :: cinid

real(sp) :: varmin
real(sp) :: varmax

real(sp), dimension(2) :: ar_cape
real(sp), dimension(2) :: ar_cin

type(timestruct) :: bt
type(timestruct) :: tt

integer :: year
integer :: mon
integer :: day
integer :: hr

real(sp) :: scale_factor
real(sp) :: add_offset

character(50) :: yyyymmdd

! -----------------------------------------------------------------------------------------------

call getarg(1,infile)

call getarg(2,outfile)

! ---

status = nf90_open(infile,nf90_nowrite,ifid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ifid,'lon',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=xlen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ifid,'lat',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=ylen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ifid,'lev',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=nlev)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ifid,'time',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=tlen)
if (status /= nf90_noerr) call handle_err(status)

! ---

allocate(lon(xlen))
allocate(lat(ylen))
allocate(lev(nlev))
allocate(time(tlen))

status = nf90_inq_varid(ifid,'lon',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,lon)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ifid,'lat',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,lat)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ifid,'lev',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,lev)
if (status /= nf90_noerr) call handle_err(status)

lm = minloc(lev,mask=lev >= 100.,dim=1)

! ---
! time relative to 1980-01-01 for time variable

bt = timestruct(1980,1,1,0,0,0.)

call ymdt2jd(bt)

call getarg(3,yyyymmdd)

read(yyyymmdd(1:4),*)year
read(yyyymmdd(5:6),*)mon
read(yyyymmdd(7:8),*)day

tt = timestruct(year,mon,day,0,0,0.)

call ymdt2jd(tt)

do t = 1,4

  hr = 6 * (t - 1)
  
  tt = timestruct(year,mon,day,hr,0,0.)
  
  call ymdt2jd(tt)

  time(t) = tt%jd - bt%jd

end do

! ---
! time index position in monthly output file

bt = timestruct(year,1,1,0,0,0.)

call ymdt2jd(bt)

tt = timestruct(year,mon,day,0,0,0.)

call ymdt2jd(tt)

otb = 4 * nint(tt%jd - bt%jd)

! write(0,*)otb,time

! ---

allocate(Tk(xlen,ylen,nlev))
allocate(QV(xlen,ylen,nlev))

allocate(Tair(nlev))
allocate(Tdew(nlev))

allocate(cape(xlen,ylen))
allocate(cin(xlen,ylen))

allocate(ovar(xlen,ylen))

! assign missing value to output variables

cape = rmissing
cin  = rmissing

! -------------------------------------------------
! get input varid

status = nf90_inq_varid(ifid,'T',tairid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ifid,tairid,'missing_value',qmissing)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ifid,'QV',humsid)
if (status /= nf90_noerr) call handle_err(status)

! -------------------------------------------------------
! open output and write coordinate variables

status = nf90_open(outfile,nf90_write,ofid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'lon',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,lon)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'lat',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,lat)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'time',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,time,start=[otb+1],count=[tlen])
if (status /= nf90_noerr) call handle_err(status)

! get varids for output

status = nf90_inq_varid(ofid,'CAPE',capeid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'CIN',cinid)
if (status /= nf90_noerr) call handle_err(status)

! get output missing, scale factor, and add offset

status = nf90_get_att(ofid,capeid,'missing_value',imissing)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ofid,capeid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ofid,capeid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

! -------------------------------------------------------
! main timestep loop, calculate one timestep at a time

do t = 1,tlen

  ot = otb + t
  
  ! ---

  status = nf90_get_var(ifid,tairid,Tk,start=[1,1,1,t],count=[xlen,ylen,nlev,1])   
  if (status /= nf90_noerr) call handle_err(status)
  
  status = nf90_get_var(ifid,humsid,QV,start=[1,1,1,t],count=[xlen,ylen,nlev,1])   
  if (status /= nf90_noerr) call handle_err(status)
  
  write(status_line,'(a,i0,a,i0,a,i0,a,i0,a,f0.3)')' working on ',year,'-',mon,'-',day,': ',ot,' ',time(t)
  call overprint(status_line)

  do y = 1,ylen
  
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,l0,Tair,Tdew,q,p)
    !$OMP DO SCHEDULE(GUIDED)

    do x = 1,xlen
    
      Tdew = rmissing
    
      do l = 1,lm

        ! there are places in the input where q = 0, which breaks the calculations, so skip these values
              
        if (QV(x,y,l) == qmissing .or. .not. QV(x,y,l) > 0.) cycle
    
        ! calculate dewpoint temperature
        
        Tair(l) = Tk(x,y,l) - tfreeze
        q       = QV(x,y,l)
        p       = lev(l)

        Tdew(l) = dewpoint(Tair(l),q,p)

      end do

      l0 = maxloc(lev,mask=Tdew /= rmissing,dim=1)

      call getcape(lev(l0:lm),Tair(l0:lm),Tdew(l0:lm),cape(x,y),cin(x,y))
      
    end do

    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

  end do

  ! ---
  ! CAPE
  
  status = nf90_get_att(ofid,capeid,'actual_range',ar_cape)
  if (status /= nf90_noerr) call handle_err(status)

  varmin = minval(cape,mask=cape/=rmissing)
  varmax = maxval(cape,mask=cape/=rmissing)

  ar_cape(1) = min(ar_cape(1),varmin)
  ar_cape(2) = max(ar_cape(2),varmax)

  status = nf90_put_att(ofid,capeid,'actual_range',ar_cape)
  if (status /= nf90_noerr) call handle_err(status)

  ! ---
  ! CIN
  
  status = nf90_get_att(ofid,cinid,'actual_range',ar_cin)
  if (status /= nf90_noerr) call handle_err(status)

  varmin = minval(cin,mask=cape/=rmissing)
  varmax = maxval(cin,mask=cape/=rmissing)

  ar_cin(1) = min(ar_cin(1),varmin)
  ar_cin(2) = max(ar_cin(2),varmax)

  status = nf90_put_att(ofid,cinid,'actual_range',ar_cin)
  if (status /= nf90_noerr) call handle_err(status)

  ! ---
  
  ovar = imissing

  where (cape /= rmissing) ovar = nint((cape + add_offset) / scale_factor,kind=i2)
  
  status = nf90_put_var(ofid,capeid,ovar,start=[1,1,ot],count=[xlen,ylen,1])
  if (status /= nf90_noerr) call handle_err(status)

  ovar = imissing

  where (cin /= rmissing) ovar = nint((cin + add_offset) / scale_factor,kind=i2)

  status = nf90_put_var(ofid,cinid,ovar,start=[1,1,ot],count=[xlen,ylen,1])
  if (status /= nf90_noerr) call handle_err(status)

end do

! write(0,*)
! write(0,*)'CAPE:',ar_cape
! write(0,*)'CIN: ',ar_cin

! -------------------------------------------------------

status = nf90_close(ifid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ofid)
if (status /= nf90_noerr) call handle_err(status)

! ----------------------------------------------------------------------------------------------------------------

contains

subroutine handle_err(status)

! Internal subroutine - checks error status after each netcdf call,
! prints out text message each time an error code is returned. 

integer, intent (in) :: status

if(status /= nf90_noerr) then 
  print *, trim(nf90_strerror(status))
  stop
end if

end subroutine handle_err

! ----------------------------------------------------------------------------------------------------------------

real(sp) function lvap(Tair)  ! (kJ kg-1)

! function to calculate the temperature-dependent latent heat of vaporization of water, based on:
! Henderson-Sellers, B. (1984). A new formula for latent heat of vaporization of water as a function of temperature.
! Quarterly Journal of the Royal Meteorological Society, 110(466), 1186-1190. doi:10.1002/qj.49711046626

use parametersmod, only : sp,tfreeze

implicit none

! argument

real(sp), intent(in) :: Tair  ! air temperature (degC)

! local variable

real(sp) :: Tk  ! air temperature (degC)

! ----

Tk = Tair + Tfreeze

lvap = 0.001 * 1.91846e6 * (Tk / (Tk - 33.91))**2  ! Henderson-Sellers (1984), eqn. 1

end function lvap

! ----------------------------------------------------------------------------------------------------------------

real(sp) function dewpoint(Tair,q,p)  ! (degC)

use parametersmod, only : sp

! arguments

real(sp), intent(in) :: Tair  ! air temperature (degC)
real(sp), intent(in) :: q     ! specific humidity (kg kg-1)
real(sp), intent(in) :: p     ! air pressure (hPa)

! parameters

real(sp), parameter :: e0 =   6.113  ! standard vapor pressure? (hPa)
real(sp), parameter :: T0 = 273.15   ! freezing point of water (K)
real(sp), parameter :: Rv = 461.     ! water vapor gas constant (J K-1 kg-1) 
real(sp), parameter :: eps =  0.622  ! molecular mass ratio of water vapor to dry air (g g-1)

real(sp), parameter :: T0i = 1. / T0

! local variables

real(sp) :: L  ! latent heat of vaporization of water (J kg-1)
real(sp) :: r  ! mixing ratio (kg water vapor kg dry air-1)
real(sp) :: e  ! vapor pressure (hPa)

! ----
! these equations are based on WMO equation 4.A.6 via the MetPy documentation

r = q / (1. - q)  

e = r / (eps + r) * p

! lvap formula in function above

L = lvap(Tair) * 1000.

! Stull formulation

dewpoint = (T0i - Rv / L * log(e / e0))**(-1) - T0

end function dewpoint

! ----------------------------------------------------------------------------------------------------------------

subroutine overprint(message)

use parametersmod, only : stderr

implicit none

! argument

character(*), intent(in) :: message

! parameter

character, parameter :: cr = char(13)

! ---

write(stderr,'(a)',advance='no')message
flush(0)
write(0,'(a1)',advance='no')cr

end subroutine overprint

! ----------------------------------------------------------------------------------------------------------------

end program gcm2cape
