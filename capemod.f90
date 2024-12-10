module capemod

! modifications to the code below done by Alex Koch in early 2022
! further edits and formatting by Jed Kaplan in late 2024

implicit none

public  :: getcape
private :: getqvs
private :: getqvi
private :: getthe

contains

! -----------------------------------------------------------------------------------------------
! the routine below is is based on
! -----------------------------------------------------------------------
! 
! getcape - a fortran90 subroutine to calculate Convective Available
!           Potential Energy (CAPE) from a sounding.
! 
! Version 1.02                Last modified:  10 October 2008
! 
! Author: George H. Bryan
!         Mesoscale and Microscale Meteorology Division
!         National Center for Atmospheric Research
!         Boulder, Colorado, USA
!         gbryan@ucar.edu
! 
! Disclaimer: This code is made available WITHOUT WARRANTY.
! 
! References: Bolton (1980, MWR, p. 1046) (constants and definitions)
!             Bryan and Fritsch (2004, MWR, p. 2421) (ice processes)
! 
! -----------------------------------------------------------------------
!
! Input:    nk - number of levels in the sounding (integer)
!
!         p_in - one-dimensional array of pressure (mb) (real)
!
!         t_in - one-dimensional array of temperature (C) (real)
!
!         d_in - one-dimensional array of dewpoint temperature (C) (real)
!
! Output: cape - Convective Available Potential Energy (J/kg) (real)
!
!          cin - Convective Inhibition (J/kg) (real)
!
! -----------------------------------------------------------------------
! 
!  ©2019 - University Corporation for Atmospheric Research
! 
!  Permission is hereby granted, free of charge, to any person obtaining a copy
!  of this software and associated documentation files (the "Software"), to deal
!  in the Software without restriction, including without limitation the rights
!  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!  copies of the Software, and to permit persons to whom the Software is furnished
!  to do so, subject to the following conditions:
! 
!  The above copyright notice and this permission notice shall be included in all
!  copies or substantial portions of the Software.
! 
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
!  THE SOFTWARE.
! 
! -----------------------------------------------------------------------

subroutine getcape(p_in,t_in,td_in,cape,cin)

use parametersmod, only : sp,r8
    
implicit none

! arguments

real(sp), dimension(:), intent(in)  :: p_in
real(sp), dimension(:), intent(in)  :: t_in
real(sp), dimension(:), intent(in)  :: td_in
real(sp),               intent(out) :: cape
real(sp),               intent(out) :: cin

! user-adjustable parameters

real(sp), parameter :: pinc = 100.  ! Pressure increment (Pa) 
                                    ! (smaller number yields more accurate results,
                                    !  larger number makes code go faster)

integer,  parameter :: source = 1   ! Source parcel: 
                                    ! 1 = surface
                                    ! 2 = most unstable (max theta-e)
                                    ! 3 = mixed-layer (specify ml_depth)

real(sp), parameter :: ml_depth =  200.  ! depth (m) of mixed layer for source = 3

integer,  parameter :: adiabat = 1  ! Formulation of moist adiabat:
                                    ! 1 = pseudoadiabatic, liquid only
                                    ! 2 = reversible, liquid only
                                    ! 3 = pseudoadiabatic, with ice
                                    ! 4 = reversible, with ice

! fixed parameters

real(sp), parameter :: g     = 9.81
real(sp), parameter :: p00   = 100000.
real(sp), parameter :: cp    = 1005.7
real(sp), parameter :: rd    = 287.04
real(sp), parameter :: rv    = 461.5
real(sp), parameter :: xlv   = 2501000.
real(sp), parameter :: xls   = 2836017.
real(sp), parameter :: t0    = 273.15
real(sp), parameter :: cpv   = 1875.
real(sp), parameter :: cpl   = 4190.
real(sp), parameter :: cpi   = 2118.636

real(sp), parameter :: lv1   = xlv + (cpl - cpv) * t0
real(sp), parameter :: lv2   = cpl - cpv
real(sp), parameter :: ls1   = xls + (cpi - cpv) * t0
real(sp), parameter :: ls2   = cpi - cpv

real(sp), parameter :: rp00  = 1. / p00
real(sp), parameter :: eps   = rd / rv
real(sp), parameter :: reps  = rv / rd
real(sp), parameter :: rddcp = rd / cp
real(sp), parameter :: cpdrd = cp / rd
real(sp), parameter :: cpdg  = cp / g

real(sp), parameter :: converge = 0.0002

integer, parameter :: debug_level = 200

! local variables

integer :: nk

logical :: doit
logical :: ice
logical :: cloud
logical :: not_converged
integer :: k
integer :: kmax
integer :: n
integer :: nloop
integer :: i
integer :: orec

real(sp), allocatable, dimension(:) :: p
real(sp), allocatable, dimension(:) :: t
real(sp), allocatable, dimension(:) :: td
real(sp), allocatable, dimension(:) :: pi
real(sp), allocatable, dimension(:) :: q
real(sp), allocatable, dimension(:) :: th
real(sp), allocatable, dimension(:) :: thv
real(sp), allocatable, dimension(:) :: z
real(sp), allocatable, dimension(:) :: pt
real(sp), allocatable, dimension(:) :: pb
real(sp), allocatable, dimension(:) :: pc
real(sp), allocatable, dimension(:) :: pn
real(sp), allocatable, dimension(:) :: ptv

real(sp) :: the
real(sp) :: maxthe
real(sp) :: parea
real(sp) :: narea
real(sp) :: lfc
real(sp) :: th1
real(sp) :: p1
real(sp) :: t1
real(sp) :: qv1
real(sp) :: ql1
real(sp) :: qi1
real(sp) :: b1
real(sp) :: pi1
real(sp) :: thv1
real(sp) :: qt
real(sp) :: dp
real(sp) :: dz
real(sp) :: ps
real(sp) :: frac
real(sp) :: th2
real(sp) :: p2
real(sp) :: t2
real(sp) :: qv2
real(sp) :: ql2
real(sp) :: qi2
real(sp) :: b2
real(sp) :: pi2
real(sp) :: thv2
real(sp) :: thlast
real(sp) :: fliq
real(sp) :: fice
real(sp) :: tbar
real(sp) :: qvbar
real(sp) :: qlbar
real(sp) :: qibar
real(sp) :: lhv
real(sp) :: lhs
real(sp) :: lhf
real(sp) :: rm
real(sp) :: cpm
real(r8) :: avgth
real(r8) :: avgqv

! real(sp) :: getqvs
! real(sp) :: getqvi
! real(sp) :: getthe

! -----------------------------------------------------------------------

nk = size(p_in)

allocate(p(nk))
allocate(t(nk))
allocate(td(nk))
allocate(pi(nk))
allocate(q(nk))
allocate(th(nk))
allocate(thv(nk))
allocate(z(nk))
allocate(pt(nk))
allocate(pb(nk))
allocate(pc(nk))
allocate(pn(nk))
allocate(ptv(nk))

! --------------

! convert p,t,td to mks units; get pi,q,th,thv 

do k=1,nk

  p(k)   = 100. * p_in(k)            ! convert pressure from kPa to Pa
  t(k)   = t0 + t_in(k)              ! convert temperature to K
  td(k)  = t0 + td_in(k)             ! convert dewpoint to K

  pi(k)  = (p(k) * rp00)**rddcp
  q(k)   = getqvs(p(k),td(k))
  th(k)  = t(k )/ pi(k)
  thv(k) = th(k) * (1. + reps * q(k)) / (1. + q(k))

end do

! get height using the hydrostatic equation

z(1) = 0.

do k = 2,nk

  dz = -cpdg * 0.5 * (thv(k) + thv(k-1)) * (pi(k) - pi(k-1))
  z(k) = z(k-1) + dz

end do

! find source parcel

select case(source)
case(1)

  ! use surface parcel

  kmax = 1

case(2)

  ! use most unstable parcel (max theta-e)
  
  if (p(1) < 50000.) then

    !  first report is above 500 mb ... just use the first level reported

    kmax = 1
    maxthe = getthe(p(1),t(1),td(1),q(1))
    
  else

    ! find max thetae below 500 mb
    maxthe = 0.

    do k = 1,nk

      if (p(k) >= 50000.) then
      
        the = getthe(p(k),t(k),td(k),q(k))

        if (the > maxthe ) then
          maxthe = the
          kmax = k

        end if
      end if
    end do
  end if

  if (debug_level >= 100) write(0,*)'  kmax,maxthe = ',kmax,maxthe

case(3)

  ! use mixed layer

  if (z(2) - z(1) > ml_depth ) then

    ! the second level is above the mixed-layer depth:  just use the lowest level

    avgth = th(1)
    avgqv = q(1)
    kmax = 1
  
  else if (z(nk) < ml_depth ) then

    ! the top-most level is within the mixed layer:  just use the upper-most level

    avgth = th(nk)
    avgqv = q(nk)
    kmax = nk
    
  else 
  
    !  calculate the mixed-layer properties:

    avgth = 0.
    avgqv = 0.
    k = 2

    if (debug_level >= 100) write(0,*)'  ml_depth = ',ml_depth
    if (debug_level >= 100) write(0,*)'  k,z,th,q:'
    if (debug_level >= 100) write(0,*)1,z(1),th(1),q(1)

    do while (z(k) <= ml_depth .and. k <= nk)

      if (debug_level >= 100) write(0,*)k,z(k),th(k),q(k)

      avgth = avgth + 0.5 * (z(k) - z(k-1)) * (th(k) + th(k-1))
      avgqv = avgqv + 0.5 * (z(k) - z(k-1)) * (q(k)  + q(k-1))

      k = k + 1

    end do

    th2 = th(k-1) + (th(k) - th(k-1)) * (ml_depth - z(k-1)) / (z(k) - z(k-1))
    qv2 =  q(k-1) + (q(k)  -  q(k-1)) * (ml_depth - z(k-1)) / (z(k) - z(k-1))

    if (debug_level >= 100) write(0,*)999,ml_depth,th2,qv2

    avgth = avgth + 0.5 * (ml_depth - z(k-1)) * (th2 + th(k-1))
    avgqv = avgqv + 0.5 * (ml_depth - z(k-1)) * (qv2 +  q(k-1))

    if (debug_level >= 100) write(0,*)k,z(k),th(k),q(k)

    avgth = avgth / ml_depth
    avgqv = avgqv / ml_depth

    kmax = 1

  end if

  if (debug_level >= 100) write(0,*)avgth,avgqv
  
case default

  write(0,*)
  write(0,*)'  Invalid selection for source parcel'
  write(0,*)'  must be 1, 2, or 3'
  write(0,*)'  source = ',source
  write(0,*)
  stop

end select

! define parcel properties at initial location

narea = 0.

if (source == 1 .or. source == 2) then

  k    = kmax
  th2  = th(kmax)
  pi2  = pi(kmax)
  p2   = p(kmax)
  t2   = t(kmax)
  thv2 = thv(kmax)
  qv2  = q(kmax)
  b2   = 0.

else

  k    = kmax
  th2  = avgth
  qv2  = avgqv
  thv2 = th2 * (1. + reps * qv2) / (1. + qv2)
  pi2  = pi(kmax)
  p2   = p(kmax)
  t2   = th2 * pi2
  b2   = g * (thv2 - thv(kmax)) / thv(kmax)

end if

ql2 = 0.
qi2 = 0.
qt  = qv2

cape = 0.
cin  = 0.
lfc  = 0.

doit = .true.

cloud = .false.

if (adiabat == 1 .or. adiabat == 2) then
  ice = .false.
else
  ice = .true.
end if

the = getthe(p2,t2,t2,qv2)

if (debug_level >= 100) write(0,*)'  the = ',the

! begin ascent of parcel

if (debug_level >= 100) then

  write(0,*)'  Start loop:'
  write(0,*)'  p2,th2,qv2 = ',p2,th2,qv2

end if

do while (doit .and. k < nk)

  k  = k + 1
  b1 = b2
  
  dp = p(k-1) - p(k)
  
  if (dp < pinc) then
  
    nloop = 1
  
  else
  
    nloop = 1 + int(dp / pinc)
  
    dp = dp / real(nloop)
  
  end if

  do n = 1,nloop

    p1   =   p2
    t1   =   t2
    pi1  =  pi2
    th1  =  th2
    qv1  =  qv2
    ql1  =  ql2
    qi1  =  qi2
    thv1 = thv2

    p2 = p2 - dp

    pi2 = (p2 * rp00)**rddcp

    thlast = th1

    i = 0

    not_converged = .true.

    do while (not_converged)

      i = i + 1

      t2 = thlast * pi2

      if (ice) then

        fliq = max(min((t2 - 233.15) / (t0 - 233.15),1.),0.)
        fice = 1. - fliq

      else

        fliq = 1.
        fice = 0.

      end if

      qv2 = min(qt,fliq * getqvs(p2,t2) + fice * getqvi(p2,t2))

      qi2 = max(fice * (qt - qv2),0.)

      ql2 = max(qt - qv2 - qi2,0.)

      tbar  = 0.5 * (t1 + t2)
      qvbar = 0.5 * (qv1 + qv2)
      qlbar = 0.5 * (ql1 + ql2)
      qibar = 0.5 * (qi1 + qi2)

      lhv = lv1 - lv2 * tbar
      lhs = ls1 - ls2 * tbar
      lhf = lhs - lhv

      rm = rd + rv * qvbar

      cpm = cp + cpv * qvbar + cpl * qlbar + cpi * qibar

      th2 = th1 * exp(lhv * (ql2 - ql1) / (cpm * tbar) + & 
                  lhs * (qi2 - qi1) / (cpm * tbar)     + &
                 (rm / cpm - rd / cp) * alog(p2 / p1))

      if (i > 90) write(0,*)i,th2,thlast,th2-thlast

      if (i > 100) then

        write(0,*)
        write(0,*)'  Error:  lack of convergence'
        write(0,*)
        write(0,*)'  ... stopping iteration '
        write(0,*)
        stop 1001

      end if
      
      if (abs(th2 - thlast) > converge) then

        thlast = thlast + 0.3 * (th2 - thlast)

      else

        not_converged = .false.

      end if

    end do

    ! Latest pressure increment is complete.  Calculate some important stuff:

    if (ql2 >= 1.e-10) cloud = .true.

    if (adiabat == 1 .or. adiabat == 3) then

      ! pseudoadiabat

      qt  = qv2
      ql2 = 0.
      qi2 = 0.

    else if (adiabat <= 0 .or. adiabat >= 5) then

      write(0,*)
      write(0,*)'  Undefined adiabat'
      write(0,*)
      stop 10000

    end if

  end do

  thv2 = th2 * (1. + reps * qv2) / (1. + qv2 + ql2 + qi2)

  b2 = g * (thv2 - thv(k)) / thv(k)

  dz = -cpdg * 0.5 * (thv(k) + thv(k-1)) * (pi(k) - pi(k-1))

  the = getthe(p2,t2,t2,qv2)

  ! Get contributions to CAPE and CIN:

  if (b2 >= 0. .and. b1 < 0.) then

    ! first trip into positive area

    ps = p(k-1)+  (p(k) - p(k-1)) * -b1 / (b2 - b1)

    frac = b2 / (b2 - b1)

    parea = 0.5 * b2 * dz * frac

    narea = narea - 0.5 * b1 * dz * (1. - frac)

    if (debug_level >= 200) then
      write(0,*)'      b1,b2 = ',b1,b2
      write(0,*)'      p1,ps,p2 = ',p(k-1),ps,p(k)
      write(0,*)'      frac = ',frac
      write(0,*)'      parea = ',parea
      write(0,*)'      narea = ',narea
    end if
    
    cin  = cin  + narea
    narea = 0.

  else if (b2 < 0. .and. b1 > 0.) then

    ! first trip into neg area

    ps = p(k-1) + (p(k) - p(k-1)) * -b1 / (b2 - b1)

    frac = b1 / (b1 - b2)

    parea =  0.5 * b1 * dz * frac
    narea = -0.5 * b2 * dz * (1. - frac)

    if (debug_level >= 200) then
      write(0,*)'      b1,b2 = ',b1,b2
      write(0,*)'      p1,ps,p2 = ',p(k-1),ps,p(k)
      write(0,*)'      frac = ',frac
      write(0,*)'      parea = ',parea
      write(0,*)'      narea = ',narea
    end if
    
    else if (b2 < 0.) then

      !  still collecting negative buoyancy

      parea = 0.
      narea = narea - 0.5 * dz * (b1 + b2)

    else

      !  still collecting positive buoyancy

      parea = 0.5 * dz * (b1 + b2)
      narea = 0.

    end if

    cape = cape + max(0.,parea)

    if (debug_level >= 200) then
      write(6,102) p2,b1,b2,cape,cin,cloud
102   format(5(f13.4),2x,l1)
    end if

    if (p(k) <= 10000. .and. b2 < 0.) then

    ! stop if b < 0 and p < 100 mb

    doit = .false.

  end if

end do

end subroutine getcape

! -----------------------------------------------------------------------

real(sp) function getqvs(p,t)

use parametersmod, only : sp,tfreeze

implicit none

! arguments

real(sp), intent(in) :: p
real(sp), intent(in) :: t

! parameter

real(sp), parameter :: eps = 287.04 / 461.5

! local variable

real(sp) :: es

! ----

es = 611.2 * exp(17.67 * (tfreeze) / (t - 29.65))

getqvs = eps * es / (p - es)

end function getqvs

! -----------------------------------------------------------------------

real(sp) function getqvi(p,t)

use parametersmod, only : sp,tfreeze

implicit none

! arguments

real(sp), intent(in) :: p
real(sp), intent(in) :: t

! parameter

real(sp), parameter :: eps = 287.04 / 461.5

! local variable

real(sp) :: es

! ----

es = 611.2 * exp(21.8745584 * (t - tfreeze) / (t - 7.66))

getqvi = eps * es / (p - es)

end function getqvi

! -----------------------------------------------------------------------

real(sp) function getthe(p,t,td,q)

use parametersmod, only : sp

implicit none

! arguments

real(sp), intent(in) :: p
real(sp), intent(in) :: t
real(sp), intent(in) :: td
real(sp), intent(in) :: q

! local variable

real :: tlcl

! ----

if (td - t >= -0.1) then

  tlcl = t

else

  tlcl = 56. + ((td - 56.)**(-1) + 0.00125 * alog(t / td))**(-1)

end if

getthe = t * ((100000. / p)**(0.2854 * (1. - 0.28 * q))) * exp(((3376. / tlcl) - 2.54) * q * (1. + 0.81 * q))

end function getthe

! -----------------------------------------------------------------------------------------------

end module capemod
