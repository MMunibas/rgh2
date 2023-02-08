module rep_ker
implicit none
real*8, parameter :: dk26f1 = 1.0d0/14.0d0, dk26f2 = 1.0d0/18.0d0, &
dk24f1 = 2.0d0/15.0d0, dk24f2 = 2.0d0/21.0d0, dk25f1 = 2.0d0/21.0d0,&
dk25f2 = 1.0d0/14.0d0, akf1=2.0d0/3.0d0

contains

function drker24(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: drker24, xl, xs

xl = x
xs = xi
if (x .lt. xi) then
  xl = xi
  xs = x
end if

drker24 = dk24f1/xl**5 - dk24f2*xs/xl**6

end function drker24

function ddrker24(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: ddrker24, xl, xs

if (x .lt. xi) then
  ddrker24 = -dk24f2/xi**6
else
  ddrker24 = -5.0d0*dk24f1/x**6 + 6.0d0*dk24f2*xi/x**7
end if

end function ddrker24

function drker25(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: drker25, xl, xs

xl = x
xs = xi
if (x .lt. xi) then
  xl = xi
  xs = x
end if

drker25 = dk25f1/xl**6 - dk25f2*xs/xl**7

end function drker25
function ddrker25(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: ddrker25

if (x .lt. xi) then
  ddrker25=-dk25f2/xi**7
 else
  ddrker25=-6.0d0*dk25f1/x**7+7.0d0*dk25f2*xi/x**8
end if

end function ddrker25

function drker26(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: drker26, xl, xs

xl = x
xs = xi
if (x .lt. xi) then
  xl = xi
  xs = x
end if

drker26 = dk26f1/xl**7 - dk26f2*xs/xl**8

end function drker26

function ddrker26(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: ddrker26

if (x .lt. xi) then
  ddrker26 = -dk26f2/xi**8
else
  ddrker26 = -7.0d0*dk26f1/x**8 + 8.0d0*dk26f2*xi/x**9
end if

end function ddrker26
function atker23(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: atker23, xl, xs

xl = x
xs = xi
if (x .lt. xi) then
  xl = xi
  xs = x
end if

atker23 = 1.0d0 + xs*xl + 2.0d0*xs**2*xl  - akf1*xs**3

end function atker23

end module rep_ker





module pess
implicit none
real*8, allocatable, dimension(:,:) :: dataarray, dataarray1
integer :: nd1, nd2
contains

subroutine calcener(capr, smlr, theta, ener)
use rep_ker
use rkhs   
!use arrays
!use param!RKHS interpolation module
implicit none
real*8 :: lambda
integer :: ii, n, ios
real*8, parameter :: pi = acos(-1.0d0), piby180 = pi/180.0d0
real*8 :: rbc,d1, d2, d3,d4,h2ener, zero, db,dc, vabc, vtot, vab, vbc, vac, rab, rac, bc, derab, derac, derbc
real*8, intent(out) :: ener
real*8, intent(in) :: capr, smlr, theta
real*8, dimension(:) :: x(3)
!real*8, parameter :: pi =acos(-1.0d0)
real*8, parameter :: ma = 19.9924401762d0, mb = 1.00782503223d0, mc = 1.00782503223d0
real*8, parameter :: cb =  mc/(mb+mc), cc = mb/(mb+mc), cp=epsilon(1.0d0)
integer ::  nii
logical, save :: stored = .false.
type(kernel), save  :: pes1  ! The kernel type is needed to set up and evaluate a RKHS model
logical, save :: kread = .false.
logical, save :: ker1 = .false.

if (.not. ker1) then
  inquire(file="pes1.kernel", exist=ker1)  !file_exists will be true if the file exists and false otherwise
end if

lambda=0.1d-19
if (.not. stored ) then
 open(unit=11, file="asymph2p.coeff", status="old")
 open(unit=9, file="asympnehp.coeff", status="old")
 
 nd1=55 !h2+
 allocate(dataarray(nd1,2))
 do nii=1,nd1
   read(11,*)dataarray(nii,1),dataarray(nii,2)
 end do
 
 nd2=52 !neh+
 allocate(dataarray1(nd2,2))
 do nii=1,nd2
   read(9,*)dataarray1(nii,1),dataarray1(nii,2)
 end do
stored = .true.

end if

if (.not. kread) then
  if (ker1) then
    call pes1%load_from_file("pes1.kernel")
    print*,"data read from pes1.kernel"
    kread = .true.
  else



 n=38220
zero= -128.859837014476d0-0.499994784584d0 !Ne + (H+) + H atomic energies
    do ii = 1, n
        open(unit=102,file='pes.dat',status='old',iostat=ios)
        open(unit=104,file='pes1.csv',iostat=ios)
        read(102,*,iostat = ios) d1, d2, d3, d4
db = d3*cb
dc = d3*cc
d1=d1*piby180
rab = sqrt(abs(d2**2+db**2-2.0d0*d2*db*cos(pi-d1)))
rac = sqrt(abs(d2**2+dc**2-2.0d0*d2*dc*cos(d1)))
rbc = d3 !smlr h2+

vabc=0.0d0
vtot=0.0d0
vab=0.0d0
vbc=0.0d0
vac=0.0d0

call hehpdpot(rab,vab, derab)
call h2pdpot(rbc,vbc, derbc) 
call hehpdpot(rac,vac, derac)

if ((d4-zero) > 0.015d0) then
 write (104,*) (1.0d0-cos(d1))/2.0d0,",", d2, ",", d3, ",", "NaN"
 else
!write(104,*)(1.0d0-cos(d1))/2.0d0,",",d2,",",d3,",", d4, d4-vac-vab-vbc-zero, d4-(d4-vac-vab-vbc+zero),vab, vbc, vac
write(104,*)(1.0d0-cos(d1))/2.0d0,",",d2,",",d3,",", d4-vbc-vab-vac-zero !,rab, vab,rbc, vbc,rac, vac   !, d4-(d4-vac-vab-vbc+zero),vab, vbc, vac
end if
end do
 close(102)
 close(104)

    
    call pes1%read_grid("pes1.csv")
    print*,"data read from pes1.csv"
    call pes1%k1d(1)%init(TAYLOR_SPLINE_N2_KERNEL)         ! choose one-dimensional kernel for dimension 1
    call pes1%k1d(2)%init(RECIPROCAL_POWER_N2_M4_KERNEL)   ! choose one-dimensional kernel for dimension 2
    call pes1%k1d(3)%init(RECIPROCAL_POWER_N2_M4_KERNEL)

!    call pes1%calculate_coefficients_slow(lambda)
    call pes1%calculate_coefficients_fast()

    call pes1%calculate_sums()

    call pes1%save_to_file("pes1.kernel")
!
    kread = .true.
  end if
end if

x(1)=(1.0d0-cos(theta))/2.0d0
x(2)=capr
x(3)=smlr

ener = 0.0d0
call pes1%evaluate_fast(x,ener)

return

end subroutine calcener


subroutine h2pdpot(r,ener,der)
use rep_ker
!use arrays
implicit none
real*8, intent(in) :: r
real*8, intent(out) :: ener, der
integer :: jj

 ener=0.0d0
 der=0.0d0
 do jj = 1, nd1
   ener = ener + drker24(r,dataarray(jj,1))*dataarray(jj,2)
   der = der+ ddrker24(r,dataarray(jj,1))*dataarray(jj,2)
 end do

return

end subroutine h2pdpot

subroutine hehpdpot(r,ener,der)
use rep_ker
!use arrays
implicit none
real*8, intent(in) :: r
real*8, intent(out) :: ener, der
integer :: kk

 ener=0.0d0
 der=0.0d0
 do kk = 1, nd2
   ener = ener + drker24(r,dataarray1(kk,1))*dataarray1(kk,2)
   der = der+ ddrker24(r,dataarray1(kk,1))*dataarray1(kk,2)
 end do

return

end subroutine hehpdpot


end module pess

module arrays

implicit none
real*8, allocatable, dimension(:,:) :: dataarray , dataarray1
integer :: nd1,nd2
end module





!This module contains a subroutine which calculate the
!potential energies for He----H2+ system
module heh2pes
use arrays
use pess
contains 

subroutine tot_pes(capr,smlr,theta,vtot)
!!!!input
!theta in radian
!R in a.u.
!r in a.u.
!!!output
!energy in hartree 
!zero is set to the energy of E(Ne)+E(H)+E(H+)
use pess
use arrays
implicit none
real*8, intent (in) :: capr,smlr,theta
real*8, intent (out) :: vtot
real*8, parameter :: pi =acos(-1.0d0)
real*8, parameter :: ma = 19.9924401762d0 , mb = 1.00782503223d0, mc = 1.00782503223d0
real*8, parameter :: cb =  mc/(mb+mc), cc = mb/(mb+mc), cp=epsilon(1.0d0)
real*8 :: rab, rbc, rac, vabc, vab, vbc, vac, db, dc, sf, vl, angle, derab, derbc, derac
real*8, parameter :: piby180 = pi/180.0d0
integer ::  nii
vl=0.0d0
vabc=0.0d0
vtot=0.0d0
vab=0.0d0
vbc=0.0d0
vac=0.0d0
db = smlr*cb
dc = smlr*cc
rab = sqrt(abs(capr**2+db**2-2.0d0*capr*db*cos(pi-theta)))
rac = sqrt(abs(capr**2+dc**2-2.0d0*capr*dc*cos(theta)))
rbc = smlr

angle=theta
if (angle>pi/2.0d0)angle=pi-angle
!  print*, capr, smlr, angle! , angle*piby180
call calcener(capr, smlr, angle,vabc)


call hehpdpot(rab,vab, derab)
call h2pdpot(rbc,vbc, derbc) !HH
call hehpdpot(rac,vac, derac)



  vtot=vabc+vbc+vab+vac

return

end subroutine tot_pes

end module
