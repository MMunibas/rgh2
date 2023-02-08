program surfNeH2 
use heh2pes
implicit none
real*8, parameter :: pi = acos(-1.0d0)
real*8 :: thetarad
real*8 :: theta, capr, smlr, V, zero
integer :: i,j,k,n
theta=30.0d0
capr=3.2d0
smlr=2.084d0

thetarad=pi/180.0d0*theta
CALL tot_pes(CAPR,SMLR,thetarad,V)
write(19,*)  theta,capr,smlr,(V)
print*,theta,thetarad, capr, smlr, V

end program

