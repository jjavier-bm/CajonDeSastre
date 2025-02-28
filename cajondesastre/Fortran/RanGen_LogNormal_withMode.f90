program RanGen
IMPLICIT NONE

INTEGER, PARAMETER :: long=500, cycles=100
REAL*8, PARAMETER :: EV=50.d0, Mode=8.d0

INTEGER :: i, j, ent, total, MeanMode, small
REAL*8 :: output, mean, variance, polydispersity, swap
REAL*8, ALLOCATABLE, DIMENSION(:,:) :: table
REAL*8, ALLOCATABLE, DIMENSION(:) :: diameters
LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASK

allocate(table(0:long,2))
allocate(diameters(cycles))
allocate(MASK(cycles))
MASK=.true.

IF (EV<Mode) STOP 'In a log-normal distribution the expected value should be greater than the mode'
MeanMode=Int(100000*EV+100*Mode) !Just to print
small=0
total=0

DO i=0,long
 table(i,1)=i
 table(i,2)=0.d0
END DO

DO i=1,cycles
 Call lognormal(EV,Mode,output)
 diameters(i)=output
 output=nint(10.d0*output)/10.d0
 IF (diameters(i)<1.d0) small=small+1
 DO j=0,long
  IF (j==nint(output)) THEN
   table(j,2)=table(j,2)+1.d0
   total=total+1
  END IF
 END DO ! j
END DO ! i

!Computation of Mean, Variance and Polydispersity Index

mean=sum(diameters)/cycles
variance=0.d0

DO i=1,cycles
 variance=variance+diameters(i)**2
END DO ! i
 variance=variance/cycles
 variance=variance-mean**2
 polydispersity=variance/mean**2

! PRINT
DO i=0,long
 WRITE (MeanMode,*) table(i,1), table(i,2)
END DO
!Sort

DO i=1,cycles
 j=minloc(diameters(i:100))
 swap=diameters(j)
 diameters(j)=diameters(i)
 diameters(i)=swap

END DO

!PRINT2
DO i=1,cycles
 print*, i, diameters(i), diameters(i)-mean
 WRITE(100,*) i, diameters(i)
 MASK(i)=.false.
END DO

print*, 'Distribution contains', total, ' out of', cycles
print*, 'There are', small, 'small NPs out of', cycles
print*, '    Mean = ', mean
print*, 'Variance = ', variance
print*, '      PI = ', polydispersity


!****************************************************
!****************************************************
!****************************************************
!****************************************************
!****************************************************
contains
!----------------------
SUBROUTINE lognormal(ExpVal,Mode,outnum)
IMPLICIT NONE
!To get a random number from a normal distribution centered on 0 and of an input std. deviation
!Taken from Frenkel and Smit's Understandint molecular simulation
REAL*8, INTENT(IN) :: ExpVal,Mode
REAL*8, INTENT(OUT) :: outnum
REAL*8 :: r, v1, v2,randnum
REAL*8 :: stdev, mean, gaussback !Gaussian

stdev = dsqrt(2.d0*dlog(ExpVal/Mode)/3.d0)
mean = dlog(ExpVal**2*Mode)/3.d0

Call gauss(stdev,mean,gaussback)

outnum = dexp(gaussback)

RETURN
END SUBROUTINE lognormal
!-----------------------

!----------------------
SUBROUTINE gauss(stdev,mean,outnum)
IMPLICIT NONE
!To get a random number from a normal distribution centered on 0 and of an input std. deviation
!Taken from Frenkel and Smit's Understandint molecular simulation
REAL*8, INTENT(IN) :: stdev, mean
REAL*8, INTENT(OUT) :: outnum
REAL*8 :: r, v1, v2,randnum

r=2.d0
DO WHILE (r>1.d0)
 Call Ran01(randnum)  
 v1=2.d0*randnum-1.d0
 Call Ran01(randnum)  
 v2=2.d0*randnum-1.d0
 r=v1**2+v2**2
END DO
outnum=v1*dsqrt(-2.d0*dlog(r)/r)
outnum=mean+stdev*outnum

RETURN
END SUBROUTINE gauss
!-----------------------


!-----------------------
  SUBROUTINE Ran01(variable)
  IMPLICIT NONE
  REAL*8, INTENT(OUT) :: variable
  integer(kind=4), dimension(8) :: val
  integer :: idummy

   call date_and_time(values=val)  ! random numbers are linked to time
   idummy=val(3)+val(6)+val(7)+val(8)
   variable=ran2(idummy)

  END SUBROUTINE Ran01

  ! ********************************************************
  ! *******Function ran2 to choose random numbers***********
  ! ********************************************************

  real function ran2(idum)


    integer  idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    real  AM,EPS,RNMX
    parameter (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, & 
         IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,  &
         IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)

    integer  idum2,j,k,iv(NTAB),iy
    save iv,iy,idum2
    data idum2/123456789/,iv/NTAB*0/,iy/0/

    if(idum.le.0) then
       idum=max(-idum,1)
       idum2=idum

       do j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if(idum.lt.0) idum=idum+IM1
          if(j.le.NTAB) iv(j)=idum
       end do

       iy=iv(1)
    end if

    k=idum/IQ1
    idum=IA1*(idum-k*IQ1)-k*IR1
    if(idum.lt.0) idum=idum+IM1
    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2
    if(idum2.lt.0) idum2=idum2+IM2
    j=1+iy/NDIV
    iy=iv(j)-idum2
    iv(j)=idum
    if(iy.lt.1) iy=iy+IMM1
    ran2=min(AM*iy,RNMX)

    return
  end function ran2

end program RanGen
