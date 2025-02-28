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
  INTEGER :: idummy

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
