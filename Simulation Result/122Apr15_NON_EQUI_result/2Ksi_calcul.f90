!!!=======================================================================================
!!! 2015.04.16: Ksi as a function of epsilon
!!!=======================================================================================
   PROGRAM max_min
   IMPLICIT NONE

   INTEGER (KIND=8),PARAMETER::nJ=81
   REAL (KIND=8),PARAMETER:: Jc=0.4365

   INTEGER (KIND=8) :: i  
   REAL (KIND=8),DIMENSION(nJ) :: J,E,M,Cv,Ksi,E2,M2,eps

!!!=======================================================================================
!!!======================================================================================= 
   OPEN(unit=12,file='E_av.txt')
   OPEN(unit=13,file='Ksi.dat')

   DO i=1,nJ
      READ(12,*)J(i),E(i),M(i),Cv(i),Ksi(i),E2(i),M2(i)
      eps(i)=abs((J(i)/Jc)-1.)
   END DO  
   
   DO i=1,nJ
      eps(i)=abs((J(i)/Jc)-1.)
      
      WRITE(13,*)J(i),eps(i),Ksi(i)
   END DO 
          
   CLOSE(12)
   CLOSE(13)
    
   END PROGRAM max_min

