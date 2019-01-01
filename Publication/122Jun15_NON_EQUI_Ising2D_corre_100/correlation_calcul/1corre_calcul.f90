   PROGRAM max_min
   IMPLICIT NONE

   INTEGER (KIND=8),PARAMETER::n=16
   REAL (KIND=8),PARAMETER:: Jc=0.44

   INTEGER (KIND=8) :: i
   REAL    (KIND=8) :: J,x,y,z,l,eps
   
!!!=======================================================================================
!!!======================================================================================= 
!!! Doc gia tri
   OPEN(unit=11,file='correlation.dat')
   OPEN(unit=12,file='correlation_new.dat')

   DO i=1,n
      READ(11,*)J,x,y,z,l
      eps=1.-J/Jc
      WRITE(12,*)J,eps,abs(1./eps),l
   END DO  
         

   CLOSE(11)
   CLOSE(12)
    
   END PROGRAM max_min

