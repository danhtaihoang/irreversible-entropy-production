
   PROGRAM chi_calcul
   IMPLICIT NONE
     
   INTEGER (KIND=8),PARAMETER::n=41
 
   INTEGER (KIND=8) :: i
   
   REAL (KIND=8),DIMENSION(n) :: eps,E,Cv,M1,M2,M3,chi1,chi2,chi3,slop_chi

!!!=======================================================================================
!!!======================================================================================= 
!!!=======================================================================================
!!!=======================================================================================  
   OPEN(unit=11,file='E_av_40_1.dat')  
   OPEN(unit=12,file='chi_slop_40_1.txt') 
 
   DO i=1,n
      READ(11,*)eps(i),E(i),Cv(i),M1(i),M2(i),M3(i),chi1(i),chi2(i),chi3(i)
   END DO  
   
   slop_chi(:)=0.
   
   DO i=1,n-1
      WRITE(*,*)chi1(i),chi1(i+1),eps(i),eps(i+1)
      slop_chi(i+1)=(log(chi1(i+1))-log(chi1(i)))/(log(abs(eps(i+1)))-log(abs(eps(i))))
            
      WRITE(12,*)abs(eps(i+1)),slop_chi(i+1)
   END DO 
  
!!!=======================================================================================
      
   CLOSE(11)
   CLOSE(12)  
    
   END PROGRAM chi_calcul

