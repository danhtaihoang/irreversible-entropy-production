   PROGRAM max_min
   IMPLICIT NONE

   INTEGER (KIND=8),PARAMETER::nS=10001
   INTEGER (KIND=8),PARAMETER::N=10000
   REAL (KIND=8),PARAMETER:: Jmp=0.44
   REAL (KIND=8),PARAMETER:: Jc=0.433, H0=0.43
   REAL (KIND=8),PARAMETER:: Smin=-1.,Smax=1,z=4.

   INTEGER (KIND=8) :: i,iS
   REAL    (KIND=8) :: T,J,S,delS,JMF
   
   REAL (KIND=8),DIMENSION(nS) :: Smp,P,F,FMF

!!!=======================================================================================
!!!======================================================================================= 

!!!! =====================================================================================
!!! F by mean field

   OPEN(11,file='FMF.dat')

   JMF=(0.25/Jc)*Jmp
  ! JMF=1./z
   
   WRITE(*,*)'JMF=',JMF

   IF (nS==1) THEN
      delS=0.
   ELSE
      delS=(Smax-Smin)/real(nS)
   END IF
   
   DO iS=1,nS
      S=Smin+delS*real(iS-0.5)
     
      FMF(iS)=-0.5*z*JMF*S*S + S*0.5*log((1.+S)/(1.-S))+0.5*log(1.-(S*S))
      
      !FMF(iS)=-log(2.)+0.5*(1.-z*JMF)*S*S+(1./12.)*(S**4.)

   WRITE(11,*)S,FMF(iS)  
   END DO


!!!! =====================================================================================
!!! F by numerical method


!!! Doc gia tri
   OPEN(unit=12,file='P0440.txt')
   OPEN(unit=13,file='F.dat')

   DO i=1,nS
      READ(12,*)T,J,Smp(i),P(i)
   END DO  
   
   DO i=1,nS
         F(i)=-log(P(i)/P(5001))/real(N) 
         !F(i)=-log(P(i)/P(1001))/real(N)
      !WRITE(12,*)Smp(i)/N,F(i)*FMF(1001)/F(1001)
      
      WRITE(13,*)Smp(i)/N,F(i),F(i)*100.
   END DO 
         

   CLOSE(11)
   CLOSE(12)
   CLOSE(13)
    
   END PROGRAM max_min

