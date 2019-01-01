   PROGRAM max_min
   IMPLICIT NONE

   INTEGER (KIND=8),PARAMETER::nH=100,nS=2000
   INTEGER (KIND=8),PARAMETER::N=10000
   REAL (KIND=8),PARAMETER:: Jmp=0.44
   REAL (KIND=8),PARAMETER:: Jc=1./2.28, H0=0.43
   REAL (KIND=8),PARAMETER:: Smin=-1.,Smax=1,z=4.

   INTEGER (KIND=8) :: i,iS
   REAL    (KIND=8) :: T,J,S,delS,JMF
   
   REAL (KIND=8),DIMENSION(nS) :: Smp,P,F,FMF

!!!=======================================================================================
!!!======================================================================================= 

!!!! =====================================================================================
!!! F by mean field

   OPEN(13,file='FMF.dat')

   !JMF=(1./z/Jc)*Jmp
   JMF=1.1/z
   
   WRITE(*,*)'JMF=',JMF

   IF (nS==1) THEN
      delS=0.
   ELSE
      delS=(Smax-Smin)/real(nS)
   END IF

   
   DO iS=1,nS
      S=Smin+delS*real(iS-0.5)

      !FMF=-0.5*z*JMF*S*S + S*(1.+exp(-2.*S))/(1.-exp(-2.*S))+0.5*log(1.-S*S)
      
      !FMF=-0.5*z*JMF*S*S + S*(exp(S)+exp(-S))/(exp(S)-exp(-S))+0.5*log(1.-(S*S))
      
      FMF(iS)=-0.5*z*JMF*S*S + S*0.5*log((1.+S)/(1.-S))+0.5*log(1.-(S*S))
      
      !FMF(iS)=-log(2.)+0.5*(1.-z*JMF)*S*S+(1./12.)*(S**4.)

   WRITE(13,*)S,FMF(iS)  
   END DO


!!!! =====================================================================================
!!! F by numerical method


!!! Doc gia tri
   OPEN(unit=11,file='P044.txt')
   OPEN(unit=12,file='F.dat')

   DO i=1,nS
      READ(11,*)T,J,Smp(i),P(i)

      F(i)=log(P(i))/real(N)           

   END DO  
   
   DO i=1,nS
      !!WRITE(12,*)Smp(i)/N,F(i)*FMF(1001)/F(1001)
      
      WRITE(12,*)Smp(i)/N,F(i)
   END DO 
         

   CLOSE(11)
   CLOSE(12)
    
   END PROGRAM max_min

