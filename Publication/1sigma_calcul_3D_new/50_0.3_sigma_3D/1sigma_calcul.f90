!!! 21.05.2015: Sua lai theo cach cua Mahn-Soo, khong can phai tinh Hamintonian

   PROGRAM sigma_calcul
   IMPLICIT NONE

   REAL (KIND=8),PARAMETER:: delta=0.3
      
   INTEGER (KIND=8),PARAMETER::nbin=2000,ndata=15
   REAL (KIND=8),PARAMETER:: del_Pc=0.0000001
 
   INTEGER (KIND=8) :: i,k
   REAL    (KIND=8) :: eps1,eps2,sigma1,sigma2,sigma,P1_normalized
   REAL    (KIND=8) :: P1c,delta1,P1_total,P2c,delta2,P2_total
   
   REAL (KIND=8),DIMENSION(nbin) :: S1,S2,P1,P2,logP1,logP2

!!!=======================================================================================
!!!======================================================================================= 

!!! Doc gia tri


   OPEN(unit=13,file='result_sigma.dat')
   OPEN(unit=21,file='P1_cut.dat')
   OPEN(unit=22,file='P2_cut.dat')
 
!!!=======================================================================================
!!!=======================================================================================  
   DO k=1,ndata
 
   !k=15
   
   IF (k==1) THEN
      OPEN(unit=11,file='P_10000015.dat')
      OPEN(unit=12,file='P_10000017.dat')
   END IF
   
   IF (k==2) THEN
      OPEN(unit=11,file='P_10000014.dat')
      OPEN(unit=12,file='P_10000018.dat')
   END IF
   
   IF (k==3) THEN
      OPEN(unit=11,file='P_10000013.dat')
      OPEN(unit=12,file='P_10000019.dat')
   END IF
   
   IF (k==4) THEN
      OPEN(unit=11,file='P_10000012.dat')
      OPEN(unit=12,file='P_10000020.dat')
   END IF
   
   IF (k==5) THEN
      OPEN(unit=11,file='P_10000011.dat')
      OPEN(unit=12,file='P_10000021.dat')
   END IF
   
   IF (k==6) THEN
      OPEN(unit=11,file='P_10000010.dat')
      OPEN(unit=12,file='P_10000022.dat')
   END IF
   
   IF (k==7) THEN
      OPEN(unit=11,file='P_10000009.dat')
      OPEN(unit=12,file='P_10000023.dat')
   END IF
   
   IF (k==8) THEN
      OPEN(unit=11,file='P_10000008.dat')
      OPEN(unit=12,file='P_10000024.dat')
   END IF
   
   IF (k==9) THEN
      OPEN(unit=11,file='P_10000007.dat')
      OPEN(unit=12,file='P_10000025.dat')
   END IF
   
   IF (k==10) THEN
      OPEN(unit=11,file='P_10000006.dat')
      OPEN(unit=12,file='P_10000026.dat')
   END IF
    
   IF (k==11) THEN
      OPEN(unit=11,file='P_10000005.dat')
      OPEN(unit=12,file='P_10000027.dat')
   END IF 
   
   IF (k==12) THEN
      OPEN(unit=11,file='P_10000004.dat')
      OPEN(unit=12,file='P_10000028.dat')
   END IF 
   
   IF (k==13) THEN
      OPEN(unit=11,file='P_10000003.dat')
      OPEN(unit=12,file='P_10000029.dat')
   END IF 
   
   IF (k==14) THEN
      OPEN(unit=11,file='P_10000002.dat')
      OPEN(unit=12,file='P_10000030.dat')
   END IF 
   
   IF (k==15) THEN
      OPEN(unit=11,file='P_10000001.dat')
      OPEN(unit=12,file='P_10000031.dat')
   END IF 
      
 
   P1(:)=0. ; P2(:)=0.

   DO i=1,nbin
      READ(11,*)eps1,S1(i),P1(i),logP1(i)
   END DO  
   
   DO i=1,nbin
      READ(12,*)eps2,S2(i),P2(i),logP2(i)
   END DO 
   WRITE(*,*)'Xong doc gia tri P'

!!!=======================================================================================   
!!! Normalize P1
   P1_normalized=0.
   DO i=1,nbin
      IF (S1(i)>=0.) THEN
         P1_normalized=P1_normalized+P1(i)
      END IF     
   END DO

   DO i=1,nbin
      P1(i)=P1(i)/P1_normalized     
   END DO

!!!=======================================================================================   
!!! Tim pc theo tolerence delta0

   P1c=0. ; delta1=0.
   DO WHILE(delta1<delta)
      P1c=P1c+del_Pc

      delta1=0.
      DO i=1,nbin
         IF (S1(i)>=0.) THEN
         
         IF (P1(i)<P1c) THEN
            delta1=delta1+P1(i)
         END IF
         
         END IF
      END DO
      
   END DO

   WRITE(*,*)P1c,delta1
   
   !!!-----------------------------------------
   P2c=0. ; delta2=0.
   DO WHILE(delta2<delta)
      P2c=P2c+del_Pc

      delta2=0.
      DO i=1,nbin
         IF (P2(i)<P2c) THEN
            delta2=delta2+P2(i)
         END IF
      
      END DO
      
   END DO

   WRITE(*,*)P2c,delta2
      
   WRITE(*,*)'Xong tim Pc '   

!!! Test Pc

   P1_total=0. ; P2_total=0.
   DO i=1,nbin
      
      IF (S1(i)>=0.) THEN
      IF (P1(i)>=P1c) THEN
        ! WRITE(21,*)S1(i),P1(i)
         P1_total=P1_total+P1(i)
      END IF
      END IF
      
      IF (P2(i)>=P2c) THEN
        ! WRITE(22,*)S2(i),P2(i)
         P2_total=P2_total+P2(i)
      END IF
      
   END DO
   
   
   WRITE(*,*)'P1, P2 conlai:',P1_total,P2_total

!!!=======================================================================================
!!!=======================================================================================
!!! Tinh tu so
   sigma1=0.
   DO i=1,nbin
   
      IF ((S1(i)>=0.).and.(P1(i)>P1c)) THEN
    
      IF (P2(i)>=0.) THEN     
      sigma1=sigma1+P2(i)        
               
      END IF
      END IF        
   END DO 
!!!=======================================================================================
!!! --------------------------------------------------------------------------------------
!!! Tinh mau so
   sigma2=0.
   DO i=1,nbin
      IF (P2(i)>P2c) THEN

      sigma2=sigma2+P2(i)
      
      END IF        
   END DO 

!!! Tinh sigma   
   sigma=sigma1/sigma2   
      
   WRITE(13,*)eps1,P1c,sigma1,50.*(eps2**0.63),P2c,sigma2,sigma
   WRITE(*,*)eps1,P1c,sigma1,eps2,P2c,sigma2,sigma    

!!!=======================================================================================
!!!=======================================================================================
   
   END DO
   
   
   CLOSE(11)
   CLOSE(12)
   CLOSE(13)
   
   CLOSE(21)
   CLOSE(22)
   
    
   END PROGRAM sigma_calcul

