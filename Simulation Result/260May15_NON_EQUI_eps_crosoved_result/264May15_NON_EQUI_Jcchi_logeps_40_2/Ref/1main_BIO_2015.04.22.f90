!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!  HOANG Danh Tai - Asia Pacific Center for Theoretical Physics, 
!  Hogil Kim Memorial Building #501 POSTECH,
!  San 31, Hyoja-dong, Namgum, Pohang, Gyeongbuk 790-784, Korea.
!  E-mail: hoangdanhtai@gmail.com  
!  https://sites.google.com/site/hoangdanhtai/
!-----------------------------------------------------------------------------------------!
!!! 2015.03.26: Histogram for Magnetization
!!! 2015.04.14: Bo sung tinh Cv, Ksi, Xoa T
!!! 2015.04.21: Sua la cach tinh Ksi: khong tinh theo M, khong theo abs(M)
!!!             Chi tinh phan M>0 cho ordered phase
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 
   PROGRAM main_BIO
   IMPLICIT NONE
 
   CHARACTER (LEN=3)  :: SAt
   CHARACTER (LEN=150):: GS
   CHARACTER (LEN=50) :: name
   CHARACTER (LEN=15) :: tmp
   
   CHARACTER (256)    :: line21,line31  
   REAL    (KIND=8),PARAMETER :: nul=0.,aaa=2.5

   INTEGER (KIND=8) :: i,j,i1,i2,nx,ny,nz,natom,n_av,i_av,i_loop,n_equi,iJ,nJ,number
   INTEGER (KIND=8) :: time_2,time_3,time_5,time_6,time_7,compt,i_histo,nM
   REAL (KIND=8)    :: n_atom,E,E1,E_av,Jmp,rdn_s,rdn_mp,M1_av,M12_av,npositive,Ksi1
   REAL (KIND=8)    :: Jmin,Jmax,delJ,Mmin,Mmax,delM,M,M1,M2,M_av,M2_av,Ksi,E2,E2_av,Cv

   REAL (KIND=8),DIMENSION(:),ALLOCATABLE :: S,M_histo,H_histo,P_histo
   INTEGER (KIND=8),DIMENSION(:),ALLOCATABLE :: nn1,nn2,nn3,nn4,x,y,z
       
!!!=======================================================================================   
!!!=======================================================================================
   CALL system('rm config_ini_3D.pdb')
   CALL system('rm -r config_3D')
   CALL system('mkdir config_3D')
   CALL system('rm -r histogram')
   CALL system('mkdir histogram')
   
   CALL system('rm *.dat*')
   
   OPEN(unit=21,file='E_av.dat')
   
   CALL ini_rdm_number()
   CALL read_input_parameter()
      
   natom=nx*ny*nz
   n_atom=real(natom)

   ALLOCATE(S(natom))
   ALLOCATE(nn1(natom),nn2(natom),nn3(natom),nn4(natom))
   ALLOCATE(x(natom),y(natom),z(natom))
   ALLOCATE(M_histo(nM),H_histo(nM),P_histo(nM))

   IF (nJ==1) THEN
      delJ=0.
   ELSE
      delJ=(Jmax-Jmin)/real(nJ-1)
   END IF

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! ====== MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM ======
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   CALL generate_config()
   CALL index_coordinate()
   CALL nearest_neighbor()
   !CALL write_config_ini_3D()
     
   compt=0
   DO iJ=1,nJ
      Jmp=Jmin+delJ*real(iJ-1)
      compt=compt+1
      WRITE(*,*)compt
         
      CALL equi_lattice1()
      CALL average_thermal()
      ! CALL write_config_3D()
   
   END DO
   
   !CALL value_thermal()
   
   WRITE(*,*)'Finishing'
   CALL computation_time()
   
   CLOSE(21)

   CONTAINS
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE init_rdm_number()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE ini_rdm_number()
   IMPLICIT NONE

   INTEGER (KIND=8),DIMENSION(8) :: time
   INTEGER (KIND=8),DIMENSION(50) :: seed

   CALL DATE_AND_TIME(values=time)     ! Get the current time
   seed(1) = time(4) * (360000*time(5) + 6000*time(6) + 100*time(7) + time(8))
   CALL RANDOM_SEED(PUT=seed)

   time_2=time(2) ; time_3=time(3) ; time_5=time(5) ; time_6=time(6) ; time_7=time(7)

   END SUBROUTINE ini_rdm_number
   
!!!=======================================================================================  
   SUBROUTINE computation_time()
   IMPLICIT NONE

   INTEGER (KIND=8),DIMENSION(8) :: time1
   INTEGER (KIND=8) :: run_date,run_hour,run_minute,run_second
   REAL    (KIND=8) :: run_time

   OPEN (90,file='time_run.dat')

   CALL DATE_AND_TIME(values=time1)     ! Get the current time
      
   run_time = (time1(2)-time_2)*1296000.+(time1(3)-time_3)*86400.&
             +(time1(5)-time_5)*3600.+(time1(6)-time_6)*60.+(time1(7)-time_7) !! second (s)
   run_date=int(run_time/86400.)
   run_hour=int(run_time/3600.-run_date*24.)
   run_minute=int(run_time/60.-run_date*1440.-run_hour*60.)
   run_second=int(run_time-run_date*86400.-run_hour*3600.-run_minute*60.)

   WRITE(90,*)'run_date  :',run_date
   WRITE(90,*)'run_hour  :',run_hour
   WRITE(90,*)'run_minute:',run_minute
   WRITE(90,*)'run_second:',run_second

   CLOSE(90)

   END SUBROUTINE computation_time   

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE read_input_parameter() 
!!! OPEN the parameter from file "parameter.in"
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE read_input_parameter()
   IMPLICIT NONE

   CHARACTER (LEN=150) :: tamp
   OPEN(11,file='1parameter.in')   
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A30,(A10))')   tamp, GS
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A30,(I5))')    tamp, nx
   READ(11, '(A30,(I5))')    tamp, ny
   READ(11, '(A30,(I5))')    tamp, nz
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A30,(F12.6))') tamp, Jmin
   READ(11, '(A30,(F12.6))') tamp, Jmax
   READ(11, '(A30,(I5))')    tamp, nJ
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A30,(I10))')   tamp, n_equi
   READ(11, '(A30,(I10))')   tamp, n_av
   READ(11, '(A50)')         tamp
   READ(11, '(A30,(F12.6))') tamp, Mmin
   READ(11, '(A30,(F12.6))') tamp, Mmax
   READ(11, '(A30,(I5))')    tamp, nM
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp

   CLOSE(11) 

   END SUBROUTINE read_input_parameter

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE index_coordinate() : transfer index i to coordinates x, y, z 
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE index_coordinate()
   IMPLICIT NONE

   DO i=1,natom
      z(i)=(i-1)/(nx*ny)+1
      i1=i-(z(i)-1)*nx*ny
      y(i)=(i1-1)/nx+1
      i2=i1-(y(i)-1)*nx
      x(i)=i2
   END DO

   END SUBROUTINE index_coordinate

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE nearest_neighbor() : find nearest neighbor j of every i
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE nearest_neighbor()
   IMPLICIT NONE

   DO i=1,natom
   DO j=1,natom
   
      IF (((x(j)==(x(i)-1)).or.(x(j)==(x(i)-1+nx))).and.(y(j)==y(i)).and.(z(j)==z(i))) THEN
         nn1(i)=j
      END IF

      IF (((x(j)==(x(i)+1)).or.(x(j)==(x(i)+1-nx))).and.(y(j)==y(i)).and.(z(j)==z(i))) THEN
         nn2(i)=j
      END IF

      IF (((y(j)==(y(i)-1)).or.(y(j)==(y(i)-1+ny))).and.(x(j)==x(i)).and.(z(j)==z(i))) THEN
         nn3(i)=j
      END IF

      IF (((y(j)==(y(i)+1)).or.(y(j)==(y(i)+1-ny))).and.(x(j)==x(i)).and.(z(j)==z(i))) THEN
         nn4(i)=j
      END IF

   END DO
   END DO

   END SUBROUTINE nearest_neighbor
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!  SUBROUTINE value_thermal()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE value_thermal()
   IMPLICIT NONE

   E=0. ; M=0.
   
   DO i=1,natom
      E1=-Jmp*S(i)*(S(nn1(i))+S(nn2(i))+S(nn3(i))+S(nn4(i)))
      
      CALL random_number(rdn_mp)         
      IF (exp(2.*E1) > rdn_mp) THEN
         S(i)=-S(i)
         E1=-E1
         M1=M1+2.*S(i)
      END IF
      E=E+E1
      M=M+S(i)     
   END DO
   
   E=E/n_atom/2. ; E2=E*E
   M2=M*M
        
   END SUBROUTINE value_thermal
   
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!  SUBROUTINE value_thermal()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE M_ini()
   IMPLICIT NONE

   M1=0.
   DO i=1,natom
      M1=M1+S(i)  
   END DO
     
   END SUBROUTINE M_ini
   
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE equi_lattice1()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE equi_lattice1()
   IMPLICIT NONE

   DO i_loop=1,n_equi
      CALL equi_lattice()
      !CALL value_thermal()
   END DO

   END SUBROUTINE equi_lattice1
   
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE equi_lattice()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE equi_lattice()
   IMPLICIT NONE

   DO i=1,natom
      E1=-Jmp*S(i)*(S(nn1(i))+S(nn2(i))+S(nn3(i))+S(nn4(i)))
      
      CALL random_number(rdn_mp)         
      IF (exp(2.*E1) > rdn_mp) THEN
         S(i)=-S(i)
      END IF   
   END DO

   END SUBROUTINE equi_lattice
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE average_thermal()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE average_thermal()
   IMPLICIT NONE   

   E_av=0. ; M_av=0. ; E2_av=0. ; M2_av=0. ; M1_av=0. ; M12_av=0.; npositive=0. 

   !!!------ Histogram ------
   P_histo(:)=0. ; M_histo(:)=0.
   delM = (Mmax-Mmin)/real(nM)
   
   CALL M_ini()

   DO i_av=1,n_av
     CALL value_thermal()
     E_av=E_av+E    
     M_av=M_av+M
     E2_av=E2_av+E2
     M2_av=M2_av+M2
     
     IF (M1>0.) THEN
         npositive=npositive+1.
         M1_av=M1_av+M1
         M12_av=M12_av+M1*M1
     END IF
     
     !!! ---- Histogram ----
      IF ((Mmin <= M1).and.(M1 <= Mmax)) THEN              
         i_histo=int((M1-Mmin)/delM)+1
         P_histo(i_histo)=P_histo(i_histo)+1.            

      END IF     
                       
   END DO

   E_av=E_av/real(n_av)
   E2_av=E2_av/real(n_av)
   M_av=M_av/real(n_av)/n_atom
   M2_av=M2_av/real(n_av)/n_atom/n_atom

   Cv= n_atom*(E2_av-E_av**2.)
   Ksi=n_atom*(M2_av-M_av**2.)
   
   M1_av=M1_av/npositive/n_atom
   M12_av=M12_av/npositive/n_atom/n_atom
    
   Ksi1=n_atom*(M12_av-M1_av**2.)

   WRITE(line21,*)Jmp,E_av,M_av,Cv,Ksi,Ksi1
   WRITE(21,'(a)') trim(line21)

   !!! ---- Histogram ----
   number=10000000+compt
   WRITE(tmp,'(I8)') number
   name='_'//TRIM(tmp)
   OPEN(unit=31,file='histogram/P'//trim(name)//'.dat')
     
   DO i_histo=1,nM
      P_histo(i_histo)=P_histo(i_histo)/real(n_av)
      M_histo(i_histo)=Mmin+(real(i_histo)-0.5)*delM

      WRITE(line31,*)Jmp,M_histo(i_histo)/n_atom,P_histo(i_histo),-log(P_histo(i_histo))
      WRITE(31,'(a)') trim(line31)
   ENDDO     
 
   CLOSE(31)
 
   END SUBROUTINE average_thermal
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!  SUBROUTINE generate_config()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE generate_config()
   IMPLICIT NONE

!!---------------------------------------------------
   IF (GS == 'UP') THEN
         S(:)=1.    
   END IF
!!---------------------------------------------------
   IF (GS == 'NO') THEN
   DO i=1,natom
      CALL random_number(rdn_s)
      IF (rdn_s<0.5) THEN
         S(i)=1.
      ELSE
         S(i)=-1.
      END IF   
   END DO     
   END IF
   
   END SUBROUTINE generate_config
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! WRITE initial position configuration in 3D
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   SUBROUTINE write_config_ini_3D()
   IMPLICIT NONE

   OPEN(unit=11,file='config_ini_3D.pdb')
      
   DO i=1,natom
      IF (int(S(i))==1) THEN 
         SAt='S'
      WRITE(11,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
         SAt,x(i)*aaa,y(i)*aaa,z(i)*aaa,nul
      ELSE
                              
      IF (int(S(i))==-1) THEN 
         SAt='C'
      WRITE(11,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
         SAt,x(i)*aaa,y(i)*aaa,z(i)*aaa,nul
      ELSE
         SAt='Cd'
      WRITE(11,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
         SAt,x(i)*aaa,y(i)*aaa,z(i)*aaa,nul

      END IF
      END IF   
            
   END DO

   CLOSE(11)

   END SUBROUTINE write_config_ini_3D      
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! WRITE initial position configuration in 3D
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE write_config_3D()
   IMPLICIT NONE
    
   OPEN(unit=12,file='config_3D/config_3D'//trim(name)//'.pdb')

   DO i=1,natom
      IF (int(S(i))==1) THEN 
         SAt='S'
      WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
         SAt,x(i)*aaa,y(i)*aaa,z(i)*aaa,nul
      ELSE
                              
      IF (int(S(i))==-1) THEN 
         SAt='C'
      WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
         SAt,x(i)*aaa,y(i)*aaa,z(i)*aaa,nul
      ELSE
         SAt='Cd'
      WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
         SAt,x(i)*aaa,y(i)*aaa,z(i)*aaa,nul

      END IF
      END IF   
            
   END DO

   CLOSE(12)

   END SUBROUTINE write_config_3D  
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   END PROGRAM

