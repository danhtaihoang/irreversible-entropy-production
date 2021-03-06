!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!  HOANG Danh Tai - Asia Pacific Center for Theoretical Physics, 
!  Hogil Kim Memorial Building #501 POSTECH,
!  San 31, Hyoja-dong, Namgum, Pohang, Gyeongbuk 790-784, Korea.
!  E-mail: hoangdanhtai@gmail.com  
!  https://sites.google.com/site/hoangdanhtai/
!-----------------------------------------------------------------------------------------!
!! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 
   PROGRAM main_BIO
   IMPLICIT NONE
 
   CHARACTER (LEN=3)  :: SAt
   CHARACTER (LEN=150):: GS
   CHARACTER (LEN=50) :: name
   CHARACTER (LEN=15) :: tmp
   
   CHARACTER (256)    :: line21,line31  
   REAL    (KIND=8),PARAMETER :: nul=0.,aaa=2.5

   INTEGER (KIND=8) :: i,j,i1,i2,nx,ny,nz,natom,nJ,n_av,i_J,i_av,i_loop,n_equi,iT,nT,number
   INTEGER (KIND=8) :: time_2,time_3,time_5,time_6,time_7,compt,i_histo,nE
   REAL (KIND=8)    :: n_atom,E,E1,E_av,Jmin,Jmax,J_mp,delJ,rdn_s,rdn_mp
   REAL (KIND=8)    :: Tmin,Tmax,delT,T,Emin,Emax,delE

   REAL (KIND=8),DIMENSION(:),ALLOCATABLE :: S,E_histo,H_histo,P_histo
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
   ALLOCATE(E_histo(nE),H_histo(nE),P_histo(nE))

   IF (nT==1) THEN
      delT=0.
   ELSE
      delT=(Tmax-Tmin)/real(nT-1)
   END IF

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
   CALL write_config_ini_3D()
     
   compt=0
   DO iT=1,nT
      T=Tmin+delT*real(iT-1)
   
      DO i_J=1,nJ
         compt=compt+1
         J_mp=Jmin+delJ*real(i_J-1)
         WRITE(*,*)compt
         
         CALL equi_lattice1()
         CALL average_thermal()
         CALL write_config_3D()
      END DO
   
   END DO
   
   !CALL value_thermal()
   
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
   READ(11, '(A30,(F12.6))') tamp, Tmin
   READ(11, '(A30,(F12.6))') tamp, Tmax
   READ(11, '(A30,(I5))')    tamp, nT
   READ(11, '(A50)')         tamp
   READ(11, '(A30,(F12.6))') tamp, Jmin
   READ(11, '(A30,(F12.6))') tamp, Jmax
   READ(11, '(A30,(I5))')    tamp, nJ
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A30,(I10))')   tamp, n_equi
   READ(11, '(A30,(I10))')   tamp, n_av
   READ(11, '(A50)')         tamp
   READ(11, '(A30,(F12.6))') tamp, Emin
   READ(11, '(A30,(F12.6))') tamp, Emax
   READ(11, '(A30,(I5))')    tamp, nE
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

   E=0.
   DO i=1,natom
      E1=-S(i)*(S(nn1(i))+S(nn2(i))+S(nn3(i))+S(nn4(i)))
      
      CALL random_number(rdn_mp)         
      IF (exp(2.*J_mp*E1/T) > rdn_mp) THEN
         S(i)=-S(i)
         E1=-E1
      END IF
      E=E+E1     
   END DO

   E=E/n_atom/2.
   
   !WRITE(*,*)E
     
   END SUBROUTINE value_thermal
   
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
      E1=-S(i)*(S(nn1(i))+S(nn2(i))+S(nn3(i))+S(nn4(i)))
      
      CALL random_number(rdn_mp)         
      IF (exp(2.*J_mp*E1/T) > rdn_mp) THEN
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

   E_av=0.

   !!!------ Histogram ------
   P_histo(:)=0. ; E_histo(:)=0.
   delE = (Emax-Emin)/real(nE)

   DO i_av=1,n_av
     CALL value_thermal()
     E_av=E_av+E    
     
      !!! ---- Histogram ----
      IF ((Emin <= E).and.(E <= Emax)) THEN              
         i_histo=int((E-Emin)/delE)+1
         P_histo(i_histo)=P_histo(i_histo)+1.            

      END IF     
                       
   END DO

   E_av=E_av/real(n_av)

   WRITE(line21,*)T,J_mp,E_av
   WRITE(21,'(a)') trim(line21)

   !!! ---- Histogram ----
   number=10000000+compt
   WRITE(tmp,'(I8)') number
   name='_'//TRIM(tmp)
   OPEN(unit=31,file='histogram/P'//trim(name)//'.dat')
   
   DO i_histo=1,nE
      P_histo(i_histo)=P_histo(i_histo)/real(n_av)
      E_histo(i_histo)=Emin+(real(i_histo)-0.5)*delE

      WRITE(line31,*)T,J_mp,E_histo(i_histo),P_histo(i_histo)
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

