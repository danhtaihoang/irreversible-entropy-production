!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!     HOANG Danh Tai - Laboratoire de Physique Théorique et Modélisation 
!     UMR 8089 CNRS-Université de Cergy-Pontoise
!     2, Avenue Adolphe Chauvin, 95032 Cergy-Pontoise Cedex, France
!----------------------------------------------------------------------------------------------------!
!     PROGRAMME: MONTE CARLO TRANSPORT OF SPINS
!     !! 16/11/2010: Tinh gia tri TB cua E, M, EM, Cv, Ksi, V theo pp histogram
!     !! 22.5.12: Chuong trinh cho DEFECT POTTS  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!            

      PROGRAM main_histogram
      IMPLICIT NONE

      CHARACTER (256) Ligne21,Ligne31
      
      CHARACTER (LEN=150) :: test_P_at_Ti

      INTEGER (KIND=4) :: natx,naty,n_atom,n_EN_histo,n_Ti_histo,i_histo,i_Ti_histo
      REAL    (KIND=8) :: delT_histo,E_moy

      REAL    (KIND=8) :: To,T,Z_histo,EN_moy_histo,Mz_moy_histo,ENd_moy_histo
      REAL    (KIND=8) :: EN_2_moy_histo,Mz_2_moy_histo,Mz_EN_moy_histo,Mz_2_EN_moy_histo
      REAL    (KIND=8) :: Cv_histo,Ksi_histo,V1_histo,V2_histo

      REAL    (KIND=8),DIMENSION(:),ALLOCATABLE :: EN_histo,Mz_histo,P_histo,Po_histo,ENd_histo
      REAL    (KIND=8),DIMENSION(:),ALLOCATABLE :: EN_2_histo,Mz_2_histo,Mz_EN_histo,Mz_2_EN_histo         

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!!!!!!!!!!!!!!!!!!!!!! CHUONG TRINH CHINH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      CALL read_input_parameter_file()

      ALLOCATE(EN_histo(n_EN_histo))
      ALLOCATE(ENd_histo(n_EN_histo))
      ALLOCATE(Mz_histo(n_EN_histo))
      ALLOCATE(P_histo(n_EN_histo))
      ALLOCATE(Po_histo(n_EN_histo))
      ALLOCATE(EN_2_histo(n_EN_histo))
      ALLOCATE(Mz_2_histo(n_EN_histo))
      ALLOCATE(Mz_EN_histo(n_EN_histo))
      ALLOCATE(Mz_2_EN_histo(n_EN_histo))


      n_atom = natx
      E_moy=-31000.

!====================================
      
      CALL open_value_histogram()
      
      IF (test_P_at_Ti=='YES') THEN            
            CALL test_P_histo_at_Ti()
      ELSE
            CALL average_thermal_histogram()
      END IF

      CONTAINS


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! Doc gia tri parameter tu file parameter.in
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
      SUBROUTINE read_input_parameter_file()
      IMPLICIT NONE

      CHARACTER (LEN=150) :: tamp
      OPEN(11,file='1parameter.in')
      

      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I5))')    tamp, natx
      READ(11, '(A30,(I5))')    tamp, naty
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp,n_EN_histo
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp,n_Ti_histo
      READ(11, '(A30,(F7.4))')  tamp,delT_histo
      READ(11, '(A30,(A10))')   tamp,test_P_at_Ti
      READ(11, '(A30,(F10.4))') tamp,T
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp

      CLOSE(11) 
      END SUBROUTINE read_input_parameter_file


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! Doc gia tri tinh duoc value_histogram tu chuong trinh main_transport_spin
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE open_value_histogram()

      OPEN(unit=12,file='value_P_at_To.dat')

      DO i_histo=1,n_EN_histo
            READ(12,*)To,EN_histo(i_histo),Mz_histo(i_histo),EN_2_histo(i_histo),Mz_2_histo(i_histo),&
                      Mz_EN_histo(i_histo),Mz_2_EN_histo(i_histo),Po_histo(i_histo),ENd_histo(i_histo)
      END DO 

      CLOSE(12)

      END SUBROUTINE open_value_histogram

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! Xem xet do tin cay cua XS tai nhiet do Ti lan can To 
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      SUBROUTINE test_P_histo_at_Ti()

      CALL system('rm value_P_at_Ti.dat')      
      
      OPEN(unit=21,file='value_P_at_Ti.dat') 
             
      Z_histo=0.
      DO i_histo=1,n_EN_histo
            Z_histo=Z_histo+Po_histo(i_histo)*exp((1./To-1./T)*(EN_histo(i_histo)-E_moy))
      END DO


      DO i_histo=1,n_EN_histo

            P_histo(i_histo)=Po_histo(i_histo)*exp((1./To-1./T)*(EN_histo(i_histo)-E_moy))/Z_histo

            WRITE(Ligne21,*)T,EN_histo(i_histo),P_histo(i_histo)
            WRITE(21,'(a)') trim(Ligne21)

      END DO

      CLOSE(21)

      END SUBROUTINE test_P_histo_at_Ti
      
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! Tinh gia tri trung binh tai nhiet do To va cac nhiet do lan can
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE average_thermal_histogram()

      OPEN(unit=31,file='average_thermal_histo.dat')

      DO i_Ti_histo=1,n_Ti_histo
            T=To-real(int(n_Ti_histo/2))*delT_histo+real(i_Ti_histo-1)*delT_histo

      Z_histo=0.     
      EN_moy_histo=0.
      ENd_moy_histo=0.
      Mz_moy_histo=0.
      EN_2_moy_histo=0.
      Mz_2_moy_histo=0.
      Mz_EN_moy_histo=0.
      Mz_2_EN_moy_histo=0.

      !CALL open_value_histogram()

!!!================================================================================================== 
!!! Tinh XS tai diem Ti lan can To

      IF (i_Ti_histo/=(int(n_Ti_histo/2)+1)) THEN

            DO i_histo=1,n_EN_histo
            Z_histo=Z_histo+Po_histo(i_histo)*exp((1./To-1./T)*(EN_histo(i_histo)-E_moy))
            END DO

            DO i_histo=1,n_EN_histo
            P_histo(i_histo)=Po_histo(i_histo)*exp((1./To-1./T)*(EN_histo(i_histo)-E_moy))/Z_histo
            END DO

      ELSE
            DO i_histo=1,n_EN_histo
            P_histo(i_histo)=Po_histo(i_histo)
            END DO      

      END IF

!!!================================================================================================== 
!!! Tinh gia tri trung binh tai diem To va lan can no
      DO i_histo=1,n_EN_histo
            EN_moy_histo=EN_moy_histo+ENd_histo(i_histo)*P_histo(i_histo)
            
           ! ENd_moy_histo=ENd_moy_histo+ENd_histo(i_histo)*P_histo(i_histo)
            Mz_moy_histo=Mz_moy_histo+Mz_histo(i_histo)*P_histo(i_histo)

            EN_2_moy_histo=EN_2_moy_histo+EN_2_histo(i_histo)*P_histo(i_histo)
            Mz_2_moy_histo=Mz_2_moy_histo+Mz_2_histo(i_histo)*P_histo(i_histo)

            Mz_EN_moy_histo=Mz_EN_moy_histo+Mz_EN_histo(i_histo)*P_histo(i_histo)         
            Mz_2_EN_moy_histo=Mz_2_EN_moy_histo+Mz_2_EN_histo(i_histo)*P_histo(i_histo)
      END DO
      
      !n_atom = natx*naty*natz*P_config

      Cv_histo=real(n_atom)*(EN_2_moy_histo-EN_moy_histo**2.)/(T**2.)
      Ksi_histo=real(n_atom)*(Mz_2_moy_histo-Mz_moy_histo**2.)/T
      V1_histo=EN_moy_histo-(Mz_EN_moy_histo/Mz_moy_histo)
      V2_histo=EN_moy_histo-(Mz_2_EN_moy_histo/Mz_2_moy_histo)
      
      WRITE(Ligne31,*)T,EN_moy_histo,Mz_moy_histo,Cv_histo,Ksi_histo,V1_histo,V2_histo,ENd_moy_histo
      WRITE(31,'(a)') trim(Ligne31)
      

      END DO
      CLOSE(31)

      WRITE(*,*)'n_atom=',n_atom
      END SUBROUTINE average_thermal_histogram


!!!--------------------------------------------------------------------------------------------------
      END PROGRAM main_histogram




