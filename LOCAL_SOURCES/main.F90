PROGRAM mhd_prog
  USE def_type_mesh
  USE initialization
  USE my_util
  USE input_data
  USE symmetric_field
  ! USE arpack_mhd
  USE fourier_to_real_for_vtu
  USE user_data
  USE post_processing_debug
  USE verbose
#include "petsc/finclude/petsc.h"
  USE petsc
  IMPLICIT NONE
  !===Navier-Stokes fields========================================================
  TYPE(mesh_type), POINTER                        :: pp_mesh, vv_mesh
  REAL(KIND=8), POINTER, DIMENSION(:,:,:)         :: un, pn
  TYPE(dyn_real_array_three), POINTER, DIMENSION(:):: der_un
  !===Maxwell fields==============================================================
  TYPE(mesh_type), POINTER                        :: H_mesh, phi_mesh
  TYPE(interface_type), POINTER                   :: interface_H_mu, interface_H_phi
  REAL(KIND=8), POINTER,      DIMENSION(:,:,:)    :: Hn, Bn, phin, vel
  REAL(KIND=8), POINTER,      DIMENSION(:)        :: sigma_field, mu_H_field
  !===Temperature field===========================================================
  TYPE(mesh_type), POINTER                        :: temp_mesh
  REAL(KIND=8), POINTER, DIMENSION(:,:,:)         :: temperature
  REAL(KIND=8), POINTER,      DIMENSION(:)        :: vol_heat_capacity_field
  REAL(KIND=8), POINTER,      DIMENSION(:)        :: temperature_diffusivity_field
  !===Concentration field===========================================================
  TYPE(mesh_type), POINTER                        :: conc_mesh
  REAL(KIND=8), POINTER, DIMENSION(:,:,:)         :: concentration
  REAL(KIND=8), POINTER,      DIMENSION(:)        :: concentration_diffusivity_field
  !===Level_set===================================================================
  REAL(KIND=8), POINTER, DIMENSION(:,:,:,:)       :: level_set
  !===Density=====================================================================
  REAL(KIND=8), POINTER, DIMENSION(:,:,:)         :: density
  !===LES=========================================================================
  REAL(KIND=8), POINTER, DIMENSION(:,:,:,:)       :: visc_LES
  REAL(KIND=8), POINTER, DIMENSION(:,:,:,:)       :: visc_LES_level
  !===Fourier modes===============================================================
  INTEGER                                         :: m_max_c
  INTEGER,      POINTER,      DIMENSION(:)        :: list_mode
  !===Time iterations=============================================================
  REAL(KIND=8)                                    :: time
  INTEGER                                         :: it
  !===Timing======================================================================
  REAL(KIND=8)                                    :: tps, tploc, tploc_max=0.d0
  !===udif===========================
  REAL(KIND=8), POINTER, DIMENSION(:,:,:)         :: udif ! for antisym multiplier
!
    !TEST DCQ
    !===TEST if you want to visualize mu
    INTEGER :: i, bloc_size, m_max_pad, nb_procs, l
    REAL(KIND = 8) :: tor_pol_norms
    !===VTU 2d======================================================================
    CHARACTER(LEN = 3) :: st_mode
    CHARACTER(LEN = 200) :: header
    CHARACTER(LEN = 30) :: name_of_field
    !==TEST to compute the time average field
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: time_average_u, time_average_p, read_DR
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: time_average_H, time_average_B, time_average_phi
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: time_avg_u_square!,int_bar_u  ! we store here the fluctuating kinetic energies
    !REAL(KIND=8), ALLOCATABLE, DIMENSION(:)                  :: ones
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :, :) :: p_dummy
    REAL(KIND = 8) :: dummy_time
    LOGICAL, SAVE :: once = .TRUE., time_average_counter_read = .FALSE.
    INTEGER, SAVE :: time_average_counter = 0
    REAL(KIND = 8), SAVE :: time_average_kenergy = 0.d0, urms = 0.d0
    REAL(KIND = 8), SAVE :: delta_t = 0.d0, avg_delta_t = 0.d0, avg_delta_t_sqr = 0.d0
    REAL(KIND = 8), SAVE :: poloidal = 0.d0, toroidal = 0.d0
    REAL(KIND = 8), SAVE :: rav_poloidal = 0.d0, rav_toroidal = 0.d0
!   REAL(KIND = 8), SAVE :: time_avg_torque = 0.d0
!   REAL(KIND = 8), SAVE :: time_avg_torque_top = 0.d0
!   REAL(KIND = 8), SAVE :: time_avg_torque_bot = 0.d0
!   REAL(KIND = 8) :: torque_t, torque_top_t, torque_bot_t
    !===Mode===========================
    INTEGER                                         :: mode ! not the physical mode
    !TEST DCQ
  !===Declare PETSC===============================================================
  PetscErrorCode :: ierr
  PetscMPIInt    :: rank
  MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d, comm_one_d_ns, comm_one_d_temp
  MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d_conc

  !===Start PETSC and MPI (mandatory)=============================================
  CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  CALL MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)

  !===User reads his/her own data=================================================
  CALL read_user_data('data')

  !===Initialize SFEMANS (mandatory)==============================================
  CALL initial(vv_mesh, pp_mesh, H_mesh, phi_mesh, temp_mesh, conc_mesh,&
       interface_H_phi, interface_H_mu, list_mode, &
       un, pn, Hn, Bn, phin, vel, &
       vol_heat_capacity_field, temperature_diffusivity_field, &
       concentration_diffusivity_field,mu_H_field, sigma_field, time, m_max_c, &
       comm_one_d, comm_one_d_ns, comm_one_d_temp, comm_one_d_conc,temperature, &
       concentration, level_set, density, &
       der_un, visc_LES, visc_LES_level)


!TEST VB 09/12/2024 (HEMERA)
!       der_un, visc_LES_level)
!TEST VB 09/12/2024

!TEST VB 02/12/2024 (RUCHE)
!       der_un, visc_LES, visc_LES_level)
!TEST VB 02/12/2024

!!$  !===Nonlinear restart========================================================
!!$  IF(user%if_nonlin) THEN
!!$   IF(inputs%type_pb == 'nst') THEN
!!$   ALLOCATE(udif(vv_mesh%np,6,m_max_c))
!!$! le 3e argument est le facteur de symetrisation eps_sa=+1 (anti), eps_sa=-1
!!$! (sym)
!!$              CALL champ_total_anti_sym(comm_one_d_ns(1), vv_mesh, list_mode, &
!!$                   +1.d0, un, udif,'u')
!!$              un = un + user%amp_nl * udif
!!$              WRITE(*,*) 'restart avec champ anti ou sym suivant eps_sa'
!!$!             CALL sauvegarde(0,freq_restart)
!!$
!!$   ELSE
!!$     DO mode=1,m_max_c
!!$!           un = user%amp_nl * un
!!$            Hn = user%amp_nl * Hn
!!$     ENDDO
!!$   ENDIF
!!$  ENDIF
  !===============================================================================
  !                        VISUALIZATION WITHOUT COMPUTING                       !
  !===============================================================================
  IF (inputs%if_just_processing) THEN
     inputs%freq_plot=1
     CALL my_post_processing(1)
     CALL error_petsc('End post_processing')
  END IF

  !===============================================================================
  !                        EIGENVALUE PROBLEMS/ARPACK                            !
  !===============================================================================
  !  IF (inputs%if_arpack) THEN
  !     !ATTENTION: m_max_c should be equal to 1, meaning each processors is dealing with 1 Fourier mode
  !     !    CALL solver_arpack_mhd(comm_one_d,H_mesh,phi_mesh,&
  !     !         inputs%dt,list_mode,mu_H_field)
  !     !===Postprocessing to check convergence
  !     IF (inputs%test_de_convergence) THEN
  !        CALL post_proc_test(vv_mesh, pp_mesh, temp_mesh, H_mesh, phi_mesh, list_mode, &
  !             un, pn, Hn, Bn, phin, temperature, level_set, mu_H_field, &
  !             time, m_max_c, comm_one_d, comm_one_d_ns, comm_one_d_temp)
  !        CALL error_Petsc('End of convergence test')
  !        !IF (rank==0) WRITE(*,*) 'End of convergence test'
  !        !RETURN
  !     END IF
  !     !=== Put your postprocessing here
  !
  !     !===End of code for ARPACK problem
  !     CALL error_Petsc('END OF ARPACK, EXITING PRGM')
  !     !IF (rank==0) WRITE(*,*) 'END OF ARPACK, EXITING PRGM'
  !     !RETURN
  !  END IF

  !===============================================================================
  !                        TIME INTEGRATION                                      !
  !===============================================================================

  IF (inputs%if_post_proc_init) THEN
     CALL my_post_processing(0)
  END IF

  !===Start time loop
  tps = user_time()
  DO it = 1, inputs%nb_iteration
     tploc =  user_time()
     time = time + inputs%dt

     CALL run_SFEMaNS(time, it)

     !===My postprocessing
     IF (.NOT.inputs%test_de_convergence) THEN
        CALL my_post_processing(it)
     END IF

     !===Write restart file
     IF (MOD(it, inputs%freq_restart) == 0) THEN
        CALL  save_run(it,inputs%freq_restart)
     ENDIF

        !TEST DCQ
        !===Save p, u average restart file
        IF (MOD(it, user%time_avg_freq_restart) == 0) THEN
            IF (inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd') THEN
                !Get the average sum of the scalar/vector average fields before writing the restart
                time_average_u = time_average_u / time_average_counter
                time_avg_u_square = time_avg_u_square / time_average_counter
                time_average_p = time_average_p / time_average_counter
                CALL my_write_restart_ns(comm_one_d_ns, vv_mesh, pp_mesh, time, &
                        list_mode, time_average_u, time_avg_u_square, time_average_p, p_dummy, &
                        p_dummy, p_dummy, inputs%file_name, it, user%time_avg_freq_restart)
                !Recover the average sum of the scalar/vector average fields
                time_average_u = time_average_u * time_average_counter
                time_avg_u_square = time_avg_u_square * time_average_counter
                time_average_p = time_average_p * time_average_counter
            END IF
            IF (inputs%type_pb=='mxw' .OR. inputs%type_pb=='mhd') THEN
                !Get avg
                time_average_H = time_average_H / time_average_counter
                time_average_B = time_average_B / time_average_counter

                CALL my_write_restart_maxwell(comm_one_d, H_mesh, phi_mesh, time, list_mode, time_average_H, time_average_H, &
                        time_average_B, time_average_B, time_average_phi, time_average_phi, inputs%file_name, &
                        it, user%time_avg_freq_restart)
                !Recover sum
                time_average_H = time_average_H * time_average_counter
                time_average_B = time_average_B * time_average_counter
            END IF
        ENDIF

        !TEST DCQ
     !===Timing
     tploc = user_time() - tploc
     IF (it>1) tploc_max = tploc_max + tploc

  ENDDO

  !===Timing======================================================================
  tps = user_time() - tps
  CALL write_verbose(rank,opt_tps=tps,opt_tploc_max=tploc_max)

  !===Postprocessing to check convergence=========================================
  IF (inputs%if_regression) THEN
     CALL regression(conc_mesh, vv_mesh, pp_mesh, temp_mesh, H_mesh, phi_mesh, list_mode, &
          un, pn, Hn, Bn, phin, temperature, level_set, concentration, mu_H_field, &
          time, m_max_c, comm_one_d, comm_one_d_ns, comm_one_d_temp, comm_one_d_conc)
     CALL error_Petsc('End of convergence test')
  END IF

    !TEST DCQ
    !CALL  my_final_processing()
    !TEST DCQ

  !===End of code=================================================================
  CALL error_Petsc('End of SFEMaNS')
CONTAINS

  SUBROUTINE my_post_processing(it)
    USE sub_plot
    USE chaine_caractere
    USE tn_axi
    USE boundary
    USE sft_parallele
    USE verbose
    USE user_data
    USE sfemans_tools
    USE vtk_viz
    USE subroutine_mass
    USE subroutine_ns_with_u
    IMPLICIT NONE
    INTEGER,                             INTENT(IN) :: it
    REAL(KIND=8)                                    :: err, norm, normr, normt, normz
    REAL(KIND=8)                                    :: tmp, normH, avg_int
    INTEGER                                         :: i, k, it_plot
    CHARACTER(LEN=3)                                :: what
    INTEGER                                         :: rank_S, rank_F
    INTEGER                                         :: rank_ns_S, rank_ns_F, nb_S
    REAL(KIND=8), DIMENSION(vv_mesh%np, 2, SIZE(list_mode)) :: sigma_fluid
    REAL(KIND=8), DIMENSION(vv_mesh%np, 2, SIZE(list_mode)) :: level_1_P2
    REAL(KIND=8), DIMENSION(pp_mesh%np, 2, SIZE(list_mode)) :: level_1_P1
    LOGICAL,      SAVE                         :: once=.TRUE., once_plot=.TRUE.
    INTEGER                                    :: my_petscworld_rank, code!, m, mode
    INTEGER                                    :: bloc_size, m_max_pad, nb_procs!, int_nb
    !REAL(KIND=8), DIMENSION(3)                 :: umax_loc, umax
    !REAL(KIND=8)                               :: interface_z_max, interface_z_max_tot
    !REAL(KIND=8)                               :: interface_z_min, interface_z_min_tot
    REAL(KIND = 8), DIMENSION(vv_mesh%np, 6, m_max_c) :: u_rotated
    REAL(KIND = 8) :: helicity
    REAL(KIND = 8) :: ibar_u_t, ibar_u_r, ibar_u_z
    REAL(KIND = 8), DIMENSION(vv_mesh%np, 6, m_max_c) :: fluct_mke_fd
    REAL(KIND = 8) :: fluct_mke
    !===Energies====================================================================
    REAL(KIND=8), DIMENSION(m_max_c)                :: e_c_u, e_cm_h, temp, conc
    CHARACTER(LEN=3)                                :: truc
    CHARACTER(LEN=100)                              :: en_file
    CHARACTER(LEN=100) , SAVE                       :: e_c_u_file, e_cm_h_file, temp_file, conc_file
    REAL(KIND=8), DIMENSION(m_max_c,6)              :: norme_u, norme_h
    CHARACTER(LEN=100)                              :: norme_u_file, norme_h_file
    REAL(KIND=8),DIMENSION(6), SAVE                 :: type_sym_u, type_sym_h
    REAL(KIND=8), DIMENSION(m_max_c)                :: e_c_u_sym, e_c_u_anti, e_cm_h_sym, e_cm_h_anti
    REAL(KIND=8), DIMENSION(m_max_c)                :: e_c_u_north, e_c_u_south
    REAL(KIND=8), DIMENSION(m_max_c)                :: e_cm_h_north, e_cm_h_south
    REAL(KIND=8), DIMENSION(m_max_c)                :: e_c_u_dif, e_c_u_sym_dif, e_c_u_anti_dif
    REAL(KIND=8), DIMENSION(m_max_c)                :: e_c_ud_north, e_c_ud_south
    !===Kinetic energy===============================================================
    REAL(KIND=8), DIMENSION(vv_mesh%np, 6, SIZE(list_mode)) :: momentum
    !===Anemometers=================================================================
    CHARACTER(LEN=200) , SAVE                       :: anemometre_h_file, anemometre_v_file, anemometre_T_file
    CHARACTER(LEN=200) , SAVE                       :: anemometre_conc_file
    INTEGER , ALLOCATABLE, DIMENSION(:) , SAVE      :: anemo_v, anemo_h, anemo_T, anemo_conc
    INTEGER                             , SAVE      :: nb_anemo_v, nb_anemo_h, nb_anemo_T, nb_anemo_conc
    INTEGER                                         :: i_x, i_y
    !===FLUCTUATING VELOCITY========================================================
!   REAL(KIND = 8), DIMENSION(vv_mesh%np, 6, SIZE(list_mode)) :: un_fluct ! real velocity fluctuation : un_fluct = un - time_average_u
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,m_max_c)    :: udif !un - re_theta=vitesse ds ref wall
    !===Errors======================================================================
    REAL(KIND=8), DIMENSION(SIZE(un,1),SIZE(un,2),SIZE(un,3))   :: u_err
    REAL(KIND=8), DIMENSION(SIZE(pn,1),SIZE(pn,2),SIZE(pn,3))   :: p_err
    REAL(KIND=8), DIMENSION(SIZE(level_set,1),SIZE(level_set,2),SIZE(level_set,3),SIZE(level_set,4)) :: level_set_err
    REAL(KIND=8), DIMENSION(SIZE(temperature,1),SIZE(temperature,2),SIZE(temperature,3))       :: temp_err
    REAL(KIND=8), DIMENSION(SIZE(concentration,1),SIZE(concentration,2),SIZE(concentration,3)) :: conc_err
    REAL(KIND=8), DIMENSION(SIZE(Hn,1),SIZE(Hn,2),SIZE(Hn,3))       :: H_err
    REAL(KIND=8), DIMENSION(SIZE(phin,1),SIZE(phin,2),SIZE(phin,3)) :: phi_err
    !concentration, magn field
    REAL(KIND=8), DIMENSION(4) :: norm_err
    REAL(KIND=8)               :: avg
    INTEGER                    :: int_nb
    !===VTU 2D======================================================================
    CHARACTER(LEN=3)   :: st_mode
    CHARACTER(LEN=200) :: header
    CHARACTER(LEN=3)   :: name_of_field
    !===VTU 3D======================================================================
    REAL(KIND=8), DIMENSION(SIZE(concentration,1),SIZE(concentration,2),SIZE(concentration,3)):: molar_fraction

    !===Summary outputs==============================================================
    !Outputs are written every freq_en iterations (to set in data)
    !Plots 3D are written every freq_plot iterations (to set in data)
    !Plots 2D disable by default. To enable, add the following line to data and set to true
    !===Create 2D plots (true/false)
    !=======Energies
    !For each variable, the energy of the Fourier mode i is saved in the following file
    !e_c_u_i (velocity field) (UNIT=20)
    !e_cm_h_i (magnetic field) (UNIT=50)
    !temp_i (temperature) (UNIT=60)
    !conc_i (concentration) (UNIT=61)
    !=======Navier-Stokes and mass conservation equations
    !fort.11 = Time, L2-norm divergence u, (L2-norm divergence u)/(H1-norm u), weak_div_L2_relative
    !fort.31 = Time, L2-norm divergence u, (L2-norm divergence u)/(H1-norm u)*vv_mesh%global_diameter
    !fort.98 = Time, kinetic energy, poloidal/toroidal energy, ... (see below)
    !fort.99 = Time, t_average_counter,delta(t),avg_delta_t,avg_delta_t_sqr,delta_2^2,avg_int,bar{u_t},bar{u_r},bar{u_z},fluct_mean_kin_energy
    !fort.86 = Time, time_average_counter,rav_polo(un),rav_toro(un),rav_gamma(un),rav_polo(bar{u}),rav_toro(bar{u}),rav_gamma(bar{u})
    !fort.88 = Time, norm L2 of pressure
    !fort.666= Time, torque-total
    !fort.665= Time, torque-total, torque-tot-top, torque-tot-bot
    !fort.667= Time, t_average_counter, time_avg_torque, time_avg_torque_top, time_avg_torque_bot
    !fort.888= Time, Enstrophy_over domain without obstacles, Volume
    !fort.889= Time, Enstrophy_over whole domain
    !fort.998= Time, Dissipation_over domain without obstacles
    !fort.999= Time, Dissipation_over whole domain

    !fort.196 = Time, kinetic energy for problem with variable density
    !fort.197 = Time, mass conservation for problem with variable density
    !======Maxwell equations
    !fort.41 = Time, L2-norm divergence B, (L2-norm divergence B)/(H1-norm B)*H_mesh%global_diameter
    !fort.51 = Time, L2-norm divergence B, (L2-norm divergence B)/(L2-norm B), L2-norm B, L2-norm H, 0.5 H.B, helicity
    !fort.52 = Time, time_average_counter, L2 norm of B_bar, L2 norm of H_bar, B_bar \cdot H_bar
    !fort.78 = Time, 0.5*(L2-norm B)**2
    !======Temperature equation
    !fort.299 = Time, L2-norm temperature
    !======Concentration equation
    !fort.199 = Time, L2-norm concentration
    !======Anenometers (disabled by default)
    !Allow to compute value of variables of interest at given positions in the domain
    !To enable, see read_my_data.F90 and look for anenometer
    !======Errors (disabled by defaut)
    !Allow to compute L2/H1 relative errors of variables of interest (velocity, ...)
    !Errors are written in the file fort.10 every freq_en iterations
    !To enable add the following line to data and set to true
    !===Compute L2/H1 relative errors (true/false)
    !===END Summary outputs==========================================================

    !===File format Energies and Anemometers
    103 FORMAT(1500(e22.9, 2x))
    104 FORMAT(e22.9, i4, 1500(e22.9, 2x))
    106 FORMAT(A, x, i4, 2x, A, x, e12.5, 2x, A, x, e12.5)
    109 FORMAT(1500(A))
    110 FORMAT(e22.9, 2x, i4, 1500(e22.9, 2x))

    !===Check ranks
    IF (vv_mesh%me /=0) THEN
       CALL MPI_Comm_rank(comm_one_d_ns(1), rank_ns_S, ierr)
       CALL MPI_Comm_rank(comm_one_d_ns(2), rank_ns_F, ierr)
    ELSE
       rank_ns_S = -1
       rank_ns_F = -1
    END IF
    CALL MPI_Comm_rank(comm_one_d(1), rank_S, ierr)
    CALL MPI_Comm_rank(comm_one_d(2), rank_F, ierr)
    CALL MPI_Comm_size(comm_one_d(1), nb_S, ierr)

    !===Check numerical stability
    IF (inputs%check_numerical_stability) THEN
       IF (inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd' &
            .OR. inputs%type_pb=='mhs') THEN
          norm = norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode, un)
       ELSE IF (inputs%type_pb=='mxw') THEN
          norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn)
       END IF
       IF (norm>1.d8.OR.isnan(norm)) THEN
          WRITE(*,*) ' Norm L2 of magnetic field (mxw) or velocity (nst/mhd/fhd/mhs) ', norm
          CALL error_petsc('From my_post_processing: numerical unstability')
       END IF
    END IF

    IF (once) THEN
       once=.FALSE.
       CALL MPI_COMM_RANK(PETSC_COMM_WORLD,my_petscworld_rank,code)
       IF (inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd' &
            .OR. inputs%type_pb=='mhs') THEN
          CALL MPI_COMM_SIZE(comm_one_d_ns(2), nb_procs, code)

! udif
          DO i = 1, m_max_c
             DO k = 1, 6
                udif(:,k,i)= un(:,k,i) - vv_exact(k,vv_mesh%rr,list_mode(i),time)

             END DO
          END DO
! udif
                ALLOCATE(time_average_u(vv_mesh%np, 6, SIZE(list_mode)))
                time_average_u = 0.d0
                ALLOCATE(time_avg_u_square(vv_mesh%np, 6, SIZE(list_mode)))
                time_avg_u_square = 0.d0
                ALLOCATE(time_average_p(pp_mesh%np, 2, SIZE(list_mode)))
                time_average_p = 0.d0
                ALLOCATE(p_dummy(pp_mesh%np, 2, SIZE(list_mode)))
                p_dummy = 0.d0
                IF(user%read_time_average_info) THEN
                    !Read u_average & p_average
                    CALL my_read_restart_ns(comm_one_d_ns, vv_mesh, pp_mesh, err, list_mode, time_average_u, time_avg_u_square, &
                            time_average_p, p_dummy, p_dummy, p_dummy, inputs%file_name)

                    !Read data in time_average_info
                    OPEN(UNIT = 12, FILE = 'time_average_info', FORM = 'formatted', STATUS = 'unknown')
                    CALL read_until(12, '===Time average counter')
                    READ(12, *) time_average_counter
                    time_average_counter_read = .TRUE.
                    CALL read_until(12, '===Average Kinetic Energy')
                    READ(12, *) time_average_kenergy
                    time_average_kenergy = time_average_kenergy * time_average_counter
                    CALL read_until(12, '===URMS')
                    READ(12, *)  urms
                    urms = urms * time_average_counter
                    CALL read_until(12, '===Avg Delta')
                    READ(12, *)  avg_delta_t
                    avg_delta_t = avg_delta_t * time_average_counter
                    CALL read_until(12, '===Avg Delta Squared')
                    READ(12, *)  avg_delta_t_sqr
                    avg_delta_t_sqr = avg_delta_t_sqr * time_average_counter
!                    CALL read_until(12, '===Avg torque')
!                    READ(12, *)  time_avg_torque
!                    time_avg_torque = time_avg_torque * time_average_counter
!                    CALL read_until(12, '===Avg torque top')
!                    READ(12, *)  time_avg_torque_top
!                    time_avg_torque_top = time_avg_torque_top * time_average_counter
!                    CALL read_until(12, '===Avg torque bot')
!                    READ(12, *)  time_avg_torque_bot
!                    time_avg_torque_bot = time_avg_torque_bot * time_average_counter

                    !Recover the sum of the scalar/vector average fields
                    time_average_u = time_average_u * time_average_counter
                    time_avg_u_square = time_avg_u_square * time_average_counter
                    time_average_p = time_average_p * time_average_counter
                END IF !time_avg_info

                IF(rank_ns_S == 0) THEN

                    !Write div info
                    WRITE(11, *) '#time, div_L2, div_L2_relative ,weak_div_L2_relative'

                    !===Write Ravelet's Poloidal and Toroidal Energies
                    WRITE(86, 109)  '#time,time_average_counter,rav_polo(un),rav_toro(un),rav_gamma(un),', &
                            'rav_polo(bar{u}),rav_toro(bar{u}),rav_gamma(bar{u})'
                    !===Write Average info in fort.97
                    WRITE(97, 109)   '#t,t_average_counter,kin_e(un),kin_e(bar{u}),delta(t), URMS(t),', &
                            'URMS(t)/Vmax, time_average_kin_energy,*helicity(t)*'
                    !===Kinetic energy, Poloidal/Toroidal, En_r, En_t, En_z
                    WRITE(98, 109) '#t, Kinetic energy, Poloidal/Toroidal, En_r, En_t, En_z'
                    !===Write Average info in fort.99
                    WRITE(99, 109)   '#t,t_average_counter,delta(t),avg_delta_t,avg_delta_t_sqr,delta_2^2,', &
                            'avg_int,bar{u_t},bar{u_r},bar{u_z},fluct_mean_kin_energy'
                    !===Write Average info in fort.667
!                   WRITE(667, 109) '#t, t_average_counter, time_avg_torque, time_avg_torque_top, time_avg_torque_bot'
                END IF

                !===Write e_c_u files
                en_file = 'e_c_u'
                DO i = 1, m_max_c
                    WRITE(truc, '(i3)') list_mode(i)
                    en_file = TRIM(ADJUSTL(en_file)) // '_' // TRIM(ADJUSTL(truc))
                END DO
                e_c_u_file = en_file
                !         DO i=1, m_max_c
                !            e_c_u(i) = 0.5*(norme_L2_champ_par(comm_one_d_ns(1), vv_mesh, list_mode(i:i), un(:,:,i:i)))**2
                !         ENDDO
                ! Central symmetry
                CALL val_ener_sym_centrale(comm_one_d_ns(1), vv_mesh, list_mode, un, &
                     e_c_u, e_c_u_sym, e_c_u_anti, 'u')
                CALL val_ener_north_south(comm_one_d_ns(1), vv_mesh, list_mode, un, e_c_u_north, e_c_u_south)
                CALL val_ener_sym_centrale(comm_one_d_ns(1), vv_mesh, list_mode, udif, e_c_u_dif, &
                     e_c_u_sym_dif, e_c_u_anti_dif, 'u')
                CALL val_ener_north_south(comm_one_d_ns(1), vv_mesh, list_mode, udif, e_c_ud_north, e_c_ud_south)
                ! Central symmetry
                IF (rank_ns_S == 0) THEN
                    OPEN(UNIT = 20, FILE = e_c_u_file, FORM = 'formatted', STATUS = 'unknown')
                    WRITE(20, "(A)")'#energy per mode'
                    WRITE(20, "(A)")'# t  ec mode'
                    !            WRITE(20,103) time, e_c_u
                    WRITE(20,103) time, e_c_u, e_c_u_sym, e_c_u_anti, e_c_u_north, e_c_u_south,&
                         e_c_u_dif, e_c_u_sym_dif, e_c_u_anti_dif, e_c_ud_north, e_c_ud_south
                    CLOSE(20)
                END IF
                !===End Write e_c_u files

                !===Write norm_comp_u
                en_file = 'norme_comp_u'
                DO i = 1, m_max_c
                    WRITE(truc, '(i3)') list_mode(i)
                    en_file = TRIM(ADJUSTL(en_file)) // '_' // TRIM(ADJUSTL(truc))
                END DO
                norme_u_file = en_file
                DO i = 1, m_max_c
                    DO l = 1, 6
                        norme_u(i, l) = norme_L2_champ_par(comm_one_d_ns(1), vv_mesh, list_mode(i : i), un(:, l : l, i : i))
                    ENDDO
                ENDDO
                IF (rank_ns_S == 0) THEN
                    OPEN(UNIT = 25, FILE = norme_u_file, FORM = 'formatted', &
                            STATUS = 'unknown')
                    WRITE(25, *)'#normes des composantes par mode'
                    WRITE(25, *)'# Temps    C1(m=0..m_nax_c) C2 C3 C4 C5 C6  '

                    WRITE(25, 103) time, norme_u(:, 1), norme_u(:, 2), norme_u(:, 3), &
                            norme_u(:, 4), norme_u(:, 5), norme_u(:, 6)
                    CLOSE(25)
                ENDIF
                !===End Write norm_comp_u

                !===Write anemo_v
                IF (inputs%if_anemo_v) THEN
                    ALLOCATE(anemo_v(SIZE(inputs%r_anemo_v) * SIZE(inputs%z_anemo_v)))
                    nb_anemo_v = 0
                    DO i_x = 1, SIZE(inputs%r_anemo_v, 1)
                        DO i_y = 1, SIZE(inputs%z_anemo_v, 1)
                            anemo_v(nb_anemo_v + 1) = find_point(vv_mesh, inputs%r_anemo_v(i_x), inputs%z_anemo_v(i_y))
                            IF (anemo_v(nb_anemo_v + 1) /= 0) THEN
                                nb_anemo_v = nb_anemo_v + 1
                            END IF
                        END DO
                    END DO
                    en_file = 'anemometre_V'
                    WRITE(truc, '(i3)') rank_ns_S
                    en_file = TRIM(ADJUSTL(en_file)) // '_' // TRIM(ADJUSTL(truc))
                    DO i = 1, m_max_c
                        WRITE(truc, '(i3)') list_mode(i)
                        en_file = TRIM(ADJUSTL(en_file)) // '_' // TRIM(ADJUSTL(truc))
                    END DO
                    anemometre_v_file = en_file
                    OPEN(UNIT = 56, FILE = anemometre_v_file, FORM = 'formatted', STATUS = 'unknown')
                    WRITE(56, "(A)")'#Anemometer for the velocity field'
                    WRITE(56, "(A)")'#Points coordinates:    '
                    DO i = 1, nb_anemo_v
                        WRITE(56, 106) '# point ', i, ' : r=', vv_mesh%rr(1, anemo_v(i)), '; z = ', vv_mesh%rr(2, anemo_v(i))
                    ENDDO
                    CLOSE(56)
                END IF  ! (inputs%if_anemo_v)
            END IF ! nst or mhd or fhd or mhs

       IF (inputs%type_pb=='mxw' .OR. inputs%type_pb=='mhd' .OR. &
            inputs%type_pb=='mxx' .OR. inputs%type_pb=='fhd' .OR. inputs%type_pb=='mhs') THEN

                ! SYMETRIE Rpi
                type_sym_h(1) = 1.d0
                type_sym_h(2) = -1.d0
                type_sym_h(3) = -1.d0
                type_sym_h(4) = 1.d0
                type_sym_h(5) = -1.d0
                type_sym_h(6) = 1.d0
                ! SYMETRIE Rpi
                !===Init time avgs for H,B  & phi
                ALLOCATE(time_average_H(H_mesh%np, 6, SIZE(list_mode)))
                ALLOCATE(time_average_B(H_mesh%np, 6, SIZE(list_mode)))
                ALLOCATE(time_average_phi(phi_mesh%np, 2, SIZE(list_mode)))   !This is dummy var, we do not get the average  of phi
                time_average_H = 0.d0
                time_average_B = 0.d0
                time_average_phi = 0.d0
                IF(rank==0) THEN
                    WRITE(51, *)  '#t,div(Bn),div(Bn)/L2(Bn),L2(Bn),L2(Hn),Magnetic energy,helicity(t)'
                    WRITE(52, *)  '#t,t_avg_counter, L2(bar{Bn}),L2(bar{Hn}), 0.5*dot (bar{H},bar{B})'
                END IF

                !===Read time_averages  H,B  & phi if specified
                IF(user%read_time_average_info) THEN
                    !Read H & B
                    CALL  my_read_restart_maxwell(comm_one_d, H_mesh, phi_mesh, time, list_mode, time_average_H, time_average_H, &
                            time_average_B, time_average_B, time_average_phi, time_average_phi, inputs%file_name)
                    IF(time_average_counter_read .EQV. .FALSE.) THEN
                        OPEN(UNIT = 12, FILE = 'time_average_info', FORM = 'formatted', STATUS = 'unknown')
                        CALL read_until(12, '===Time average counter')
                        READ(12, *) time_average_counter
                        time_average_counter_read = .TRUE.
                    ENDIF
                    !Recover the average sum of the scalar/vector average fields
                    time_average_B = time_average_B * time_average_counter
                    time_average_H = time_average_H * time_average_counter
                ENDIF
                !===End Read time_averages  H,B  & phi if specified

                !===Write e_cm_h files
                en_file = 'e_cm_h'
                DO i = 1, m_max_c
                    WRITE(truc, '(i3)') list_mode(i)
                    en_file = TRIM(ADJUSTL(en_file)) // '_' // TRIM(ADJUSTL(truc))
                END DO
                e_cm_h_file = en_file
                !         DO i=1, m_max_c
                !            e_cm_h(i) = 0.5*(norme_L2_champ_par(comm_one_d(1), H_mesh, list_mode(i:i), Hn(:,:,i:i)))**2
                !         ENDDO

                ! Central symmetry
                CALL val_ener_sym_centrale(comm_one_d(1), H_mesh, list_mode, Hn, &
                     e_cm_h, e_cm_h_sym, e_cm_h_anti, 'h')
                CALL val_ener_north_south(comm_one_d(1), H_mesh, list_mode, Hn, e_cm_h_north, e_cm_h_south)
                ! Central symmetry
                IF (rank_S == 0) THEN
                   OPEN(UNIT=50,FILE=e_cm_h_file, FORM='formatted', STATUS='unknown')
                   WRITE(50,"(A)")'#energie par mode dans le conducteur '
                   WRITE(50,"(A)")'# t  em mode'
!                  WRITE(50,103) time, e_cm_h
                   WRITE(50,103) time, e_cm_h, e_cm_h_sym, e_cm_h_anti, e_cm_h_north, e_cm_h_south
                   CLOSE(50)
                END IF
                !===End Write e_cm_h files

                !===Write norme_comp_h
                en_file = 'norme_comp_h'
                DO i = 1, m_max_c
                    WRITE(truc, '(i3)') list_mode(i)
                    en_file = TRIM(ADJUSTL(en_file)) // '_' // TRIM(ADJUSTL(truc))
                END DO
                norme_h_file = en_file
                DO i = 1, m_max_c
                    DO l = 1, 6
                        norme_h(i, l) = norme_L2_champ_par(comm_one_d(1), H_mesh, list_mode(i : i), Hn(:, l : l, i : i))
                    ENDDO
                ENDDO
                IF (rank_S == 0) THEN
                    OPEN(UNIT = 55, FILE = norme_h_file, FORM = 'formatted', STATUS = 'unknown')
                    WRITE(55, *)'#norms of components per mode for H in the conductor'
                    WRITE(55, *)'# Time    C1(m=0..m_nax_c) C2 C3 C4 C5 C6  '
                    WRITE(55, 103) time, norme_h(:, 1), norme_h(:, 2), norme_h(:, 3), &
                            norme_h(:, 4), norme_h(:, 5), norme_h(:, 6)
                    CLOSE(55)
                ENDIF
                !===End Write norme_comp_h

                !===Write anemo_h
                IF (inputs%if_anemo_h) THEN
                    ALLOCATE(anemo_h(SIZE(inputs%r_anemo_h) * SIZE(inputs%z_anemo_h)))
                    nb_anemo_h = 0
                    DO i_x = 1, SIZE(inputs%r_anemo_h, 1)
                        DO i_y = 1, SIZE(inputs%z_anemo_h, 1)
                            anemo_h(nb_anemo_h + 1) = find_point(H_mesh, inputs%r_anemo_h(i_x), inputs%z_anemo_h(i_y))
                            IF (anemo_h(nb_anemo_h + 1) /= 0) THEN
                                nb_anemo_h = nb_anemo_h + 1
                            END IF
                        END DO
                    END DO
                    en_file = 'anemometre_H'
                    WRITE(truc, '(i3)') rank_S
                    en_file = TRIM(ADJUSTL(en_file)) // '_' // TRIM(ADJUSTL(truc))
                    DO i = 1, m_max_c
                        WRITE(truc, '(i3)') list_mode(i)
                        en_file = TRIM(ADJUSTL(en_file)) // '_' // TRIM(ADJUSTL(truc))
                    END DO
                    anemometre_h_file = en_file
                    OPEN(UNIT = 57, FILE = anemometre_h_file, FORM = 'formatted', STATUS = 'unknown')
                    WRITE(57, "(A)") '#Anemometer for the magnetic field'
                    WRITE(57, "(A)") '#points coordinates :    '
                    DO i = 1, nb_anemo_h
                        WRITE(57, 106) '# point ', i, 'r=', H_mesh%rr(1, anemo_h(i)), '; z = ', H_mesh%rr(2, anemo_h(i))
                    ENDDO
                    CLOSE(57)
                END IF ! (inputs%if_anemo_h)
            ENDIF ! mxw or mhd or mxx or fhd or mhs
!       END IF !end of once is after level_set and concentration
        !TEST DCQ
!!$          en_file = 'e_c_u'
!!$          DO i = 1, m_max_c
!!$             WRITE(truc,'(i3)') list_mode(i)
!!$             en_file = trim(adjustl(en_file))//'_'//trim(adjustl(truc))
!!$          END DO
!!$          e_c_u_file = en_file
!!$          DO i=1, m_max_c
!!$             e_c_u(i) = 0.5*(norme_L2_champ_par(comm_one_d_ns(1), vv_mesh, list_mode(i:i), un(:,:,i:i)))**2
!!$          ENDDO
!!$
!!$          IF (rank_ns_S == 0) THEN
!!$             OPEN(UNIT=20,FILE=e_c_u_file, FORM='formatted', STATUS='unknown')
!!$             WRITE(20,"(A)")'#energie par mode'
!!$             WRITE(20,"(A)")'# t  ec mode'
!!$             WRITE(20,103) time, e_c_u
!!$             CLOSE(20)
!!$          END IF
!!$          IF (inputs%if_anemo_v) THEN
!!$             ALLOCATE(anemo_v(SIZE(inputs%r_anemo_v)*SIZE(inputs%z_anemo_v)))
!!$             nb_anemo_v = 0
!!$             DO i_x = 1, SIZE(inputs%r_anemo_v,1)
!!$                DO i_y = 1, SIZE(inputs%z_anemo_v,1)
!!$                   anemo_v(nb_anemo_v+1) = find_point(vv_mesh,inputs%r_anemo_v(i_x),inputs%z_anemo_v(i_y))
!!$                   IF (anemo_v(nb_anemo_v+1) /= 0) THEN
!!$                      nb_anemo_v = nb_anemo_v + 1
!!$                   END IF
!!$                END DO
!!$             END DO
!!$
!!$             en_file = 'anemometre_V'
!!$             WRITE(truc, '(i3)') rank_ns_S
!!$             en_file = trim(adjustl(en_file))//'_'//trim(adjustl(truc))
!!$             DO i = 1, m_max_c
!!$                WRITE(truc,'(i3)') list_mode(i)
!!$                en_file = trim(adjustl(en_file))//'_'//trim(adjustl(truc))
!!$             END DO
!!$             anemometre_v_file = en_file
!!$             OPEN(UNIT=56,FILE=anemometre_v_file, FORM='formatted', STATUS='unknown')
!!$             WRITE(56,"(A)")'#Anemometre pour le champ de vitesse'
!!$             WRITE(56,"(A)")'#Coordonnees des points :    '
!!$             DO i = 1, nb_anemo_v
!!$                WRITE(56,106) '# point ',i, ' : r=', vv_mesh%rr(1,anemo_v(i)), '; z = ', vv_mesh%rr(2,anemo_v(i))
!!$             ENDDO
!!$             CLOSE(56)
!!$          END IF ! (inputs%if_anemo_v)
!!$       ENDIF ! nst or mhd or fhd or mhs
!!$
!!$       IF (inputs%type_pb=='mxw' .OR. inputs%type_pb=='mhd' .OR. &
!!$            inputs%type_pb=='mxx' .OR. inputs%type_pb=='fhd' .OR. inputs%type_pb=='mhs') THEN
!!$          en_file = 'e_cm_h'
!!$          DO i = 1, m_max_c
!!$             WRITE(truc,'(i3)') list_mode(i)
!!$             en_file = trim(adjustl(en_file))//'_'//trim(adjustl(truc))
!!$          END DO
!!$          e_cm_h_file = en_file
!!$          DO i=1, m_max_c
!!$             e_cm_h(i) = 0.5*(norme_L2_champ_par(comm_one_d(1), H_mesh, list_mode(i:i), Hn(:,:,i:i)))**2
!!$          ENDDO
!!$          IF (rank_S == 0) THEN
!!$             OPEN(UNIT=50,FILE=e_cm_h_file, FORM='formatted', STATUS='unknown')
!!$             WRITE(50,"(A)")'#energie par mode dans le conducteur '
!!$             WRITE(50,"(A)")'# t  em mode'
!!$             WRITE(50,103) time, e_cm_h
!!$             CLOSE(50)
!!$          END IF
!!$
!!$          IF (inputs%if_anemo_h) THEN
!!$             ALLOCATE(anemo_h(SIZE(inputs%r_anemo_h)*SIZE(inputs%z_anemo_h)))
!!$             nb_anemo_h = 0
!!$             DO i_x = 1, SIZE(inputs%r_anemo_h,1)
!!$                DO i_y = 1, SIZE(inputs%z_anemo_h,1)
!!$                   anemo_h(nb_anemo_h+1) = find_point(H_mesh,inputs%r_anemo_h(i_x),inputs%z_anemo_h(i_y))
!!$                   IF (anemo_h(nb_anemo_h+1) /= 0) THEN
!!$                      nb_anemo_h = nb_anemo_h + 1
!!$                   END IF
!!$                END DO
!!$             END DO
!!$             en_file = 'anemometre_H'
!!$             WRITE(truc, '(i3)') rank_S
!!$             en_file = trim(adjustl(en_file))//'_'//trim(adjustl(truc))
!!$             DO i = 1, m_max_c
!!$                WRITE(truc,'(i3)') list_mode(i)
!!$                en_file = trim(adjustl(en_file))//'_'//trim(adjustl(truc))
!!$             END DO
!!$             anemometre_h_file = en_file
!!$             OPEN(UNIT=57,FILE=anemometre_h_file, FORM='formatted', STATUS='unknown')
!!$             WRITE(57,"(A)") '#Anemometre pour le champ magnetique'
!!$             WRITE(57,"(A)") '#Coordonnees des points :    '
!!$             DO i = 1, nb_anemo_h
!!$                WRITE(57,106) '# point ',i, 'r=', H_mesh%rr(1,anemo_h(i)), '; z = ', H_mesh%rr(2,anemo_h(i))
!!$             ENDDO
!!$             CLOSE(57)
!!$          END IF ! (inputs%if_anemo_h)
!!$       ENDIF ! mxw or mhd or mxx or fhd or mhs

       IF (inputs%if_temperature) THEN
          en_file = 'temp'
          DO i = 1, m_max_c
             WRITE(truc,'(i3)') list_mode(i)
             en_file = trim(adjustl(en_file))//'_'//trim(adjustl(truc))
          END DO
          temp_file = en_file
          DO i=1, m_max_c
             temp(i) = norm_S(comm_one_d, 'L2', temp_mesh, list_mode(i:i), temperature(:,:,i:i))
          ENDDO
          IF (rank_S == 0) THEN
             OPEN(UNIT=60,FILE=temp_file, FORM='formatted', STATUS='unknown')
             WRITE(60,"(A)")'#temperature par mode'
             WRITE(60,"(A)")'# t  temp mode'
             WRITE(60,103) time, temp
             CLOSE(60)
          END IF
          IF (inputs%if_anemo_T) THEN
             ALLOCATE(anemo_T(SIZE(inputs%r_anemo_T)*SIZE(inputs%z_anemo_T)))
             nb_anemo_T = 0
             DO i_x = 1, SIZE(inputs%r_anemo_T,1)
                DO i_y = 1, SIZE(inputs%z_anemo_T,1)
                   anemo_T(nb_anemo_T+1) = find_point(temp_mesh,inputs%r_anemo_T(i_x),inputs%z_anemo_T(i_y))
                   IF (anemo_T(nb_anemo_T+1) /= 0) THEN
                      nb_anemo_T = nb_anemo_T + 1
                   END IF
                END DO
             END DO
             en_file = 'anemometre_T'
             WRITE(truc, '(i3)') rank_S
             en_file = trim(adjustl(en_file))//'_'//trim(adjustl(truc))
             DO i = 1, m_max_c
                WRITE(truc,'(i3)') list_mode(i)
                en_file = trim(adjustl(en_file))//'_'//trim(adjustl(truc))
             END DO
             anemometre_T_file = en_file
             OPEN(UNIT=58,FILE=anemometre_T_file, FORM='formatted', STATUS='unknown')
             WRITE(58,"(A)") '#Anemometre pour la temperature'
             WRITE(58,"(A)") '#Coordonnees des points :    '
             DO i = 1, nb_anemo_T
                WRITE(58,106) '# point ',i, 'r=', temp_mesh%rr(1,anemo_T(i)), '; z = ', temp_mesh%rr(2,anemo_T(i))
             ENDDO
             CLOSE(58)
          END IF ! (inputs%if_anemo_T)
       ENDIF ! temperature

       IF (inputs%if_concentration) THEN
          en_file = 'conc'
          DO i = 1, m_max_c
             WRITE(truc,'(i3)') list_mode(i)
             en_file = trim(adjustl(en_file))//'_'//trim(adjustl(truc))
          END DO
          conc_file = en_file
          DO i=1, m_max_c
             conc(i) = norm_S(comm_one_d, 'L2', conc_mesh, list_mode(i:i), concentration(:,:,i:i))
          ENDDO
          IF (rank_S == 0) THEN
             OPEN(UNIT=61,FILE=conc_file, FORM='formatted', STATUS='unknown')
             WRITE(61,"(A)")'#concentration par mode'
             WRITE(61,"(A)")'# t  conc mode'
             WRITE(61,103) time, conc
             CLOSE(61)
          END IF
          IF (inputs%if_anemo_conc) THEN
             ALLOCATE(anemo_conc(SIZE(inputs%r_anemo_conc)*SIZE(inputs%z_anemo_conc)))
             nb_anemo_conc = 0
             DO i_x = 1, SIZE(inputs%r_anemo_conc,1)
                DO i_y = 1, SIZE(inputs%z_anemo_conc,1)
                   anemo_conc(nb_anemo_conc+1) = find_point(conc_mesh,inputs%r_anemo_conc(i_x),inputs%z_anemo_conc(i_y))
                   IF (anemo_conc(nb_anemo_conc+1) /= 0) THEN
                      nb_anemo_conc = nb_anemo_conc + 1
                   END IF
                END DO
             END DO
             en_file = 'anemometre_conc'
             WRITE(truc, '(i3)') rank_S
             en_file = trim(adjustl(en_file))//'_'//trim(adjustl(truc))
             DO i = 1, m_max_c
                WRITE(truc,'(i3)') list_mode(i)
                en_file = trim(adjustl(en_file))//'_'//trim(adjustl(truc))
             END DO
             anemometre_conc_file = en_file
             OPEN(UNIT=59,FILE=anemometre_conc_file, FORM='formatted', STATUS='unknown')
             WRITE(59,"(A)") '#Anemometre pour la concentration'
             WRITE(59,"(A)") '#Coordonnees des points :    '
             DO i = 1, nb_anemo_conc
                WRITE(59,106) '# point ',i, 'r=', conc_mesh%rr(1,anemo_conc(i)), '; z = ', conc_mesh%rr(2,anemo_conc(i))
             ENDDO
             CLOSE(59)
          END IF ! (inputs%if_anemo_conc)
       ENDIF ! concentration
    ENDIF ! end once

    !===Put your postprocessing stuff here
    IF (MOD(it,inputs%freq_en) == 0) THEN !Frequency for energies & Frequency for Average Velocity and Pressure
! DCQ
            !===Check divergence of fields
            IF (inputs%verbose_divergence) THEN
                IF (inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd') THEN
                    norm = norm_SF(comm_one_d, 'H1', vv_mesh, list_mode, un)
                    err = norm_SF(comm_one_d, 'div', vv_mesh, list_mode, un)
                    tmp = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un)
                    IF (rank==0) THEN
                        WRITE(*, '(5(A,e10.3))') ' Time = ', time, &
                                ', ||div(un)||_L2/||un||_H1 = ', err / norm, &
                                ', ||div(un)||_L2 = ', err, &
                                ', ||un||_H1 = ', norm, &
                                ', ||un||_L2 = ', tmp
                    END IF
                END IF
                IF (inputs%type_pb/='nst') THEN
                    err = norm_SF(comm_one_d, 'div', H_mesh, list_mode, Bn)
                    norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Bn)
                    normH = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn)
                    IF (rank==0) THEN
                        WRITE(*, '(5(A,e10.3))') ' Time = ', time, &
                                ', ||div(Bn)||_L2/||Bn||_L2 = ', err / norm, &
                                ', ||div(Bn)||_L2  = ', err, &
                                ', ||Bn||_L2 = ', norm, &
                                ', ||Hn||_L2 = ', normH
                    END IF
                END IF
            END IF
! DCQ
       !===Verbose
       CALL write_verbose(rank)

       IF (inputs%if_temperature) THEN
          !CALL tempmax_tempmin(time, comm_one_d_temp, list_mode, temperature)
       ENDIF

       IF (inputs%if_concentration) THEN
          !CALL concmax_concmin(time, comm_one_d_conc, list_mode, concentration)
          !CALL xmax_xmin(time, comm_one_d_conc, list_mode, concentration)
       ENDIF

            !===Keep the counter
            time_average_counter = time_average_counter + 1


       IF (inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd' &
            .OR. inputs%type_pb=='mhs') THEN
          IF (inputs%if_compute_momentum_pseudo_force) THEN
             !===Compute the term -2/pi*integral((1-chi)*(u-u_solid).e_z/dt)
             !===chi is the penalty function, u_solid the velocity of the solid, dt the time step
             !==Output written in the file fort.12
             !CALL FORCES_AND_MOMENTS(time, vv_mesh, comm_one_d_ns, list_mode, un, &
             !       torque_t, torque_top_t, torque_bot_t)

             CALL compute_enstrophy(time, vv_mesh, comm_one_d_ns, list_mode, un)
          END IF

          err = norm_SF(comm_one_d_ns, 'div', vv_mesh, list_mode, un)
          norm = norm_SF(comm_one_d_ns, 'H1', vv_mesh, list_mode, un)
          IF (rank == 0) THEN
             !===Divergence of velocity
             WRITE(31, 103) time, err, (err/norm)*vv_mesh%global_diameter
             WRITE(11, 103) time, err, err / norm, talk_to_me%weak_div_L2
          END IF
          normr = norm_SF(comm_one_d_ns, 'L2', vv_mesh,list_mode,un(:,1:2,:))
          normt = norm_SF(comm_one_d_ns, 'L2', vv_mesh,list_mode,un(:,3:4,:))
          normz = norm_SF(comm_one_d_ns, 'L2', vv_mesh,list_mode,un(:,5:6,:))
          norm = norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode, un)
          IF (rank == 0) THEN
             !===Kinetic energy, Poloidal/Toroidal, En_rt, En_z
             WRITE(98,103) time, 0.5*norm**2, SQRT((normr**2+normz**2)/normt**2), &
                  0.5*normr**2, 0.5*normt**2, 0.5*normz**2
             WRITE(*,*) 'norm L2 of velocity', time, norm
          END IF

         !===Write e_c_u
!!$          DO i=1, m_max_c
!!$             e_c_u(i) = 0.5*(norme_L2_champ_par(comm_one_d_ns(1), vv_mesh, list_mode(i:i), un(:,:,i:i)))**2
!!$          ENDDO
!!$          IF (rank_ns_S == 0) THEN
!!$             OPEN(UNIT=20,FILE=e_c_u_file, FORM='formatted', POSITION = 'append', &
!!$                  STATUS='unknown')
!!$             WRITE(20,103) time, e_c_u

! udif
          DO i = 1, m_max_c
             DO k = 1, 6
                udif(:,k,i)= un(:,k,i) - vv_exact(k,vv_mesh%rr,list_mode(i),time)
             END DO
          END DO
! udif
!
          CALL val_ener_sym_centrale(comm_one_d_ns(1), vv_mesh, list_mode, un, &
               e_c_u, e_c_u_sym, e_c_u_anti, 'u')
          CALL val_ener_north_south(comm_one_d_ns(1), vv_mesh, list_mode, un, e_c_u_north, e_c_u_south)
          CALL val_ener_sym_centrale(comm_one_d_ns(1), vv_mesh, list_mode, udif, e_c_u_dif, &
               e_c_u_sym_dif, e_c_u_anti_dif, 'u')
          CALL val_ener_north_south(comm_one_d_ns(1), vv_mesh, list_mode, udif, e_c_ud_north, e_c_ud_south)
          !
          IF (rank_ns_S == 0) THEN
             OPEN(UNIT=20,FILE=e_c_u_file, FORM='formatted', POSITION = 'append', &
                  STATUS='unknown')
!            WRITE(20,103) time, e_c_u
             WRITE(20,103) time, e_c_u, e_c_u_sym, e_c_u_anti, e_c_u_north, e_c_u_south,&
                  e_c_u_dif, e_c_u_sym_dif, e_c_u_anti_dif, e_c_ud_north, e_c_ud_south

             CLOSE(20)
          ENDIF
         !===End Write e_c_u

         !===Write norm_comp_u
                en_file = 'norme_comp_u'
                DO i = 1, m_max_c
                    WRITE(truc, '(i3)') list_mode(i)
                    en_file = TRIM(ADJUSTL(en_file)) // '_' // TRIM(ADJUSTL(truc))
                END DO
                norme_u_file = en_file
                DO i = 1, m_max_c
                    DO l = 1, 6
                        norme_u(i, l) = norme_L2_champ_par(comm_one_d_ns(1), vv_mesh, list_mode(i : i), un(:, l : l, i : i))
                    ENDDO
                ENDDO
                IF (rank_ns_S == 0) THEN
                    OPEN(UNIT = 25, FILE = norme_u_file, FORM = 'formatted', POSITION = 'append', &
                            STATUS = 'unknown')
                    WRITE(25, 103) time, norme_u(:, 1), norme_u(:, 2), norme_u(:, 3), &
                            norme_u(:, 4), norme_u(:, 5), norme_u(:, 6)
                    CLOSE(25)
                ENDIF
         !===End Write norm_comp_u

         !===Write anemo_v
          IF (inputs%if_anemo_v) THEN
             IF (nb_anemo_v /=0) THEN
                OPEN(UNIT=56,FILE=anemometre_v_file, FORM='formatted', POSITION = 'append', &
                     STATUS='unknown')
                WRITE(56,103) time, un(anemo_v(1:nb_anemo_v),1,:), un(anemo_v(1:nb_anemo_v),2,:), &
                     un(anemo_v(1:nb_anemo_v),3,:), un(anemo_v(1:nb_anemo_v),4,:), &
                     un(anemo_v(1:nb_anemo_v),5,:), un(anemo_v(1:nb_anemo_v),6,:)
                CLOSE(56)
             ENDIF
          ENDIF
         !===End Write anemo_v

         !===Compute Time Average Stuff
                !u time average
                DO i = 1, m_max_c
                    DO k = 1, 6
                        time_average_u(:, k, i) = un(:, k, i) + time_average_u(:, k, i)
                    END DO
                END DO
                !p time average
                DO i = 1, m_max_c
                    DO k = 1, 2
                        time_average_p(:, k, i) = pn(:, k, i) + time_average_p(:, k, i)
                    END DO
                END DO
                !===End Compute Time Average Stuff


!                !TEST LC
!                time_avg_torque = time_avg_torque + torque_t
!                time_avg_torque_top = time_avg_torque_top + torque_top_t
!                time_avg_torque_bot = time_avg_torque_bot + torque_bot_t
!                !TEST LC

                norm = norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode, un)
                time_average_kenergy = time_average_kenergy + 0.5 * norm**2

                !===compute URMS= SQRT(int(u**2dV) / int(dV))
                !urms = urms + SQRT(2.0)*norm/SQRT(user%omega_Vol)
                !TEST CN 24 july 2015
                urms = urms + norm / SQRT(user%omega_Vol)
                !TEST CN 24 july 2015

                time_average_u = time_average_u / time_average_counter  !Get average_u
                tmp = norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode, time_average_u)

                !===Toroidal and poloidal,Ravelet's Definition
                rav_toroidal = norm_S_L1_zero_mode(comm_one_d, vv_mesh, list_mode, time_average_u(:, 3 : 4, :))
                rav_toroidal = rav_toroidal / user%omega_Vol !norm_S_L1_zero_mode multiplies 2*pi already
                rav_poloidal = norm_S_L1_zero_mode(comm_one_d_ns, vv_mesh, list_mode, &
                        SQRT(time_average_u(:, 1 : 2, :)**2 + time_average_u(:, 5 : 6, :)**2))
                rav_poloidal = rav_poloidal / user%omega_Vol

                !===Instantaneous delta: delta(t)
                delta_t = norm**2 / tmp**2
                avg_delta_t = avg_delta_t + delta_t!avg
                avg_delta_t_sqr = avg_delta_t_sqr + delta_t**2!avg

                !===Time average integral
                !int (bar{u} DOmega)/int( DOmega)
                avg_int = mean_int_S(comm_one_d_ns, vv_mesh, list_mode, time_average_u)

                !===Fluctuating toroidal kinetic energy  in all modes -> u_fluct(:,3:4,:)
                !(v_theta - temporal_average(v_theta))**2
                !===Fluctuating poloidal kinetic energy  in all modes -> u_fluct(:,1:2,:)
                !(V_r**2 + V_z**2) - [ (temporal_average(V_r))**2 +(temporal_average(V_z))**2 ]

                DO i = 1, m_max_c
                    !HF 5/12/2019 this ufluct update wrong
                    !u_fluct(:,3:4,i) = u_fluct(:,3:4,i) + un(:,3:4,i)**2  - time_average_u(:,3:4,i)**2
                    !u_fluct(:,1:2,i) = u_fluct(:,1:2,i) + (un(:,1:2,i)**2 + un(:,5:6,i)**2) &
                    !     -(time_average_u(:,1:2,i)**2 +  time_average_u(:,5:6,i)**2)
                    time_avg_u_square(:, :, i) = time_avg_u_square(:, :, i) + un(:, :, i)**2
                    !HF 5/12/2019
                END DO

                !===Integral components of  bar{u}, only  zero mode
                ibar_u_t = user%omega_Vol * mean_int_S(comm_one_d, vv_mesh, list_mode, time_average_u(:, 1 : 2, :))
                ibar_u_r = user%omega_Vol * mean_int_S(comm_one_d, vv_mesh, list_mode, time_average_u(:, 3 : 4, :))
                ibar_u_z = user%omega_Vol * mean_int_S(comm_one_d, vv_mesh, list_mode, time_average_u(:, 5 : 6, :))



                !===Fluctuational mean kinetic energy
                DO i = 1, m_max_c
                    DO k = 1, 6
                        fluct_mke_fd(:, k, i) = un(:, k, i) - time_average_u(:, k, i)
                    ENDDO
                ENDDO
                fluct_mke = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, fluct_mke_fd)
                fluct_mke = fluct_mke**2 / user%omega_Vol

                !===Kinetic helicity
                CALL compute_helicity_Vel(comm_one_d_ns, vv_mesh, list_mode, un, helicity)

                !===Update time_average_u
                time_average_u = time_average_u * time_average_counter  !Recover sum for average_u

                IF (rank == 0) THEN
                    WRITE(97, 110)  time, time_average_counter, &
                            0.5 * norm**2, & !KE of un
                            0.5 * tmp**2, &  !KE of u time average
                            delta_t, & !delta(t)
                            urms / time_average_counter, & !URMS
                            urms / time_average_counter, & !URMS
                            !!urms / (time_average_counter * user%disk_radius * user%solid_vel(:)), & ! URMS/Vmax
                            time_average_kenergy / time_average_counter, & ! average_kin_energy
                            helicity   ! helicity
                    WRITE(99, 110)  time, time_average_counter, &
                            delta_t, & !delta(t)
                            avg_delta_t / time_average_counter, & !bar{delta(t)}
                            avg_delta_t_sqr / time_average_counter, & !bar{delta(t)^2}
                            avg_delta_t_sqr / time_average_counter - avg_delta_t**2 / time_average_counter**2, & !delta_2^2
                            avg_int, & !int (bar{u} DOmega)/int( DOmega)
                            ibar_u_t, ibar_u_r, ibar_u_z, &
                            fluct_mke
!                   WRITE(667, 110) time, time_average_counter, time_avg_torque / time_average_counter, &
!                           time_avg_torque_top / time_average_counter, time_avg_torque_bot / time_average_counter
                END IF



                !===Ravelet's Poloidal and Toroidal instantenuous energies
                poloidal = norm_S_L1_zero_mode(comm_one_d_ns, vv_mesh, list_mode, &
                        SQRT(un(:, 1 : 2, :)**2 + un(:, 5 : 6, :)**2))
                poloidal = poloidal / user%omega_Vol
                toroidal = norm_S_L1_zero_mode(comm_one_d_ns, vv_mesh, list_mode, un(:, 3 : 4, :))
                toroidal = toroidal / user%omega_Vol

                !===Write Ravelet's Poloidal and Toroidal Energies
                IF (rank == 0) THEN
                    WRITE(86, 110)  time, time_average_counter, &
                            poloidal, &  ! Poloidal(un)
                            toroidal, &  ! Toroidal(un)
                            poloidal / toroidal, & !Gamma(un)
                            rav_poloidal, &  ! Poloidal(bar{u})
                            rav_toroidal, &  ! Toroidal(bar{u})
                            rav_poloidal / rav_toroidal !Gamma(bar{u})
                END IF



!!$                !===Complete Poloidal and Toroidal  Energies
!!$                normr = norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode, un(:, 1 : 2, :))
!!$                normt = norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode, un(:, 3 : 4, :))
!!$                normz = norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode, un(:, 5 : 6, :))
!!$                IF (rank == 0) THEN
!!$                    !===Kinetic energy, Poloidal/Toroidal, En_r, En_t, En_z,
!!$                    WRITE(98, 103) time, 0.5 * norm**2, SQRT((normr**2 + normz**2) / normt**2), &
!!$                            0.5 * normr**2, 0.5 * normt**2, 0.5 * normz**2
!!$                END IF

                err = norm_SF(comm_one_d, 'L2', pp_mesh, list_mode, pn)
                IF (rank == 0) THEN
                   WRITE(*,*) 'norm L2 of pressure', time, err
                   WRITE(88,*) time, err
                END IF

       !===DCQ

          IF (inputs%if_level_set) THEN
             !===Compute the term integral(level_set-level_set_t=0)/integral(level_set_t=0)
             !===Output written in file fort.197
             IF (inputs%if_level_set_P2) THEN
                !CALL compute_level_set_conservation(time,vv_mesh,comm_one_d_ns,list_mode,level_set)
             ELSE
                !CALL compute_level_set_conservation(time,pp_mesh,comm_one_d_ns,list_mode,level_set)
             END IF

             !===Compute kinetic energy
             CALL MPI_COMM_SIZE(comm_one_d_ns(2), nb_procs, code)
             bloc_size = SIZE(density,1)/nb_procs+1
             m_max_pad = 3*SIZE(list_mode)*nb_procs/2
             CALL FFT_SCALAR_VECT_DCL(comm_one_d_ns(2), un, density, momentum, 1, nb_procs, &
                  bloc_size, m_max_pad)

             norm=dot_product_SF(comm_one_d_ns, vv_mesh, list_mode, un, momentum)
             IF (rank == 0) THEN
                WRITE(*,*) 'Kinetic energy (integral of momentum times velocity)', time, err
                WRITE(196,*) time, 0.5*err
             END IF
             !===End compute kinetic energy
          END IF ! end if level_set

          !TESTJLG to get the interface height on the frontier
          !IF (inputs%if_level_set) THEN
          !   CALL interface_mpr(rank)
          !ENDIF
          IF (inputs%if_level_set) THEN
             !CALL vmax_height_m0(time,vv_mesh,comm_one_d_ns,list_mode,un,level_set)
          ELSE
             !CALL vmax_height_m0(time,vv_mesh,comm_one_d_ns,list_mode,un)
          ENDIF
          ! TESTJLG to get the interface height on the frontier

       END IF ! end nst or mhd or fhd or mhs CN 02/2024

       IF (inputs%if_temperature) THEN
          err = norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, temperature)
          IF (rank == 0) THEN
             WRITE(299,*) time, err
             WRITE(*,*) 'norm L2 of temperature', time, err
          END IF

          DO i=1, m_max_c
             temp(i) = norm_S(comm_one_d, 'L2', temp_mesh, list_mode(i:i), temperature(:,:,i:i))
          ENDDO
          IF (rank_S == 0) THEN
             OPEN(UNIT=60,FILE=temp_file, FORM='formatted', POSITION = 'append', &
                  STATUS='unknown')
             WRITE(60,103) time, temp
             CLOSE(60)
          ENDIF

          IF (inputs%if_anemo_T) THEN
             IF (nb_anemo_T /=0) THEN
                OPEN(UNIT=58,FILE=anemometre_T_file, FORM='formatted', POSITION = 'append', &
                     STATUS='unknown')
                WRITE(58,103) time, temperature(anemo_T(1:nb_anemo_T),1,:), temperature(anemo_T(1:nb_anemo_T),2,:)
                CLOSE(58)
             ENDIF
          ENDIF
       END IF ! end if_temperature

       IF (inputs%if_concentration) THEN
          err = norm_SF(comm_one_d_conc, 'L2', conc_mesh, list_mode, concentration)
          IF (rank == 0) THEN
             WRITE(199,*) time, err
             WRITE(*,*) 'norm L2 of concentration', time, err
          END IF

          DO i=1, m_max_c
             conc(i) = norm_S(comm_one_d, 'L2', conc_mesh, list_mode(i:i), concentration(:,:,i:i))
          ENDDO
          IF (rank_S == 0) THEN
             OPEN(UNIT=61,FILE=conc_file, FORM='formatted', POSITION = 'append', &
                  STATUS='unknown')
             WRITE(61,103) time, conc
             CLOSE(61)
          ENDIF

          IF (inputs%if_anemo_conc) THEN
             IF (nb_anemo_conc /=0) THEN
                OPEN(UNIT=59,FILE=anemometre_conc_file, FORM='formatted', POSITION = 'append', &
                     STATUS='unknown')
                WRITE(59,103) time, concentration(anemo_conc(1:nb_anemo_conc),1,:), &
                     concentration(anemo_conc(1:nb_anemo_conc),2,:)
                CLOSE(59)
             ENDIF
          ENDIF
          !CALL potential_on_axis(time, H_mesh, conc_mesh, comm_one_d, comm_one_d_conc,list_mode, Hn, concentration)
       END IF ! end if_concentration

       IF (inputs%type_pb=='mxw' .OR. inputs%type_pb=='mhd' &
            .OR. inputs%type_pb=='fhd' .OR. inputs%type_pb=='mhs' ) THEN
          err = norm_SF(comm_one_d, 'div', H_mesh, list_mode, Hn)
          norm = norm_SF(comm_one_d, 'H1', H_mesh, list_mode, Hn)
          IF (rank == 0) THEN
             !===L2 norm of magnetic field
             WRITE(41,*) time, err, (err/norm)*H_mesh%global_diameter
          END IF

!!$       err = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Bn)
!!$       norm = norm_SF(comm_one_d, 'H1', H_mesh, list_mode, Bn)
          !===Get energies and magnetic energy
          err = norm_SF(comm_one_d, 'div', H_mesh, list_mode, Bn)
          norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Bn)
          normr = dot_product_SF(comm_one_d, H_mesh, list_mode, Hn, Bn)
          normH = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn)
          CALL compute_helicity_Mag(comm_one_d, H_mesh, list_mode, Hn, Bn, helicity)

          IF (rank == 0) THEN
              WRITE(51, 103) time, err, err / norm, norm, normH, 0.5 * (normr), helicity
!!$          !===L2 norm of div(Bn)
!!$          WRITE(51,*) time, err, norm
             !===Magnetic energy 0.5*B**2
             WRITE(78,*) time, 0.5*norm**2
             WRITE(*,*) 'norm L2 of magnetic field', time, norm
             WRITE(*,*) 'div(Bn) relative of magnetic field Bn', time, (err/norm)*H_mesh%global_diameter
             WRITE(*,*) '=========================================================='
          END IF

!!$          DO i=1, m_max_c
!!$             e_cm_h(i) = 0.5*(norme_L2_champ_par(comm_one_d(1), H_mesh, list_mode(i:i), Hn(:,:,i:i)))**2
!!$          ENDDO
!!$          IF (rank_S == 0) THEN
!!$             OPEN(UNIT=50,FILE=e_cm_h_file, FORM='formatted',POSITION = 'append',&
!!$                  STATUS='unknown')
!!$             WRITE(50,103) time, e_cm_h
!!$             CLOSE(50)
!!$          END IF
!!$
!!$          IF (inputs%if_anemo_h) THEN
!!$             IF (nb_anemo_h /=0) THEN
!!$                OPEN(UNIT=57,FILE=anemometre_h_file, FORM='formatted', POSITION = 'append', &
!!$                     STATUS='unknown')
!!$                WRITE(57,103) time, Hn(anemo_h(1:nb_anemo_h),1,:), Hn(anemo_h(1:nb_anemo_h),2,:), &
!!$                     Hn(anemo_h(1:nb_anemo_h),3,:), Hn(anemo_h(1:nb_anemo_h),4,:), &
!!$                     Hn(anemo_h(1:nb_anemo_h),5,:), Hn(anemo_h(1:nb_anemo_h),6,:)
!!$                CLOSE(57)
!!$             ENDIF
!!$          ENDIF
!!$       END IF ! end mxw OR mhd OR fhd or mhs

                !===Time Average Stuff
                !H,B time average
                DO i = 1, m_max_c
                    DO k = 1, 6
                        time_average_B(:, k, i) = Bn(:, k, i) + time_average_B(:, k, i)
                        time_average_H(:, k, i) = Hn(:, k, i) + time_average_H(:, k, i)
                    END DO
                END DO
                !Get the avg
                time_average_B = time_average_B / time_average_counter
                time_average_H = time_average_H / time_average_counter
                norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, time_average_B)
                normr = dot_product_SF(comm_one_d, H_mesh, list_mode, time_average_H, time_average_B)
                normH = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, time_average_H)
                IF (rank == 0) THEN
                    WRITE(52, 104)  time, time_average_counter, &
                            norm, &  !L2 norm of B_bar
                            normH, & !L2 norm of H_bar
                            0.5 * normr    !B_bar \cdot H_bar
                END IF

                !Recover  sum
                time_average_B = time_average_B * time_average_counter
                time_average_H = time_average_H * time_average_counter

                !===Write e_cm_h
                !         DO i=1, m_max_c
                !            e_cm_h(i) = 0.5*(norme_L2_champ_par(comm_one_d(1), H_mesh, list_mode(i:i), Hn(:,:,i:i)))**2
                !         ENDDO
                ! Central symmetry
                CALL val_ener_sym_centrale(comm_one_d(1), H_mesh, list_mode, Hn, &
                     e_cm_h, e_cm_h_sym, e_cm_h_anti, 'h')
                CALL val_ener_north_south(comm_one_d(1), H_mesh, list_mode, Hn, e_cm_h_north, e_cm_h_south)
                IF (rank_S == 0) THEN
                   OPEN(UNIT=50,FILE=e_cm_h_file, FORM='formatted',POSITION = 'append',&
                        STATUS='unknown')
!                  WRITE(50,103) time, e_cm_h
                   WRITE(50,103) time, e_cm_h, e_cm_h_sym, e_cm_h_anti, e_cm_h_north, e_cm_h_south
                   CLOSE(50)
                 END IF
                !===End Write e_cm_h

                !===Write norme_comp_h
                en_file = 'norme_comp_h'
                DO i = 1, m_max_c
                    WRITE(truc, '(i3)') list_mode(i)
                    en_file = TRIM(ADJUSTL(en_file)) // '_' // TRIM(ADJUSTL(truc))
                END DO
                norme_h_file = en_file
                DO i = 1, m_max_c
                    DO l = 1, 6
                        norme_h(i, l) = norme_L2_champ_par(comm_one_d(1), H_mesh, list_mode(i : i), Hn(:, l : l, i : i))
                    ENDDO
                ENDDO
                IF (rank_S == 0) THEN
                    OPEN(UNIT = 55, FILE = norme_h_file, FORM = 'formatted', POSITION = 'append', &
                            STATUS = 'unknown')
                    WRITE(55, 103) time, norme_h(:, 1), norme_h(:, 2), norme_h(:, 3), &
                            norme_h(:, 4), norme_h(:, 5), norme_h(:, 6)
                    CLOSE(55)
                ENDIF
                !===End Write norme_comp_h

                !===Write anemo_h
                IF (inputs%if_anemo_h) THEN
                    IF (nb_anemo_h /=0) THEN
                        OPEN(UNIT = 57, FILE = anemometre_h_file, FORM = 'formatted', POSITION = 'append', &
                                STATUS = 'unknown')
                        WRITE(57, 103) time, Hn(anemo_h(1 : nb_anemo_h), 1, :), Hn(anemo_h(1 : nb_anemo_h), 2, :), &
                                Hn(anemo_h(1 : nb_anemo_h), 3, :), Hn(anemo_h(1 : nb_anemo_h), 4, :), &
                                Hn(anemo_h(1 : nb_anemo_h), 5, :), Hn(anemo_h(1 : nb_anemo_h), 6, :)
                        CLOSE(57)
                    ENDIF
                ENDIF
                !===End Write anemo_h
          END IF ! end mxw OR mhd OR fhd or mhs

       IF (inputs%if_compute_error) THEN
          !Error Velocity-pressure
          IF (vv_mesh%np>0) THEN
             DO i = 1, m_max_c
                DO k = 1, 6
                   u_err(:,k,i) = un(:,k,i) - vv_exact(k,vv_mesh%rr,list_mode(i),time)
                END DO
                DO k = 1, 2
                   p_err(:,k,i) = pn(:,k,i) - pp_exact(k,pp_mesh%rr,list_mode(i),time)
                END DO
                IF (list_mode(i) == 0)  THEN
                   CALL Moy(comm_one_d(1),pp_mesh, p_err(:,1,i),avg)
                   p_err(:,1,i) = p_err(:,1,i) - avg
                END IF
             END DO
             norm_err(1) = SQRT(dot_product_SF(comm_one_d_ns,vv_mesh, list_mode, u_err, u_err))
             norm_err(2) = norm_SF(comm_one_d_ns, 'sH1', vv_mesh, list_mode, u_err)
             norm_err(3) = norm_SF(comm_one_d_ns, 'div', vv_mesh, list_mode, un)
             norm_err(4) = norm_SF(comm_one_d_ns, 'L2', pp_mesh, list_mode, p_err)
             IF (rank==0) THEN
                WRITE(10,*) 'Velocity field   #####################'
                WRITE(10,*) 'L2 error on velocity  = ', norm_err(1)
                WRITE(10,*) 'H1 error on velocity  = ', norm_err(2)
                WRITE(10,*) 'L2 norm of divergence = ', norm_err(3)
                WRITE(10,*) 'Pressure field   #####################'
                WRITE(10,*) 'L2 error on pressure  = ', norm_err(4)
             END IF
          END IF

          !Error level set
          IF (inputs%if_level_set) THEN
             IF (inputs%if_level_set_P2) THEN
                DO i = 1, m_max_c
                   DO k = 1, 2
                      DO int_nb = 1, inputs%nb_fluid - 1
                         level_set_err(int_nb,:,k,i) = level_set(int_nb,:,k,i) &
                              -level_set_exact(int_nb,k,vv_mesh%rr,list_mode(i),time)
                      END DO
                   END DO
                END DO
                norm_err(1) = norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode,level_set_err(1,:,:,:))
             ELSE
                DO i = 1, m_max_c
                   DO k = 1, 2
                      DO int_nb = 1, inputs%nb_fluid - 1
                         level_set_err(int_nb,:,k,i) = level_set(int_nb,:,k,i) &
                              -level_set_exact(int_nb,k,pp_mesh%rr,list_mode(i),time)
                      END DO
                   END DO
                END DO
                norm_err(1) = norm_SF(comm_one_d_ns, 'L2', pp_mesh, list_mode,level_set_err(1,:,:,:))
             END IF
             IF (rank==0) THEN
                WRITE(10,*) 'Level set field#####################'
                WRITE(10,*) 'L2 error on level set', norm_err(1)
             END IF
          END IF
          
          !Error temperature
          IF (temp_mesh%np>0) THEN
             DO i = 1, m_max_c
                DO k = 1, 2
                   temp_err(:,k,i) = temperature(:,k,i) - temperature_exact(k,temp_mesh%rr, list_mode(i), time)
                END DO
             END DO
             norm_err(1) = norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, temp_err)
             IF (rank==0) THEN
                WRITE(10,*) 'Temperature field#####################'
                WRITE(10,*) 'L2 error on temperature = ', norm_err(1)
             END IF
          END IF
          
          !Error concentration
          IF (conc_mesh%np>0) THEN
             DO i = 1, m_max_c
                DO k = 1, 2
                   conc_err(:,k,i) = concentration(:,k,i) - concentration_exact(k,conc_mesh%rr, list_mode(i), time)
                END DO
             END DO
             norm_err(1) = norm_SF(comm_one_d_conc, 'L2', conc_mesh, list_mode, conc_err)
             IF (rank==0) THEN
                WRITE(10,*) 'Concentration field#####################'
                WRITE(10,*) 'L2 error on concentration = ', norm_err(1)
             END IF  
          END IF
          
          !Error Magnetic Field
          IF (H_mesh%np>0) THEN
             DO k = 1, 6
                DO i = 1, SIZE(list_mode)
                   H_err(:,k,i) = Hn(:,k,i) - Hexact(H_mesh,k, H_mesh%rr,list_mode(i), mu_H_field, time)
                END DO
             END DO
             norm_err(1) = SQRT(dot_product_SF(comm_one_d, H_mesh, list_mode, H_err, H_err))
             norm_err(2) = norm_SF(comm_one_d, 'sH1', H_mesh, list_mode, H_err)
             norm_err(3) = norm_SF(comm_one_d, 'div', H_mesh, list_mode, Hn)
             IF (rank==0) THEN
                WRITE(10,*) 'Magnetic field   #####################'
                WRITE(10,*) 'L2 error on H  = ', norm_err(1)
                WRITE(10,*) 'H1 error on H  = ', norm_err(2)
                WRITE(10,*) 'L2 norm of divergence = ', norm_err(3)
             END IF
          END IF

          !Error magnetic scalar potential
          IF (phi_mesh%np>0) THEN
             DO i = 1, m_max_c
                DO k = 1, 2
                   phi_err(:,k,i) = phin(:,k,i) - Phiexact(k, phi_mesh%rr, list_mode(i), inputs%mu_phi,time)
                END DO
             END DO
             norm_err(1) = norm_SF(comm_one_d, 'L2', phi_mesh, list_mode, phi_err)
             IF (rank==0) THEN
                WRITE(10,*) 'Concentration field#####################'
                WRITE(10,*) 'L2 error on magnetic scalar potential = ', norm_err(1)
             END IF
          END IF
       END IF ! if_compute_error
    END IF ! end freq_en

    IF (MOD(it,inputs%freq_plot) == 0) THEN
       !===Plot whatever you want here
       IF (once_plot) THEN
          once_plot=.FALSE.
          what = 'new'
       ELSE
          what = 'old'
       END IF
       it_plot = it/inputs%freq_plot

       !===3D/2D Plots for level set and density
       IF (inputs%if_level_set) THEN
          IF (inputs%if_level_set_P2) THEN
             level_1_P2=level_set(1,:,:,:)
             CALL vtu_3d(comm_one_d,level_1_P2, 'vv_mesh', 'Level_1', 'level_1', what, opt_it=it_plot)
             IF (inputs%nb_fluid.GE.3) THEN
                level_1_P2=level_set(2,:,:,:)
                CALL vtu_3d(comm_one_d,level_1_P2, 'vv_mesh', 'Level_2', 'level_2', what, opt_it=it_plot)
             END IF
          ELSE
             level_1_P1=level_set(1,:,:,:)
             CALL vtu_3d(comm_one_d,level_1_P1, 'pp_mesh', 'Level_1', 'level_1', what, opt_it=it_plot)
             IF (inputs%nb_fluid.GE.3) THEN
                level_1_P1=level_set(2,:,:,:)
                CALL vtu_3d(comm_one_d,level_1_P1, 'pp_mesh', 'Level_2', 'level_2', what, opt_it=it_plot)
             END IF
          END IF
          CALL vtu_3d(comm_one_d,density, 'vv_mesh', 'Density', 'density', what, opt_it=it_plot)
          IF (inputs%if_plot_2D) THEN
             !===Proceed as follows to make 2D plots in the Fourier space (using Hn
             !for instance)
             DO i = 1, m_max_c
                WRITE(st_mode,'(I3)') list_mode(i)
                header = 'Ln_'//'mode_'//trim(adjustl(st_mode))
                name_of_field = 'Ln'
                IF (inputs%if_level_set_P2) THEN
                   CALL make_vtu_file_2D(comm_one_d_ns(1), vv_mesh, header, level_1_P2(:,:,i), name_of_field, &
                        what, opt_it=it_plot)
                ELSE
                   CALL make_vtu_file_2D(comm_one_d_ns(1), pp_mesh, header, level_1_P1(:,:,i), name_of_field, &
                        what, opt_it=it_plot)
                END IF
                header = 'Dn_'//'mode_'//trim(adjustl(st_mode))
                name_of_field = 'Dn'
                CALL make_vtu_file_2D(comm_one_d_ns(1), vv_mesh, header, density(:,:,i), name_of_field, &
                     what, opt_it=it_plot)
             END DO
          END IF
          IF (inputs%variation_sigma_fluid) THEN
             CALL reconstruct_variable(comm_one_d_ns, list_mode, pp_mesh, vv_mesh, level_set, &
                  inputs%sigma_fluid, sigma_fluid)
             CALL vtu_3d(comm_one_d,sigma_fluid, 'vv_mesh', 'Sigma', 'sigma', what, opt_it=it_plot)
          END IF
       END IF

       !===3D/2D Plots for velocity and pressure
       IF (inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd' .OR. inputs%type_pb=='mhs') THEN
          CALL vtu_3d(comm_one_d, un, 'vv_mesh', 'Velocity', 'vel', what, opt_it=it_plot)
!         CALL vtu_3d(comm_one_d, pn, 'pp_mesh', 'Pressure', 'pre', what, opt_it=it_plot)
          IF (inputs%if_plot_2D) THEN
             !===2D plots for each mode of the velocity
             DO i = 1, m_max_c
                WRITE(st_mode,'(I3)') list_mode(i)  !=== (CHARACTER(LEN=3)   :: st_mode)
                header = 'Vn_'//'mode_'//trim(adjustl(st_mode)) !=== (CHARACTER(LEN=200) :: header)
                name_of_field = 'Vn' !===(for instance) (CHARACTER(LEN=3)   :: name_of_field)
                CALL make_vtu_file_2D(comm_one_d_ns(1), vv_mesh, header, un(:,:,i), name_of_field, what, opt_it=it_plot)
             END DO
          END IF ! if_plot_2D

                !===Averages
                time_average_u = time_average_u / time_average_counter
                time_average_p = time_average_p / time_average_counter
                time_avg_u_square = time_avg_u_square / time_average_counter
                !Plots
                ! TEST HF-CN 140119
                !         CALL vtu_3d(comm_one_d, time_average_u, 'vv_mesh', 't_avg_u', 't_avg_u',what, opt_it=it_plot)
                !         CALL vtu_3d(comm_one_d, time_average_p, 'pp_mesh', 't_avg_p', 't_avg_p',what, opt_it=it_plot)

                !===Integral components of  bar{u}, only  zero mode
                !ibar_u_t = user%omega_Vol*mean_int_S(comm_one_d, vv_mesh, list_mode,  time_average_u(:,1:2,:))
                !ibar_u_r = user%omega_Vol*mean_int_S(comm_one_d, vv_mesh, list_mode,  time_average_u(:,3:4,:))
                !ibar_u_z = user%omega_Vol*mean_int_S(comm_one_d, vv_mesh, list_mode,  time_average_u(:,5:6,:))
                !int_bar_u(:,1,1)=ones(:)*ibar_u_t
                !int_bar_u(:,3,1)=ones(:)*ibar_u_r
                !int_bar_u(:,5,1)=ones(:)*ibar_u_z
                !===Plot mode 0
                DO i = 1, 1
                    !===Proceed as follows to make 2D plots in the Fourier space (using un for instance)
                    IF(list_mode(i) == 0) THEN
                        WRITE(st_mode, '(I3)') list_mode(i)
                        header = 't_avg_u_' // 'mode_' // TRIM(ADJUSTL(st_mode))
                        name_of_field = 't_avg_u'
                        CALL make_vtu_file_2D(comm_one_d(1), vv_mesh, header, time_average_u(:, :, i), &
                                       name_of_field, what, opt_it = it_plot)
                        header = 'time_avg_u_square_' // 'mode_' // TRIM(ADJUSTL(st_mode))
                        name_of_field = 't_avg_u2'
                        CALL make_vtu_file_2D(comm_one_d(1), vv_mesh, header, time_avg_u_square(:, :, i), &
                                       name_of_field, what, opt_it = it_plot)
                        !header = 'int_bar_u_'//'mode_'//TRIM(ADJUSTL(st_mode))
                        !name_of_field = 'int_bar_u'
                        !CALL make_vtu_file_2D(comm_one_d(1), vv_mesh, header, int_bar_u(:,:,1), name_of_field,  what, opt_it=it_plot)
                    ENDIF
                END DO
                !===Now compute the field having as a system frame the bottom propeller
                DO i = 1, m_max_c
                    DO k = 1, 6
                        u_rotated(:, k, i) = time_average_u(:, k, i) &
                                     - vv_exact(k, vv_mesh%rr, list_mode(i), time)
                    END DO
                END DO
                !Plot
                ! TEST HF-CN 140119
                    CALL vtu_3d(comm_one_d, u_rotated, 'vv_mesh', 't_avg_u_RFB', 't_avg_u_RFB',what, opt_it=it_plot)
                IF (inputs%if_plot_2D) THEN
                   !===2D plots for each mode of the velocity
                   DO i = 1, m_max_c
                      WRITE(st_mode,'(I3)') list_mode(i)  !=== (CHARACTER(LEN=3)   :: st_mode)
                      header = 'Vdif_avg_'//'mode_'//trim(adjustl(st_mode)) !=== (CHARACTER(LEN=200) :: header)
                      name_of_field = 'Vda' !===(for instance) (CHARACTER(LEN=3)   :: name_of_field)
                      CALL make_vtu_file_2D(comm_one_d_ns(1), vv_mesh, header, u_rotated(:,:,i), &
                      name_of_field, what, opt_it=it_plot)
                   END DO
                END IF ! if_plot_2D
                !Recover sum in u,time_avg_u_square and p
                time_average_u = time_average_u * time_average_counter
                time_average_p = time_average_p * time_average_counter
                time_avg_u_square = time_avg_u_square * time_average_counter

       END IF ! 'nst' .OR. 'mhd' .OR. 'fhd' .OR. 'mhs'

       !===3D/2D Plots for temperature
       IF (inputs%if_temperature) THEN
          CALL vtu_3d(comm_one_d,temperature, 'temp_mesh', 'Temperature', 'temp', what, opt_it=it_plot)
          IF (inputs%if_plot_2D) THEN
             !===2D plots for each mode of the temperature
             DO i = 1, m_max_c
                WRITE(st_mode,'(I3)') list_mode(i)
                header = 'Tn_'//'mode_'//trim(adjustl(st_mode))
                name_of_field = 'Tn'
                CALL make_vtu_file_2D(comm_one_d_temp(1), temp_mesh, header, temperature(:,:,i), name_of_field, &
                     what, opt_it=it_plot)
             END DO
          END IF
       END IF

       !===3D/2D Plots for concentration
       IF (inputs%if_concentration) THEN
          CALL vtu_3d(comm_one_d,concentration, 'conc_mesh', 'Concentration', 'conc', what, opt_it=it_plot)
          IF (inputs%if_plot_2D) THEN
             !===2D plots for each mode of the concentration
             DO i = 1, m_max_c
                WRITE(st_mode,'(I3)') list_mode(i)
                header = 'Cn_'//'mode_'//trim(adjustl(st_mode))
                name_of_field = 'Cn'
                CALL make_vtu_file_2D(comm_one_d_conc(1), conc_mesh, header, concentration(:,:,i), name_of_field, &
                     what, opt_it=it_plot)
             END DO
          END IF

          !TEST LC check if we keep or if we add a if coupling_Hx
          !===molar fraction
          CALL MPI_COMM_SIZE(comm_one_d_conc(2), nb_procs, code)
          bloc_size = SIZE(concentration,1)/nb_procs+1
          m_max_pad = 3*SIZE(list_mode)*nb_procs/2
          DO k=1,2
             DO i=1,SIZE(list_mode)
                molar_fraction(:,k,i) = concentration(:,k,i)
             END DO
          END DO
          CALL FFT_PAR_SCAL_FUNCT(comm_one_d_conc(2), molar_fraction,&
               molar_fraction_from_concentration,nb_procs, bloc_size, m_max_pad)
          CALL vtu_3d(comm_one_d,molar_fraction, 'conc_mesh', 'Molar_Fraction', 'mol', what, opt_it=it_plot)
          IF (inputs%if_plot_2D) THEN
             !===2D plots for each mode of the molar_fraction
             DO i = 1, m_max_c
                WRITE(st_mode,'(I3)') list_mode(i)
                header = 'Xn_'//'mode_'//trim(adjustl(st_mode))
                name_of_field = 'Xn'
                CALL make_vtu_file_2D(comm_one_d_conc(1), conc_mesh, header, molar_fraction(:,:,i), name_of_field,&
                     what, opt_it=it_plot)
             END DO
          END IF
          !TEST LC check if we keep or if we add a if coupling_Hx
       END IF

       !===3D/2D Plots for magnetic field
       IF (inputs%type_pb=='mxw' .OR. inputs%type_pb=='mxx' .OR. inputs%type_pb=='mhd' &
            .OR. inputs%type_pb=='fhd' .OR. inputs%type_pb=='mhs') THEN
          CALL vtu_3d(comm_one_d,Hn, 'H_mesh', 'MagField', 'mag', what, opt_it=it_plot)
          CALL vtu_3d(comm_one_d,Hn, 'H_mesh', 'Current', 'cur', what, opt_it=it_plot, opt_grad_curl='curl_h',opt_2D=.TRUE.,&
               opt_mesh_in=H_mesh)
          IF (inputs%nb_dom_phi>0) THEN
             CALL vtu_3d(comm_one_d, phin, 'phi_mesh', 'ScalPot', 'phi', what, opt_it=it_plot)
          END IF

          IF (inputs%if_plot_2D) THEN
             !===2D plots for each mode of the magnetic field
             DO i = 1, m_max_c
                WRITE(st_mode,'(I3)') list_mode(i)
                header = 'Hn_'//'mode_'//trim(adjustl(st_mode))
                name_of_field = 'Hn'
                CALL make_vtu_file_2D(comm_one_d(1), H_mesh, header, Hn(:,:,i), name_of_field, &
                     what, opt_it=it_plot)
             END DO
          END IF ! if_plot_2D

                !Get Averages
                time_average_H = time_average_H / time_average_counter
                time_average_B = time_average_B / time_average_counter
                CALL vtu_3d(comm_one_d, time_average_H, 'H_mesh', 't_avg_mag', 't_avg_mag', what, opt_it = it_plot)
                CALL vtu_3d(comm_one_d, time_average_B, 'H_mesh', 't_avg_mag_in', 't_avg_mgi', what, opt_it = it_plot)

                !Plot mode 0
                DO i = 1, 1
                    !===Proceed as follows to make 2D plots in the Fourier space (using un for instance)
                    IF(list_mode(i) == 0) THEN
                        WRITE(st_mode, '(I3)') list_mode(i)
                        header = 't_avg_h_' // 'mode_' // TRIM(ADJUSTL(st_mode))
                        name_of_field = 't_avg_h'
                        CALL make_vtu_file_2D(comm_one_d(1), H_mesh, header, time_average_H(:, :, i), &
                                       name_of_field, what, opt_it = it_plot)
                        header = 't_avg_b_' // 'mode_' // TRIM(ADJUSTL(st_mode))
                        name_of_field = 't_avg_b'
                        CALL make_vtu_file_2D(comm_one_d(1), H_mesh, header, time_average_B(:, :, i), &
                                       name_of_field, what, opt_it = it_plot)
                    ENDIF
                END DO

                !Recover sum
                time_average_H = time_average_H * time_average_counter
                time_average_B = time_average_B * time_average_counter



                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!$          !===Plot  blades
                !!$          CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs, ierr)
                !!$          m_max_pad = 3*SIZE(list_mode)*nb_procs/2
                !!$          bloc_size = SIZE(mu_field_dummy,1)/nb_procs+1
                !!$
                !!$          mu_field_dummy = 0.d0
                !!$          DO i = 1, SIZE(list_mode)
                !!$             IF (list_mode(i)==0) mu_field_dummy(:,5:6,i) = 1.d0
                !!$          END DO
                !!$          CALL FFT_PAR_VAR_ETA_PROD_T_DCL(comm_one_d(2), mu_in_real_space, &
                !!$               vv_mesh, &
                !!$               mu_field_dummy, mu_field_dummy, nb_procs, bloc_size, m_max_pad,time)
                !!$          !To visualize mu plt
                !!$          CALL vtu_3d(comm_one_d, mu_field_dummy(:,5:6,:), 'vv_mesh', 'mu', 'mu', what, opt_it=it_plot)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       END IF ! 'mxw' .OR. 'mxx' .OR. 'mhd' & .OR. 'fhd' .OR. 'mhs'
    ENDIF ! end freq_plot

  END SUBROUTINE my_post_processing

    SUBROUTINE my_final_processing()
        USE user_data
        USE sub_plot
        USE chaine_caractere
        USE tn_axi
        USE boundary
        IMPLICIT NONE
        REAL(Kind = 8) :: delta, energy
        INTEGER :: i, k
        INTEGER :: bloc_size, m_max_pad, nb_procs
        REAL(KIND = 8), DIMENSION(vv_mesh%np, 6, m_max_c) :: mu_field_dummy

        109 FORMAT(1500(A))
        110 FORMAT(e22.9, 2x, I5, 1500(e22.9, 2x))

        IF (inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd') THEN
            time_average_u = time_average_u / time_average_counter
            delta = time_average_kenergy / time_average_counter
            energy = 0.5 * norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, time_average_u) ** 2
            delta = delta / energy
            urms = urms / time_average_counter
            IF (rank == 0) THEN
!!$          !Assume that Rdisk=0.9
!!$          WRITE(97,109)   '##Final: t,time_average_counter,kin_e(bar{u}),delta, URMS,URMS/Vmax: '
!!$          WRITE(97,110)    time,time_average_counter, energy,delta,  &
!!$               urms, urms/(disk_r*user%solid_vel) 
              WRITE(97,109)   '##Final: t,time_average_counter,kin_e(bar{u}),delta: '
              WRITE(97,110)    time,time_average_counter, energy,delta
            END IF
        END IF !nst or mhd
    END SUBROUTINE my_final_processing


    !DCQ, for saving u_average and p_average
    SUBROUTINE my_write_restart_ns(communicator, vv_mesh, pp_mesh, time, list_mode, &
            un, un_m1, pn, pn_m1, incpn, incpn_m1, filename, it, freq_restart, &
            opt_tempn, opt_tempn_m1, opt_level_set, opt_level_set_m1, opt_max_vel, opt_mono)! HF june 2020!, opt_DR)

        USE def_type_mesh
        USE chaine_caractere
        IMPLICIT NONE
!!        INCLUDE 'mpif.h'
        TYPE(mesh_type), TARGET :: vv_mesh, pp_mesh
        REAL(KIND = 8), INTENT(IN) :: time
        INTEGER, DIMENSION(:), INTENT(IN) :: communicator
        INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: un, un_m1
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: pn, pn_m1, incpn, incpn_m1
        REAL(KIND = 8), DIMENSION(:, :, :), OPTIONAL, INTENT(IN) :: opt_tempn, opt_tempn_m1
        REAL(KIND = 8), DIMENSION(:, :, :, :), OPTIONAL, INTENT(IN) :: opt_level_set, opt_level_set_m1
        REAL(KIND = 8), OPTIONAL, INTENT(IN) :: opt_max_vel
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_mono
        !HF june 2020 REAL(KIND=8), DIMENSION(:,:,:),   OPTIONAL,     INTENT(IN) :: opt_DR
        CHARACTER(len = 200), INTENT(IN) :: filename
        INTEGER, INTENT(IN) :: it, freq_restart
        INTEGER :: code, n, i, rang_S, rang_F, nb_procs_S, nb_procs_F
        INTEGER :: l, lblank
        CHARACTER(len = 3) :: tit, tit_S
        LOGICAL :: mono = .FALSE.
        LOGICAL :: skip
        CHARACTER(len = 250) :: out_name

        CALL MPI_COMM_SIZE(communicator(1), nb_procs_S, code)
        CALL MPI_COMM_SIZE(communicator(2), nb_procs_F, code)
        CALL MPI_COMM_RANK(communicator(1), rang_S, code)
        CALL MPI_COMM_RANK(communicator(2), rang_F, code)

        WRITE(tit, '(i3)') it / freq_restart
        lblank = eval_blank(3, tit)
        DO l = 1, lblank - 1
            tit(l : l) = '0'
        END DO
        WRITE(tit_S, '(i3)') rang_S
        lblank = eval_blank(3, tit_S)
        DO l = 1, lblank - 1
            tit_S(l : l) = '0'
        END DO

        IF (PRESENT(opt_mono)) THEN
            mono = opt_mono
        END IF

        IF (mono) THEN
            !DCQ
            out_name = 'suite_VKS_ns_I' // tit // '.' // filename
        ELSE
            out_name = 'suite_VKS_ns_S' // tit_S // '_I' // tit // '.' // filename
        END IF

        skip = (mono .AND. rang_S /= 0)

        DO n = 1, nb_procs_F
            IF ((rang_F == n - 1) .AND. (.NOT. skip)) THEN
                IF (rang_F == 0) THEN
                    OPEN(UNIT = 10, FILE = out_name, POSITION = 'append', &
                            FORM = 'unformatted', STATUS = 'replace')
                    IF (mono) THEN
                        WRITE(10) time, vv_mesh%np, pp_mesh%np, nb_procs_F, SIZE(list_mode)
                    ELSE
                        WRITE(10) time, nb_procs_S, nb_procs_F, SIZE(list_mode)
                    END IF
                ELSE
                    OPEN(UNIT = 10, FILE = out_name, POSITION = 'append', &
                            FORM = 'unformatted', STATUS = 'unknown')
                END IF

                DO i = 1, SIZE(list_mode)
                    WRITE(10) list_mode(i)
                    WRITE(10) un(:, :, i)
                    WRITE(10) un_m1(:, :, i)
                    WRITE(10) pn(:, :, i)
                    WRITE(10) pn_m1(:, :, i)
                    WRITE(10) incpn(:, :, i)
                    WRITE(10) incpn_m1(:, :, i)
                    IF (PRESENT(opt_tempn) .AND. PRESENT(opt_tempn_m1)) THEN
                        WRITE(10) opt_tempn(:, :, i)
                        WRITE(10) opt_tempn_m1(:, :, i)
                    END IF
                    IF (PRESENT(opt_level_set) .AND. PRESENT(opt_level_set_m1)) THEN
                        WRITE(10) opt_level_set(:, :, :, i)
                        WRITE(10) opt_level_set_m1(:, :, :, i)
                        WRITE(10) opt_max_vel
                    END IF
                    ! HF june 2020
                    !IF (PRESENT(opt_DR)) THEN
                    !   WRITE(10) opt_DR(:,:,i)
                    !END IF
                END DO
                CLOSE(10)
            END IF
            CALL MPI_BARRIER(communicator(2), code)
        END DO

    END SUBROUTINE my_write_restart_ns

    !DCQ, for reading u_average and p_average
    SUBROUTINE my_read_restart_ns(communicator, vv_mesh, pp_mesh, time, list_mode, &
            un, un_m1, pn, pn_m1, incpn, incpn_m1, filename, val_init, interpol, &
            opt_tempn, opt_tempn_m1, opt_level_set, opt_level_set_m1, opt_max_vel, opt_mono)!HF june 2020, opt_DR)

        USE def_type_mesh
        USE chaine_caractere
        USE my_util
        IMPLICIT NONE
!!        INCLUDE 'mpif.h'
        TYPE(mesh_type), TARGET :: vv_mesh, pp_mesh
        REAL(KIND = 8), INTENT(OUT) :: time
        INTEGER, DIMENSION(:), INTENT(IN) :: communicator
        INTEGER, DIMENSION(:) :: list_mode
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(OUT) :: un, un_m1
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(OUT) :: pn, pn_m1, incpn, incpn_m1
        REAL(KIND = 8), DIMENSION(:, :, :), OPTIONAL, INTENT(OUT) :: opt_tempn, opt_tempn_m1
        REAL(KIND = 8), DIMENSION(:, :, :, :), OPTIONAL, INTENT(OUT) :: opt_level_set, opt_level_set_m1
        REAL(KIND = 8), OPTIONAL, INTENT(OUT) :: opt_max_vel
        CHARACTER(len = 200), INTENT(IN) :: filename
        REAL(KIND = 8), OPTIONAL, INTENT(IN) :: val_init
        LOGICAL, OPTIONAL, INTENT(IN) :: interpol
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_mono
        !HF june 2020 REAL(KIND=8), DIMENSION(:,:,:),   OPTIONAL,     INTENT(OUT):: opt_DR
        REAL(KIND = 8) :: max_vel_loc
        INTEGER :: code, n, i, mode, j, rang_S, nb_procs_S, rang_F, nb_procs_F, nlignes, rank
        INTEGER :: m_max_cr, nb_procs_r, nb_procs_Sr
        INTEGER :: m_max_c, nb_mode_r, mode_cherche
        LOGICAL :: trouve, okay
        INTEGER :: npv, npp
        INTEGER :: l, lblank
        CHARACTER(len = 3) :: tit_S
        LOGICAL :: mono = .FALSE.
        CHARACTER(len = 250) :: in_name
        CALL MPI_COMM_RANK(communicator(2), rang_F, code)
        CALL MPI_COMM_SIZE(communicator(2), nb_procs_F, code)
        CALL MPI_COMM_RANK(communicator(1), rang_S, code)
        CALL MPI_COMM_SIZE(communicator(1), nb_procs_S, code)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)

        max_vel_loc = 0.d0

        nlignes = 6
        IF (PRESENT(opt_tempn) .AND. PRESENT(opt_tempn_m1)) THEN
            nlignes = nlignes + 2
        END IF
        IF (PRESENT(opt_level_set) .AND. PRESENT(opt_level_set_m1)) THEN
            nlignes = nlignes + 3
        END IF
        !HF june 2020
        !IF (PRESENT(opt_DR)) THEN
        !    nlignes = nlignes + 1
        !END IF

        WRITE(tit_S, '(i3)') rang_S
        lblank = eval_blank(3, tit_S)
        DO l = 1, lblank - 1
            tit_S(l : l) = '0'
        END DO

        IF (PRESENT(opt_mono)) THEN
            mono = opt_mono
        END IF
        !DCQ
        IF (mono) THEN
            in_name = 'suite_VKS_ns.' // filename
        ELSE
            in_name = 'suite_VKS_ns_S' // tit_S // '.' // filename
        END IF

        WRITE(*, *) 'restart VKS Navier-Stokes'
        OPEN(UNIT = 10, FILE = in_name, FORM = 'unformatted', STATUS = 'unknown')

        IF (mono) THEN
            READ(10) time, npv, npp, nb_procs_r, m_max_cr
        ELSE
            READ(10) time, nb_procs_Sr, nb_procs_r, m_max_cr
        END IF
        CLOSE(10)

        IF ((nb_procs_Sr /= nb_procs_S) .AND. (.NOT. mono)) THEN
            CALL error_petsc('BUG in read_restart: nb_procs_Sr /= nb_procs_S')
            !STOP
        END IF

        IF (rang_F == 0) THEN
            WRITE(*, *) 'File name', TRIM(ADJUSTL(in_name))
            WRITE(*, *) 'Time = ', time
            WRITE(*, *) 'Number of processors from restart file = ', nb_procs_r
            WRITE(*, *) 'Number of modes per processor from restart file = ', m_max_cr
        ENDIF

        m_max_c = SIZE(list_mode)      !nombre de modes par proc pour le calcul
        nb_mode_r = nb_procs_r * m_max_cr  !nombre total de modes contenus dans le suite

        !June 7 2007, JLG
        IF (nb_procs_F * m_max_c /= nb_mode_r) THEN
            !CALL error_petsc('Bug in read_restart_ns: nb_procs_F*m_max_c /= nb_mode_r')
            WRITE(*, *) 'Warning in read_restart_ns: nb_procs_F*m_max_c /= nb_mode_r'
            !STOP
        END IF

        okay = .FALSE.
        IF (PRESENT(interpol)) THEN
            IF (interpol) THEN
                okay = .TRUE.
            END IF
        END IF
        !June 7 2007, JLG

        IF (rank==0) THEN
            WRITE(*, *) 'Reading Navier-Stokes modes ...'
        END IF
        DO i = 1, m_max_c                  !pour tout les modes du processeur courant
            !ouverture du fichier
            OPEN(UNIT = 10, FILE = in_name, FORM = 'unformatted', STATUS = 'unknown')
            !on saute la premiere ligne du fichier qui contient des donnees
            READ(10)
            mode_cherche = list_mode(i)
            !recherche du bon mode
            trouve = .FALSE.
            DO j = 1, nb_mode_r             !pour tout les modes ecris dans le suite.
                !lecture du mode
                READ(10) mode
                !June 7 2007, JLG
                IF (okay) THEN
                    IF (j/=rang_F * m_max_c + i) THEN
                        DO n = 1, nlignes
                            READ(10)
                        ENDDO
                        CYCLE
                    ELSE
                        list_mode(i) = mode
                        mode_cherche = mode
                    END IF
                END IF
                !June 7 2007, JLG
                IF (mode == mode_cherche) THEN   !on a trouve le bon mode
                    READ(10) un(:, :, i)
                    READ(10) un_m1(:, :, i)
                    READ(10) pn(:, :, i)
                    READ(10) pn_m1(:, :, i)
                    READ(10) incpn(:, :, i)
                    READ(10) incpn_m1(:, :, i)
                    IF (PRESENT(opt_tempn) .AND. PRESENT(opt_tempn_m1)) THEN
                        READ(10) opt_tempn(:, :, i)
                        READ(10) opt_tempn_m1(:, :, i)
                    END IF
                    IF (PRESENT(opt_level_set) .AND. PRESENT(opt_level_set_m1)) THEN
                        READ(10) opt_level_set(:, :, :, i)
                        READ(10) opt_level_set_m1(:, :, :, i)
                        READ(10) max_vel_loc
                    END IF
                    !HF june 2020
                    !IF (PRESENT(opt_DR)) THEN
                    !   READ(10) opt_DR(:,:,i)
                    !END IF
                    WRITE(*, '(A,i4,A)') 'mode ns ', mode_cherche, ' found '
                    trouve = .TRUE.
                    EXIT                        !car on a trouve le bon mode
                ELSE                             !on passe au mode suivant en sautant 6 lignes
                    DO n = 1, nlignes
                        READ(10)
                    ENDDO
                ENDIF
            ENDDO

            IF (.NOT.trouve) THEN               !mode_cherche non trouve
                IF (PRESENT(val_init)) THEN ! not implemented yet
                    un(:, :, i) = val_init  ; un_m1(:, :, i) = val_init
                    pn(:, :, i) = val_init ; pn_m1(:, :, i) = val_init
                    incpn(:, :, i) = val_init ; incpn_m1(:, :, i) = val_init
                    IF (PRESENT(opt_tempn) .AND. PRESENT(opt_tempn_m1)) THEN
                        opt_tempn(:, :, i) = val_init
                        opt_tempn_m1(:, :, i) = val_init
                    END IF
                    IF (PRESENT(opt_level_set) .AND. PRESENT(opt_level_set_m1)) THEN
                        opt_level_set(:, :, :, i) = val_init
                        opt_level_set_m1(:, :, :, i) = val_init
                        max_vel_loc = val_init
                    END IF
                    WRITE(*, '(A,i4,A)') 'mode ns', mode_cherche, ' not found'
                ELSE
                    un(:, :, i) = 0.d0 ; un_m1(:, :, i) = 0.d0
                    pn(:, :, i) = 0.d0 ; pn_m1(:, :, i) = 0.d0
                    incpn(:, :, i) = 0.d0 ; incpn_m1(:, :, i) = 0.d0
                    IF (PRESENT(opt_tempn) .AND. PRESENT(opt_tempn_m1)) THEN
                        opt_tempn(:, :, i) = 0.d0
                        opt_tempn_m1(:, :, i) = 0.d0
                    END IF
                    IF (PRESENT(opt_level_set) .AND. PRESENT(opt_level_set_m1)) THEN
                        opt_level_set(:, :, :, i) = 0.d0
                        opt_level_set_m1(:, :, :, i) = 0.d0
                    END IF
                    WRITE(*, *) 'mode ns', mode_cherche, ' not found'
                ENDIF
            ENDIF
            CLOSE(10)                          !fermeture du fichier suite
        ENDDO

        IF (PRESENT(opt_max_vel)) THEN
            CALL MPI_ALLREDUCE(max_vel_loc, opt_max_vel, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, communicator(2), code)
        END IF

    END SUBROUTINE my_read_restart_ns

    SUBROUTINE my_write_restart_maxwell(communicator, H_mesh, phi_mesh, time, list_mode, Hn, Hn1, Bn, Bn1, phin, phin1, &
            filename, it, freq_restart, opt_mono)

        USE def_type_mesh
        USE chaine_caractere
        IMPLICIT NONE
!!        INCLUDE 'mpif.h'
        TYPE(mesh_type), TARGET :: H_mesh, phi_mesh
        INTEGER, DIMENSION(:), INTENT(IN) :: communicator
        REAL(KIND = 8), INTENT(IN) :: time
        INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: Hn, Hn1, Bn, Bn1
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: phin, phin1
        CHARACTER(len = 200), INTENT(IN) :: filename
        INTEGER, INTENT(IN) :: it, freq_restart
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_mono

        INTEGER :: rang_S, rang_F, code, nb_procs_S, nb_procs_F, n, i
        INTEGER :: l, lblank
        CHARACTER(len = 3) :: tit, tit_S
        CHARACTER(len = 250) :: out_name
        LOGICAL :: mono = .FALSE.
        LOGICAL :: skip

        CALL MPI_COMM_SIZE(communicator(1), nb_procs_S, code)
        CALL MPI_COMM_SIZE(communicator(2), nb_procs_F, code)
        CALL MPI_COMM_RANK(communicator(1), rang_S, code)
        CALL MPI_COMM_RANK(communicator(2), rang_F, code)

        WRITE(tit, '(i3)') it / freq_restart
        lblank = eval_blank(3, tit)
        DO l = 1, lblank - 1
            tit(l : l) = '0'
        END DO
        WRITE(tit_S, '(i3)') rang_S
        lblank = eval_blank(3, tit_S)
        DO l = 1, lblank - 1
            tit_S(l : l) = '0'
        END DO

        IF (PRESENT(opt_mono)) THEN
            mono = opt_mono
        END IF

        IF (mono) THEN
            out_name = 'suite_VKS_maxwell_I' // tit // '.' // filename
        ELSE
            out_name = 'suite_VKS_maxwell_S' // tit_S // '_I' // tit // '.' // filename
        END IF
        skip = (mono .AND. rang_S /= 0)

        DO n = 1, nb_procs_F
            IF ((rang_F == n - 1) .AND. (.NOT. skip)) THEN
                IF (rang_F == 0) THEN
                    OPEN(UNIT = 10, FILE = out_name, POSITION = 'append', &
                            FORM = 'unformatted', STATUS = 'replace')
                    IF (mono) THEN
                        WRITE(10) time, H_mesh%np, phi_mesh%np, nb_procs_F, SIZE(list_mode)
                    ELSE
                        WRITE(10) time, nb_procs_S, nb_procs_F, SIZE(list_mode)
                    END IF
                ELSE
                    OPEN(UNIT = 10, FILE = out_name, POSITION = 'append', &
                            FORM = 'unformatted', STATUS = 'unknown')
                END IF
                DO i = 1, SIZE(list_mode)
                    WRITE(10) list_mode(i)
                    IF (H_mesh%me /=0) THEN
                        WRITE(10) Hn(:, :, i)
                        WRITE(10) Hn1(:, :, i)
                        WRITE(10) Bn(:, :, i)
                        WRITE(10) Bn1(:, :, i)
                    ELSE
                        WRITE(10) 1
                        WRITE(10) 1
                        WRITE(10) 1
                        WRITE(10) 1
                    END IF
                    IF (phi_mesh%me /=0) THEN
                        WRITE(10) phin(:, :, i)
                        WRITE(10) phin1(:, :, i)
                    ELSE
                        WRITE(10) 1
                        WRITE(10) 1
                    END IF
                END DO
                CLOSE(10)
            END IF
            CALL MPI_BARRIER(communicator(2), code)
        END DO

    END SUBROUTINE my_write_restart_maxwell

    SUBROUTINE my_read_restart_maxwell(communicator, H_mesh, phi_mesh, time, list_mode, Hn, Hn1, Bn, Bn1, phin, phin1, &
            filename, val_init, interpol, opt_mono)
        USE def_type_mesh
        USE chaine_caractere
        USE my_util
        IMPLICIT NONE
!!        INCLUDE 'mpif.h'
        TYPE(mesh_type), TARGET :: H_mesh, phi_mesh
        INTEGER, DIMENSION(:), INTENT(IN) :: communicator
        REAL(KIND = 8), INTENT(OUT) :: time
        INTEGER, DIMENSION(:) :: list_mode
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(OUT) :: Hn, Hn1, Bn, Bn1
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(OUT) :: phin, phin1
        CHARACTER(len = 200), INTENT(IN) :: filename
        REAL(KIND = 8), OPTIONAL, INTENT(IN) :: val_init
        LOGICAL, OPTIONAL, INTENT(IN) :: interpol
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_mono
        INTEGER :: code, n, i, mode, j, rang_S, rang_F, nb_procs_F, nb_procs_S
        INTEGER :: m_max_cr, nb_procs_r, nb_procs_Sr
        INTEGER :: m_max_c, nb_mode_r, mode_cherche
        LOGICAL :: trouve, okay
        INTEGER :: nph, npp
        INTEGER :: l, lblank
        CHARACTER(len = 3) :: tit_S
        CHARACTER(len = 250) :: in_name
        LOGICAL :: mono = .FALSE.

        CALL MPI_COMM_RANK(communicator(2), rang_F, code)
        CALL MPI_COMM_SIZE(communicator(2), nb_procs_F, code)
        CALL MPI_COMM_RANK(communicator(1), rang_S, code)
        CALL MPI_COMM_SIZE(communicator(1), nb_procs_S, code)

        WRITE(tit_S, '(i3)') rang_S
        lblank = eval_blank(3, tit_S)
        DO l = 1, lblank - 1
            tit_S(l : l) = '0'
        END DO
        IF (PRESENT(opt_mono)) THEN
            mono = opt_mono
        END IF

        IF (mono) THEN
            in_name = 'suite_VKS_maxwell.' // filename
        ELSE
            in_name = 'suite_VKS_maxwell_S' // tit_S // '.' // filename
        END IF

        WRITE(*, *) 'restart VKS Maxwell'
        OPEN(UNIT = 10, FILE = in_name, FORM = 'unformatted', STATUS = 'unknown')
        IF (mono) THEN
            READ(10) time, nph, npp, nb_procs_r, m_max_cr
        ELSE
            READ(10) time, nb_procs_Sr, nb_procs_r, m_max_cr
        END IF

        IF ((nb_procs_Sr /= nb_procs_S) .AND. (.NOT. mono)) THEN
            CALL error_petsc('BUG in read_restart: nb_procs_Sr /= nb_procs_S')
            !STOP
        END IF

        CLOSE(10)

        IF (rang_F == 0) THEN
            WRITE(*, *) 'proprietes fichier ', in_name
            WRITE(*, *) 'time =', time
            WRITE(*, *) 'nombre de processeurs = ', nb_procs_r
            WRITE(*, *) 'nombre de modes par processeur = ', m_max_cr
        ENDIF

        m_max_c = SIZE(list_mode)      !nombre de modes par proc pour le calcul
        nb_mode_r = nb_procs_r * m_max_cr  !nombre total de modes contenus dans le suite

        !June 7 2007, JLG
        IF (nb_procs_F * m_max_c /= nb_mode_r) THEN
            WRITE(*, *) ' BUG '
            !STOP
        END IF

        okay = .FALSE.
        IF (PRESENT(interpol)) THEN
            IF (interpol) THEN
                okay = .TRUE.
            END IF
        END IF
        !June 7 2007, JLG

        WRITE(*, *) 'Reading Maxwell modes ...'
        DO i = 1, m_max_c                  !pour tout les modes du processeur courant
            !ouverture du fichier
            OPEN(UNIT = 10, FILE = in_name, FORM = 'unformatted', STATUS = 'unknown')
            !on saute la premier ligne du fichier qui contient des donnes
            READ(10)
            mode_cherche = list_mode(i)
            !recherche du bon mode
            trouve = .FALSE.
            DO j = 1, nb_mode_r             !pour tout les modes ecris dans le suite.
                !lecture du mode
                READ(10) mode
                !June 7 2007, JLG
                IF (okay) THEN
                    IF (j/=rang_F * m_max_c + i) THEN
                        DO n = 1, 6
                            READ(10)
                        ENDDO
                        CYCLE
                    ELSE
                        list_mode(i) = mode
                        mode_cherche = mode
                    END IF
                END IF
                !June 7 2007, JLG
                IF (mode == mode_cherche) THEN   !on a trouve le bon mode
                    IF (H_mesh%me /=0) THEN
                        READ(10) Hn(:, :, i)
                        READ(10) Hn1(:, :, i)
                        READ(10) Bn(:, :, i)
                        READ(10) Bn1(:, :, i)
                    ELSE
                        READ(10)
                        READ(10)
                        READ(10)
                        READ(10)
                    END IF
                    IF (phi_mesh%me /=0) THEN
                        READ(10) phin(:, :, i)
                        READ(10) phin1(:, :, i)
                    ELSE
                        READ(10)
                        READ(10)
                    END IF
                    WRITE(*, *) 'mode maxwell', mode_cherche, ' trouve '
                    trouve = .TRUE.
                    EXIT                        !car on a trouve le bon mode
                ELSE                             !on passe au mode suivant en sautant 4 lignes
                    DO n = 1, 6
                        READ(10)
                    ENDDO
                ENDIF
            ENDDO
            IF (.NOT.trouve) THEN               !mode_cherche non trouve
                IF (PRESENT(val_init)) THEN
                    Hn(:, :, i) = val_init ; Hn1(:, :, i) = val_init
                    phin(:, :, i) = val_init ; phin1(:, :, i) = val_init
                    WRITE(*, *) 'mode maxwell', mode_cherche, ' non trouve'
                ELSE
                    Hn(:, :, i) = 0.d0 ; Hn1(:, :, i) = 0.d0
                    phin(:, :, i) = 0.d0 ; phin1(:, :, i) = 0.d0
                    WRITE(*, *) 'mode maxwell', mode_cherche, ' non trouve'
                ENDIF
            ENDIF
            CLOSE(10)                          !fermeture du fichier suite
        ENDDO

    END SUBROUTINE my_read_restart_maxwell

    SUBROUTINE compute_helicity_Vel(comm_one_d_ns, vv_mesh, list_mode, un, helicity, opt_rr, opt_phys)
        USE tn_axi
        USE boundary
        USE user_data
        USE sfemans_tools
        USE sft_parallele
        USE my_util
        IMPLICIT NONE
        MPI_Comm, DIMENSION(:), POINTER :: comm_one_d_ns
        TYPE(mesh_type), TARGET :: vv_mesh
        INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: un
        REAL(KIND = 8), DIMENSION(:, :), OPTIONAL, INTENT(IN) :: opt_rr
        REAL(KIND = 8), DIMENSION(:, :, :), OPTIONAL, INTENT(INOUT) :: opt_phys
        REAL(KIND = 8), INTENT(OUT) :: helicity
        REAL(KIND = 8), DIMENSION(3) :: umax_loc, umax
        INTEGER :: l, i, k
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%me, 6, SIZE(list_mode)) :: RotV, V_gauss
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%me, 2, SIZE(list_mode)) :: u_rotu_gauss
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%n_w, 6) :: Vs
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%k_d, vv_mesh%gauss%n_w) :: dw_loc
        REAL(KIND = 8), DIMENSION(H_mesh%gauss%k_d, H_mesh%gauss%n_w) :: dw_loc_h
        INTEGER :: index
        REAL(KIND = 8) :: inte_hel_loc, inte_hel_tot_F, inte_hel_tot
        REAL(KIND = 8) :: ray
        INTEGER, DIMENSION(vv_mesh%gauss%n_w) :: j_loc
        INTEGER, SAVE :: nb_procs, bloc_size_gauss
        INTEGER :: code, m, mode, rank_S, rank_F
        INTEGER :: bloc_size, m_max_pad, int_nb, n, np_gauss_for_phys
        REAL(KIND = 8), PARAMETER :: pi = 3.14159265358979323846d0

        CALL MPI_Comm_rank(comm_one_d_ns(1), rank_S, code)
        CALL MPI_Comm_rank(comm_one_d_ns(2), rank_F, code)

        !==Compute helicity = integration on volume of U.ROT(U)
        DO i = 1, SIZE(list_mode)
            mode = list_mode(i)
            index = 0
            DO m = 1, vv_mesh%me
                j_loc = vv_mesh%jj(:, m)
                DO k = 1, 6
                    Vs(:, k) = un(j_loc, k, i)
                END DO

                DO l = 1, vv_mesh%gauss%l_G
                    index = index + 1
                    dw_loc = vv_mesh%gauss%dw(:, :, l, m)

                    !===Compute radius of Gauss point
                    ray = SUM(vv_mesh%rr(1, j_loc) * vv_mesh%gauss%ww(:, l))

                    !-----------------vitesse sur les points de Gauss---------------------------
                    V_gauss(index, 1, i) = SUM(Vs(:, 1) * vv_mesh%gauss%ww(:, l))
                    V_gauss(index, 3, i) = SUM(Vs(:, 3) * vv_mesh%gauss%ww(:, l))
                    V_gauss(index, 5, i) = SUM(Vs(:, 5) * vv_mesh%gauss%ww(:, l))

                    V_gauss(index, 2, i) = SUM(Vs(:, 2) * vv_mesh%gauss%ww(:, l))
                    V_gauss(index, 4, i) = SUM(Vs(:, 4) * vv_mesh%gauss%ww(:, l))
                    V_gauss(index, 6, i) = SUM(Vs(:, 6) * vv_mesh%gauss%ww(:, l))
                    !-----------------rotational sur les points de Gauss---------------------------
                    !coeff sur les cosinus
                    RotV(index, 1, i) = mode / ray * V_gauss(index, 6, i) &
                            - SUM(Vs(:, 3) * dw_loc(2, :))
                    RotV(index, 4, i) = SUM(Vs(:, 2) * dw_loc(2, :)) &
                            - SUM(Vs(:, 6) * dw_loc(1, :))
                    RotV(index, 5, i) = 1 / ray * V_gauss(index, 3, i) &
                            + SUM(Vs(:, 3) * dw_loc(1, :)) &
                            - mode / ray * V_gauss(index, 2, i)
                    !coeff sur les sinus
                    RotV(index, 2, i) = -mode / ray * V_gauss(index, 5, i) &
                            - SUM(Vs(:, 4) * dw_loc(2, :))
                    RotV(index, 3, i) = SUM(Vs(:, 1) * dw_loc(2, :)) &
                            - SUM(Vs(:, 5) * dw_loc(1, :))
                    RotV(index, 6, i) = 1 / ray * V_gauss(index, 4, i) &
                            + SUM(Vs(:, 4) * dw_loc(1, :))&
                            + mode / ray * V_gauss(index, 1, i)
                ENDDO
            ENDDO
        END DO

        CALL MPI_COMM_SIZE(comm_one_d_ns(2), nb_procs, code)
        bloc_size = SIZE(RotV, 1) / nb_procs + 1
        m_max_pad = 3 * SIZE(list_mode) * nb_procs / 2
        CALL FFT_PAR_DOT_PROD_DCL(comm_one_d_ns(2), RotV, V_gauss, u_rotu_gauss, nb_procs, bloc_size, m_max_pad)
        bloc_size_gauss = vv_mesh%gauss%l_G * vv_mesh%dom_me / nb_procs + 1
        IF (PRESENT(opt_phys)) THEN
            IF (.NOT.(PRESENT(opt_rr))) CALL error_petsc('Wrong optionals in compute helicity_Vel')
            !CALL FFT_PAR_REAL_GAUSS_PDF(comm_one_d_ns(1), comm_one_d_ns(2), u_rotu_gauss, opt_rr, &
            !            opt_phys, 'hel', nb_procs, bloc_size_gauss)!, rank_ns_S, rank_ns_F)
!TEST LC DEBUG 2024 not in fft
!            CALL FFT_REAL(comm_one_d_ns(2), u_rotu_gauss, opt_phys, np_gauss_for_phys, &
!                    nb_procs, bloc_size_gauss)
!TEST LC DEBUG 2024 not in fft
            !CALL WRITE_REAL_GAUSS_PDF(rank_S, rank_F, opt_phys, &
            !        np_gauss_for_phys, opt_rr, 'hel', nb_procs, bloc_size_gauss)
        ENDIF
        !===Integration on the volume of the helicity
        inte_hel_loc = 0.d0
        DO i = 1, SIZE(list_mode)
            IF (list_mode(i)==0) THEN
                index = 0
                DO m = 1, vv_mesh%me
                    j_loc = vv_mesh%jj(:, m)
                    DO l = 1, vv_mesh%gauss%l_G
                        index = index + 1
                        !===Compute radius of Gauss point
                        ray = SUM(vv_mesh%rr(1, j_loc) * vv_mesh%gauss%ww(:, l))
                        inte_hel_loc = inte_hel_loc + u_rotu_gauss(index, 1, i) * ray * vv_mesh%gauss%rj(l, m)
                    END DO
                END DO
            END IF
        END DO
        inte_hel_loc = inte_hel_loc * 2 * pi
        CALL MPI_ALLREDUCE(inte_hel_loc, inte_hel_tot_F, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                comm_one_d_ns(2), code)
        CALL MPI_ALLREDUCE(inte_hel_tot_F, inte_hel_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                comm_one_d_ns(1), code)

        !helicity
        helicity = inte_hel_tot
    END SUBROUTINE compute_helicity_Vel

    SUBROUTINE compute_helicity_Mag(comm_one_d, H_mesh, list_mode, Hn, Bn, helicity)
        USE tn_axi
        USE boundary
        USE user_data
        USE sfemans_tools
        USE sft_parallele
        USE my_util
        IMPLICIT NONE
        MPI_Comm, DIMENSION(:), POINTER :: comm_one_d
        TYPE(mesh_type), TARGET :: H_mesh
        INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: Hn, Bn
        REAL(KIND = 8), INTENT(OUT) :: helicity
        REAL(KIND = 8), DIMENSION(3) :: umax_loc, umax
        INTEGER :: l, i, k
        REAL(KIND = 8), DIMENSION(H_mesh%gauss%l_G * H_mesh%me, 6, SIZE(list_mode)) :: RotH, H_gauss, B_gauss
        REAL(KIND = 8), DIMENSION(H_mesh%gauss%l_G * H_mesh%me, 2, SIZE(list_mode)) :: B_rotH_gauss
        REAL(KIND = 8), DIMENSION(H_mesh%gauss%n_w, 6) :: Hs, Bs
        REAL(KIND = 8), DIMENSION(H_mesh%gauss%k_d, H_mesh%gauss%n_w) :: dw_loc
        REAL(KIND = 8), DIMENSION(H_mesh%gauss%k_d, H_mesh%gauss%n_w) :: dw_loc_h
        INTEGER :: index
        REAL(KIND = 8) :: inte_hel_loc, inte_hel_tot_F, inte_hel_tot
        REAL(KIND = 8) :: ray
        INTEGER, DIMENSION(H_mesh%gauss%n_w) :: j_loc
        INTEGER, SAVE :: nb_procs
        INTEGER :: code, m, mode
        INTEGER :: bloc_size, m_max_pad, int_nb, n
        REAL(KIND = 8), PARAMETER :: pi = 3.14159265358979323846d0

        !==Compute helicity = integration on volume of B.ROT(H)
        DO i = 1, SIZE(list_mode)
            mode = list_mode(i)
            index = 0
            DO m = 1, H_mesh%me
                j_loc = H_mesh%jj(:, m)
                DO k = 1, 6
                    Hs(:, k) = Hn(j_loc, k, i)
                    Bs(:, k) = Bn(j_loc, k, i)
                END DO
                DO l = 1, H_mesh%gauss%l_G
                    index = index + 1
                    dw_loc = H_mesh%gauss%dw(:, :, l, m)
                    !===Compute radius of Gauss point
                    ray = SUM(H_mesh%rr(1, j_loc) * H_mesh%gauss%ww(:, l))
                    !-----------------vitesse sur les points de Gauss---------------------------
                    H_gauss(index, 1, i) = SUM(Hs(:, 1) * H_mesh%gauss%ww(:, l))
                    H_gauss(index, 3, i) = SUM(Hs(:, 3) * H_mesh%gauss%ww(:, l))
                    H_gauss(index, 5, i) = SUM(Hs(:, 5) * H_mesh%gauss%ww(:, l))

                    H_gauss(index, 2, i) = SUM(Hs(:, 2) * H_mesh%gauss%ww(:, l))
                    H_gauss(index, 4, i) = SUM(Hs(:, 4) * H_mesh%gauss%ww(:, l))
                    H_gauss(index, 6, i) = SUM(Hs(:, 6) * H_mesh%gauss%ww(:, l))

                    B_gauss(index, 1, i) = SUM(Bs(:, 1) * H_mesh%gauss%ww(:, l))
                    B_gauss(index, 3, i) = SUM(Bs(:, 3) * H_mesh%gauss%ww(:, l))
                    B_gauss(index, 5, i) = SUM(Bs(:, 5) * H_mesh%gauss%ww(:, l))

                    B_gauss(index, 2, i) = SUM(Bs(:, 2) * H_mesh%gauss%ww(:, l))
                    B_gauss(index, 4, i) = SUM(Bs(:, 4) * H_mesh%gauss%ww(:, l))
                    B_gauss(index, 6, i) = SUM(Bs(:, 6) * H_mesh%gauss%ww(:, l))
                    !-----------------rotational sur les points de Gauss---------------------------
                    !coeff sur les cosinus
                    RotH(index, 1, i) = mode / ray * H_gauss(index, 6, i) &
                            - SUM(Hs(:, 3) * dw_loc(2, :))
                    RotH(index, 4, i) = SUM(Hs(:, 2) * dw_loc(2, :)) &
                            - SUM(Hs(:, 6) * dw_loc(1, :))
                    RotH(index, 5, i) = 1 / ray * H_gauss(index, 3, i) &
                            + SUM(Hs(:, 3) * dw_loc(1, :)) &
                            - mode / ray * H_gauss(index, 2, i)
                    !coeff sur les sinus
                    RotH(index, 2, i) = -mode / ray * H_gauss(index, 5, i) &
                            - SUM(Hs(:, 4) * dw_loc(2, :))
                    RotH(index, 3, i) = SUM(Hs(:, 1) * dw_loc(2, :)) &
                            - SUM(Hs(:, 5) * dw_loc(1, :))
                    RotH(index, 6, i) = 1 / ray * H_gauss(index, 4, i) &
                            + SUM(Hs(:, 4) * dw_loc(1, :))&
                            + mode / ray * H_gauss(index, 1, i)
                ENDDO
            ENDDO
        END DO

        CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs, code)
        bloc_size = SIZE(RotH, 1) / nb_procs + 1
        m_max_pad = 3 * SIZE(list_mode) * nb_procs / 2
        CALL FFT_PAR_DOT_PROD_DCL(comm_one_d(2), RotH, B_gauss, B_rotH_gauss, nb_procs, bloc_size, m_max_pad)
        !===Integration on the volume of the helicity
        inte_hel_loc = 0.d0
        DO i = 1, SIZE(list_mode)
            IF (list_mode(i)==0) THEN
                index = 0
                DO m = 1, H_mesh%me
                    j_loc = H_mesh%jj(:, m)
                    DO l = 1, H_mesh%gauss%l_G
                        index = index + 1
                        !===Compute radius of Gauss point
                        ray = SUM(H_mesh%rr(1, j_loc) * H_mesh%gauss%ww(:, l))
                        inte_hel_loc = inte_hel_loc + B_rotH_gauss(index, 1, i) * ray * H_mesh%gauss%rj(l, m)
                    END DO
                END DO
            END IF
        END DO
        inte_hel_loc = inte_hel_loc * 2 * pi
        CALL  MPI_ALLREDUCE(inte_hel_loc, inte_hel_tot_F, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                comm_one_d(2), code)
        CALL MPI_ALLREDUCE(inte_hel_tot_F, inte_hel_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                comm_one_d(1), code)
        !helicity
        helicity = inte_hel_tot
    END SUBROUTINE compute_helicity_Mag


    SUBROUTINE compute_enstrophy(time, vv_mesh, communicator, list_mode, un)
        !over domain without obstacles
        USE def_type_mesh
        USE input_data
        USE boundary
        USE sft_parallele
#include "petsc/finclude/petsc.h"
        USE petsc
        IMPLICIT NONE
        TYPE(mesh_type), INTENT(IN) :: vv_mesh
        INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: un
        REAL(KIND = 8), INTENT(IN) :: time
        REAL(KIND = 8), DIMENSION(2, vv_mesh%gauss%l_G * vv_mesh%dom_me) :: rr_gauss
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 6, SIZE(list_mode)) :: RotV, W
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 2, SIZE(list_mode)) :: RotV_sqr
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 2, SIZE(list_mode)) :: RotV_sqr_penal
        INTEGER, DIMENSION(vv_mesh%gauss%n_w) :: j_loc
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%k_d, vv_mesh%gauss%n_w) :: dw_loc
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%n_w, 6) :: Vs
        REAL(KIND = 8) :: enstrophy, enstrophy_tot, enstrophy_loc
        REAL(KIND = 8) :: vol, vol_tot, vol_loc
        REAL(KIND = 8) :: ray
        INTEGER :: m, l, i, mode, index, k, TYPE, nb_procs, m_max_pad, bloc_size
        !v5.0!#include "petsc/finclude/petsc.h"
        PetscErrorCode                   :: ierr
        MPI_Comm,DIMENSION(2) :: communicator


        DO i = 1, SIZE(list_mode)
            mode = list_mode(i)
            index = 0
            DO m = 1, vv_mesh%dom_me
                j_loc = vv_mesh%jj(:, m)
                DO k = 1, 6
                    Vs(:, k) = un(j_loc, k, i)
                END DO
                DO l = 1, vv_mesh%gauss%l_G
                    index = index + 1

                    dw_loc = vv_mesh%gauss%dw(:, :, l, m)

                    rr_gauss(1, index) = SUM(vv_mesh%rr(1, j_loc) * vv_mesh%gauss%ww(:, l))
                    rr_gauss(2, index) = SUM(vv_mesh%rr(2, j_loc) * vv_mesh%gauss%ww(:, l))

                    ray = rr_gauss(1, index)

                    !-----------------vitesse sur les points de Gauss---------------------------
                    W(index, 1, i) = SUM(Vs(:, 1) * vv_mesh%gauss%ww(:, l))
                    W(index, 3, i) = SUM(Vs(:, 3) * vv_mesh%gauss%ww(:, l))
                    W(index, 5, i) = SUM(Vs(:, 5) * vv_mesh%gauss%ww(:, l))

                    W(index, 2, i) = SUM(Vs(:, 2) * vv_mesh%gauss%ww(:, l))
                    W(index, 4, i) = SUM(Vs(:, 4) * vv_mesh%gauss%ww(:, l))
                    W(index, 6, i) = SUM(Vs(:, 6) * vv_mesh%gauss%ww(:, l))
                    !-----------------vecteur rotation sur les points
                    RotV(index, 1, i) = mode / ray * W(index, 6, i) &
                            - SUM(Vs(:, 3) * dw_loc(2, :))
                    RotV(index, 4, i) = SUM(Vs(:, 2) * dw_loc(2, :)) &
                            - SUM(Vs(:, 6) * dw_loc(1, :))
                    RotV(index, 5, i) = 1 / ray * W(index, 3, i) &
                            + SUM(Vs(:, 3) * dw_loc(1, :)) &
                            - mode / ray * W(index, 2, i)

                    RotV(index, 2, i) = -mode / ray * W(index, 5, i) &
                            - SUM(Vs(:, 4) * dw_loc(2, :))
                    RotV(index, 3, i) = SUM(Vs(:, 1) * dw_loc(2, :)) &
                            - SUM(Vs(:, 5) * dw_loc(1, :))
                    RotV(index, 6, i) = 1 / ray * W(index, 4, i) &
                            + SUM(Vs(:, 4) * dw_loc(1, :))&
                            + mode / ray * W(index, 1, i)
                END DO
            END DO
        END DO

        CALL MPI_COMM_SIZE(communicator(2), nb_procs, ierr)
        m_max_pad = 3 * SIZE(list_mode) * nb_procs / 2
        bloc_size = SIZE(RotV, 1) / nb_procs + 1

        CALL FFT_PAR_DOT_PROD_DCL(communicator(2), RotV, RotV, RotV_sqr, nb_procs, &
                bloc_size, m_max_pad)
        CALL FFT_PAR_VAR_ETA_PROD_GAUSS_DCL(communicator(2), penal_in_real_space, vv_mesh, &
                RotV_sqr, RotV_sqr_penal, nb_procs, bloc_size, m_max_pad, rr_gauss, time)

        enstrophy = 0.d0
        enstrophy_tot = 0.d0
        enstrophy_loc = 0.d0
        vol = 0.d0
        vol_tot = 0.d0
        vol_loc = 0.d0
        DO i = 1, SIZE(list_mode)
            mode = list_mode(i)
            IF (mode/=0) THEN
                CYCLE
            ELSE
                index = 0
                DO m = 1, vv_mesh%dom_me
                    DO l = 1, vv_mesh%gauss%l_G
                        index = index + 1
                        !===Enstrophy is int_domain chi*rotv**2 dx; chi=1 in fluid, 0 in solid
                        enstrophy_loc = enstrophy_loc + RotV_sqr_penal(index, 1, i) &
                                * rr_gauss(1, index) * vv_mesh%gauss%rj(l, m)
                        vol_loc = vol_loc + 1.d0 * rr_gauss(1, index) * vv_mesh%gauss%rj(l, m)
                    END DO
                END DO
            END IF
        END DO
        enstrophy_loc = 2 * ACOS(-1.d0) * enstrophy_loc
        vol_loc = 2 * ACOS(-1.d0) * vol_loc
        CALL MPI_ALLREDUCE(enstrophy_loc, enstrophy_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, communicator(1), ierr)
        CALL MPI_ALLREDUCE(enstrophy_tot, enstrophy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, communicator(2), ierr)
        CALL MPI_ALLREDUCE(vol_loc, vol_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, communicator(1), ierr)
        CALL MPI_ALLREDUCE(vol_tot, vol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, communicator(2), ierr)

        IF (rank==0) THEN
            WRITE(*, *) ' Enstrophy ', time, enstrophy
            WRITE(888, *) time, enstrophy, vol
        END IF

    END SUBROUTINE compute_enstrophy

    SUBROUTINE compute_enstrophy_ALL(time, vv_mesh, communicator, list_mode, un)
        ! over the whole domain
        ! enstrophy = \int (rotv)^2, be careful if you want the coefficient 1/2
        USE def_type_mesh
        USE input_data
        USE boundary
        USE sft_parallele
#include "petsc/finclude/petsc.h"
        USE petsc
        IMPLICIT NONE
        TYPE(mesh_type), INTENT(IN) :: vv_mesh
        INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: un
        REAL(KIND = 8), INTENT(IN) :: time
        REAL(KIND = 8), DIMENSION(2, vv_mesh%gauss%l_G * vv_mesh%dom_me) :: rr_gauss
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 6, SIZE(list_mode)) :: RotV, W
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 2, SIZE(list_mode)) :: RotV_sqr
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 2, SIZE(list_mode)) :: RotV_sqr_penal
        INTEGER, DIMENSION(vv_mesh%gauss%n_w) :: j_loc
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%k_d, vv_mesh%gauss%n_w) :: dw_loc
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%n_w, 6) :: Vs
        REAL(KIND = 8) :: enstrophy, enstrophy_tot, enstrophy_loc
        !REAL(KIND=8)      :: vol, vol_tot, vol_loc
        REAL(KIND = 8) :: ray
        INTEGER :: m, l, i, mode, index, k, TYPE, nb_procs, m_max_pad, bloc_size
        !v5.0!#include "petsc/finclude/petsc.h"
        PetscErrorCode                   :: ierr
        MPI_Comm,DIMENSION(2) :: communicator


        DO i = 1, SIZE(list_mode)
            mode = list_mode(i)
            index = 0
            DO m = 1, vv_mesh%dom_me
                j_loc = vv_mesh%jj(:, m)
                DO k = 1, 6
                    Vs(:, k) = un(j_loc, k, i)
                END DO
                DO l = 1, vv_mesh%gauss%l_G
                    index = index + 1

                    dw_loc = vv_mesh%gauss%dw(:, :, l, m)

                    rr_gauss(1, index) = SUM(vv_mesh%rr(1, j_loc) * vv_mesh%gauss%ww(:, l))
                    rr_gauss(2, index) = SUM(vv_mesh%rr(2, j_loc) * vv_mesh%gauss%ww(:, l))

                    ray = rr_gauss(1, index)

                    !-----------------vitesse sur les points de Gauss---------------------------
                    W(index, 1, i) = SUM(Vs(:, 1) * vv_mesh%gauss%ww(:, l))
                    W(index, 3, i) = SUM(Vs(:, 3) * vv_mesh%gauss%ww(:, l))
                    W(index, 5, i) = SUM(Vs(:, 5) * vv_mesh%gauss%ww(:, l))

                    W(index, 2, i) = SUM(Vs(:, 2) * vv_mesh%gauss%ww(:, l))
                    W(index, 4, i) = SUM(Vs(:, 4) * vv_mesh%gauss%ww(:, l))
                    W(index, 6, i) = SUM(Vs(:, 6) * vv_mesh%gauss%ww(:, l))
                    !-----------------vecteur rotation sur les points
                    RotV(index, 1, i) = mode / ray * W(index, 6, i) &
                            - SUM(Vs(:, 3) * dw_loc(2, :))
                    RotV(index, 4, i) = SUM(Vs(:, 2) * dw_loc(2, :)) &
                            - SUM(Vs(:, 6) * dw_loc(1, :))
                    RotV(index, 5, i) = 1 / ray * W(index, 3, i) &
                            + SUM(Vs(:, 3) * dw_loc(1, :)) &
                            - mode / ray * W(index, 2, i)

                    RotV(index, 2, i) = -mode / ray * W(index, 5, i) &
                            - SUM(Vs(:, 4) * dw_loc(2, :))
                    RotV(index, 3, i) = SUM(Vs(:, 1) * dw_loc(2, :)) &
                            - SUM(Vs(:, 5) * dw_loc(1, :))
                    RotV(index, 6, i) = 1 / ray * W(index, 4, i) &
                            + SUM(Vs(:, 4) * dw_loc(1, :))&
                            + mode / ray * W(index, 1, i)
                END DO
            END DO
        END DO

        CALL MPI_COMM_SIZE(communicator(2), nb_procs, ierr)
        m_max_pad = 3 * SIZE(list_mode) * nb_procs / 2
        bloc_size = SIZE(RotV, 1) / nb_procs + 1

        CALL FFT_PAR_DOT_PROD_DCL(communicator(2), RotV, RotV, RotV_sqr, nb_procs, &
                bloc_size, m_max_pad)
        !CALL FFT_PAR_VAR_ETA_PROD_GAUSS_DCL(communicator(2), penal_in_real_space, vv_mesh, &
        !     RotV_sqr, RotV_sqr_penal, nb_procs, bloc_size, m_max_pad, rr_gauss, time)
        !TEST LC: no ETA_FFT (we compute enstrophy in all domain)
        RotV_sqr_penal = RotV_sqr

        enstrophy = 0.d0
        enstrophy_tot = 0.d0
        enstrophy_loc = 0.d0
        !vol = 0.d0
        !vol_tot=0.d0
        !vol_loc=0.d0
        DO i = 1, SIZE(list_mode)
            mode = list_mode(i)
            IF (mode/=0) THEN
                CYCLE
            ELSE
                index = 0
                DO m = 1, vv_mesh%dom_me
                    DO l = 1, vv_mesh%gauss%l_G
                        index = index + 1
                        !===Enstrophy is int_domain chi*rotv**2 dx; chi=1 in fluid, 0 in solid
                        enstrophy_loc = enstrophy_loc + RotV_sqr_penal(index, 1, i) &
                                * rr_gauss(1, index) * vv_mesh%gauss%rj(l, m)
                        !vol_loc =  vol_loc + 1.d0*rr_gauss(1,index)*vv_mesh%gauss%rj(l,m)
                    END DO
                END DO
            END IF
        END DO
        enstrophy_loc = 2 * ACOS(-1.d0) * enstrophy_loc
        !vol_loc=2*ACOS(-1.d0)*vol_loc
        CALL MPI_ALLREDUCE(enstrophy_loc, enstrophy_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, communicator(1), ierr)
        CALL MPI_ALLREDUCE(enstrophy_tot, enstrophy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, communicator(2), ierr)
        !CALL MPI_ALLREDUCE(vol_loc, vol_tot,1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator(1), ierr)
        !CALL MPI_ALLREDUCE(vol_tot, vol, 1,MPI_DOUBLE_PRECISION,MPI_SUM, communicator(2), ierr)

        IF (rank==0) THEN
            WRITE(*, *) ' Enstrophy ', time, enstrophy
            WRITE(889, *) time, enstrophy!, vol
        END IF

    END SUBROUTINE compute_enstrophy_ALL

    SUBROUTINE compute_dissipation(time, vv_mesh, communicator, list_mode, un)
        ! over domain without obstacles
        ! we compute the gradient of velocity and use each line for computing the
        ! dissipation
        USE def_type_mesh
        USE input_data
        USE boundary
        USE sft_parallele
#include "petsc/finclude/petsc.h"
        USE petsc
        IMPLICIT NONE
        TYPE(mesh_type), INTENT(IN) :: vv_mesh
        INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: un
        REAL(KIND = 8), INTENT(IN) :: time
        REAL(KIND = 8), DIMENSION(2, vv_mesh%gauss%l_G * vv_mesh%dom_me) :: rr_gauss
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 6, SIZE(list_mode)) :: grad1_vel
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 6, SIZE(list_mode)) :: grad2_vel
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 6, SIZE(list_mode)) :: grad3_vel
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 2, SIZE(list_mode)) :: grad1_vel_sqr
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 2, SIZE(list_mode)) :: grad2_vel_sqr
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 2, SIZE(list_mode)) :: grad3_vel_sqr
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 2, SIZE(list_mode)) :: grad1_vel_sqr_penal
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 2, SIZE(list_mode)) :: grad2_vel_sqr_penal
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 2, SIZE(list_mode)) :: grad3_vel_sqr_penal
        INTEGER, DIMENSION(vv_mesh%gauss%n_w) :: j_loc
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%k_d, vv_mesh%gauss%n_w) :: dw_loc
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%n_w, 6) :: vel_loc
        REAL(KIND = 8) :: dissipation, dissipation_tot, dissipation_loc
        REAL(KIND = 8) :: ray
        INTEGER :: m, l, i, mode, index, k, TYPE, nb_procs, m_max_pad, bloc_size
        PetscErrorCode                   :: ierr
        MPI_Comm,DIMENSION(2) :: communicator

        DO i = 1, SIZE(list_mode)
            mode = list_mode(i)
            index = 0
            DO m = 1, vv_mesh%dom_me
                j_loc = vv_mesh%jj(:, m)
                DO k = 1, 6
                    vel_loc(:, k) = un(j_loc, k, i)
                END DO
                DO l = 1, vv_mesh%gauss%l_G
                    index = index + 1

                    dw_loc = vv_mesh%gauss%dw(:, :, l, m)

                    !-----------------gauss radius----------------------------------------------
                    rr_gauss(1, index) = SUM(vv_mesh%rr(1, j_loc) * vv_mesh%gauss%ww(:, l))
                    rr_gauss(2, index) = SUM(vv_mesh%rr(2, j_loc) * vv_mesh%gauss%ww(:, l))
                    ray = rr_gauss(1, index)

                    !-----------------Grad u_r on Gauss points------------------------------------
                    grad1_vel(index, 1, i) = SUM(vel_loc(:, 1) * dw_loc(1, :))
                    grad1_vel(index, 2, i) = SUM(vel_loc(:, 2) * dw_loc(1, :))
                    grad1_vel(index, 3, i) = (mode * SUM(vel_loc(:, 2) * vv_mesh%gauss%ww(:, l)) &
                            - SUM(vel_loc(:, 3) * vv_mesh%gauss%ww(:, l))) / ray
                    grad1_vel(index, 4, i) = (-mode * SUM(vel_loc(:, 1) * vv_mesh%gauss%ww(:, l)) &
                            - SUM(vel_loc(:, 4) * vv_mesh%gauss%ww(:, l))) / ray
                    grad1_vel(index, 5, i) = SUM(vel_loc(:, 1) * dw_loc(2, :))
                    grad1_vel(index, 6, i) = SUM(vel_loc(:, 2) * dw_loc(2, :))

                    !-----------------Grad u_th on Gauss points-----------------------------------
                    grad2_vel(index, 1, i) = SUM(vel_loc(:, 3) * dw_loc(1, :))
                    grad2_vel(index, 2, i) = SUM(vel_loc(:, 4) * dw_loc(1, :))
                    grad2_vel(index, 3, i) = (mode * SUM(vel_loc(:, 4) * vv_mesh%gauss%ww(:, l)) &
                            + SUM(vel_loc(:, 1) * vv_mesh%gauss%ww(:, l))) / ray
                    grad2_vel(index, 4, i) = (-mode * SUM(vel_loc(:, 3) * vv_mesh%gauss%ww(:, l)) &
                            + SUM(vel_loc(:, 2) * vv_mesh%gauss%ww(:, l))) / ray
                    grad2_vel(index, 5, i) = SUM(vel_loc(:, 3) * dw_loc(2, :))
                    grad2_vel(index, 6, i) = SUM(vel_loc(:, 4) * dw_loc(2, :))

                    !-----------------Grad u_z on Gauss points------------------------------------
                    grad3_vel(index, 1, i) = SUM(vel_loc(:, 5) * dw_loc(1, :))
                    grad3_vel(index, 2, i) = SUM(vel_loc(:, 6) * dw_loc(1, :))
                    grad3_vel(index, 3, i) = mode * SUM(vel_loc(:, 6) * vv_mesh%gauss%ww(:, l)) / ray
                    grad3_vel(index, 4, i) = -mode * SUM(vel_loc(:, 5) * vv_mesh%gauss%ww(:, l)) / ray
                    grad3_vel(index, 5, i) = SUM(vel_loc(:, 5) * dw_loc(2, :))
                    grad3_vel(index, 6, i) = SUM(vel_loc(:, 6) * dw_loc(2, :))

                END DO
            END DO
        END DO

        CALL MPI_COMM_SIZE(communicator(2), nb_procs, ierr)
        m_max_pad = 3 * SIZE(list_mode) * nb_procs / 2
        bloc_size = SIZE(grad1_vel, 1) / nb_procs + 1

        CALL FFT_PAR_DOT_PROD_DCL(communicator(2), grad1_vel, grad1_vel, grad1_vel_sqr, nb_procs, &
                bloc_size, m_max_pad)
        CALL FFT_PAR_DOT_PROD_DCL(communicator(2), grad2_vel, grad2_vel, grad2_vel_sqr, nb_procs, &
                bloc_size, m_max_pad)
        CALL FFT_PAR_DOT_PROD_DCL(communicator(2), grad3_vel, grad3_vel, grad3_vel_sqr, nb_procs, &
                bloc_size, m_max_pad)

        CALL FFT_PAR_VAR_ETA_PROD_GAUSS_DCL(communicator(2), penal_in_real_space, vv_mesh, &
                grad1_vel_sqr, grad1_vel_sqr_penal, nb_procs, bloc_size, m_max_pad, rr_gauss, time)
        CALL FFT_PAR_VAR_ETA_PROD_GAUSS_DCL(communicator(2), penal_in_real_space, vv_mesh, &
                grad2_vel_sqr, grad2_vel_sqr_penal, nb_procs, bloc_size, m_max_pad, rr_gauss, time)
        CALL FFT_PAR_VAR_ETA_PROD_GAUSS_DCL(communicator(2), penal_in_real_space, vv_mesh, &
                grad3_vel_sqr, grad3_vel_sqr_penal, nb_procs, bloc_size, m_max_pad, rr_gauss, time)

        dissipation = 0.d0
        dissipation_tot = 0.d0
        dissipation_loc = 0.d0
        DO i = 1, SIZE(list_mode)
            mode = list_mode(i)
            IF (mode/=0) THEN
                CYCLE
            ELSE
                index = 0
                DO m = 1, vv_mesh%dom_me
                    DO l = 1, vv_mesh%gauss%l_G
                        index = index + 1
                        !===Compute dissipation
                        dissipation_loc = dissipation_loc + &
                                (grad1_vel_sqr_penal(index, 1, i) + &
                                 grad2_vel_sqr_penal(index, 1, i) + &
                                 grad3_vel_sqr_penal(index, 1, i)) &
                                        * rr_gauss(1, index) * vv_mesh%gauss%rj(l, m)
                    END DO
                END DO
            END IF
        END DO
        dissipation_loc = 2 * ACOS(-1.d0) * dissipation_loc
        CALL MPI_ALLREDUCE(dissipation_loc, dissipation_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, communicator(1), ierr)
        CALL MPI_ALLREDUCE(dissipation_tot, dissipation, 1, MPI_DOUBLE_PRECISION, MPI_SUM, communicator(2), ierr)

        IF (rank==0) THEN
            WRITE(*, *) ' Dissipation ', time, dissipation
            WRITE(998, *) time, dissipation
        END IF

    END SUBROUTINE compute_dissipation

    SUBROUTINE compute_dissipation_ALL(time, vv_mesh, communicator, list_mode, un)
        !over whole domain
        USE def_type_mesh
        USE input_data
        USE boundary
        USE sft_parallele
#include "petsc/finclude/petsc.h"
        USE petsc
        IMPLICIT NONE
        TYPE(mesh_type), INTENT(IN) :: vv_mesh
        INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
        REAL(KIND = 8), DIMENSION(:, :, :), INTENT(IN) :: un
        REAL(KIND = 8), INTENT(IN) :: time
        REAL(KIND = 8), DIMENSION(2, vv_mesh%gauss%l_G * vv_mesh%dom_me) :: rr_gauss
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 6, SIZE(list_mode)) :: grad1_vel
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 6, SIZE(list_mode)) :: grad2_vel
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 6, SIZE(list_mode)) :: grad3_vel
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 2, SIZE(list_mode)) :: grad1_vel_sqr
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 2, SIZE(list_mode)) :: grad2_vel_sqr
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 2, SIZE(list_mode)) :: grad3_vel_sqr
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 2, SIZE(list_mode)) :: grad1_vel_sqr_penal
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 2, SIZE(list_mode)) :: grad2_vel_sqr_penal
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%l_G * vv_mesh%dom_me, 2, SIZE(list_mode)) :: grad3_vel_sqr_penal
        INTEGER, DIMENSION(vv_mesh%gauss%n_w) :: j_loc
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%k_d, vv_mesh%gauss%n_w) :: dw_loc
        REAL(KIND = 8), DIMENSION(vv_mesh%gauss%n_w, 6) :: vel_loc
        REAL(KIND = 8) :: dissipation, dissipation_tot, dissipation_loc
        REAL(KIND = 8) :: ray
        INTEGER :: m, l, i, mode, index, k, TYPE, nb_procs, m_max_pad, bloc_size
        PetscErrorCode                   :: ierr
        MPI_Comm,DIMENSION(2) :: communicator

        DO i = 1, SIZE(list_mode)
            mode = list_mode(i)
            index = 0
            DO m = 1, vv_mesh%dom_me
                j_loc = vv_mesh%jj(:, m)
                DO k = 1, 6
                    vel_loc(:, k) = un(j_loc, k, i)
                END DO
                DO l = 1, vv_mesh%gauss%l_G
                    index = index + 1

                    dw_loc = vv_mesh%gauss%dw(:, :, l, m)

                    !-----------------gauss radius----------------------------------------------
                    rr_gauss(1, index) = SUM(vv_mesh%rr(1, j_loc) * vv_mesh%gauss%ww(:, l))
                    rr_gauss(2, index) = SUM(vv_mesh%rr(2, j_loc) * vv_mesh%gauss%ww(:, l))
                    ray = rr_gauss(1, index)

                    !-----------------Grad u_r on Gauss points------------------------------------
                    grad1_vel(index, 1, i) = SUM(vel_loc(:, 1) * dw_loc(1, :))
                    grad1_vel(index, 2, i) = SUM(vel_loc(:, 2) * dw_loc(1, :))
                    grad1_vel(index, 3, i) = (mode * SUM(vel_loc(:, 2) * vv_mesh%gauss%ww(:, l)) &
                            - SUM(vel_loc(:, 3) * vv_mesh%gauss%ww(:, l))) / ray
                    grad1_vel(index, 4, i) = (-mode * SUM(vel_loc(:, 1) * vv_mesh%gauss%ww(:, l)) &
                            - SUM(vel_loc(:, 4) * vv_mesh%gauss%ww(:, l))) / ray
                    grad1_vel(index, 5, i) = SUM(vel_loc(:, 1) * dw_loc(2, :))
                    grad1_vel(index, 6, i) = SUM(vel_loc(:, 2) * dw_loc(2, :))

                    !-----------------Grad u_th on Gauss points-----------------------------------
                    grad2_vel(index, 1, i) = SUM(vel_loc(:, 3) * dw_loc(1, :))
                    grad2_vel(index, 2, i) = SUM(vel_loc(:, 4) * dw_loc(1, :))
                    grad2_vel(index, 3, i) = (mode * SUM(vel_loc(:, 4) * vv_mesh%gauss%ww(:, l)) &
                            + SUM(vel_loc(:, 1) * vv_mesh%gauss%ww(:, l))) / ray
                    grad2_vel(index, 4, i) = (-mode * SUM(vel_loc(:, 3) * vv_mesh%gauss%ww(:, l)) &
                            + SUM(vel_loc(:, 2) * vv_mesh%gauss%ww(:, l))) / ray
                    grad2_vel(index, 5, i) = SUM(vel_loc(:, 3) * dw_loc(2, :))
                    grad2_vel(index, 6, i) = SUM(vel_loc(:, 4) * dw_loc(2, :))

                    !-----------------Grad u_z on Gauss points------------------------------------
                    grad3_vel(index, 1, i) = SUM(vel_loc(:, 5) * dw_loc(1, :))
                    grad3_vel(index, 2, i) = SUM(vel_loc(:, 6) * dw_loc(1, :))
                    grad3_vel(index, 3, i) = mode * SUM(vel_loc(:, 6) * vv_mesh%gauss%ww(:, l)) / ray
                    grad3_vel(index, 4, i) = -mode * SUM(vel_loc(:, 5) * vv_mesh%gauss%ww(:, l)) / ray
                    grad3_vel(index, 5, i) = SUM(vel_loc(:, 5) * dw_loc(2, :))
                    grad3_vel(index, 6, i) = SUM(vel_loc(:, 6) * dw_loc(2, :))

                END DO
            END DO
        END DO

        CALL MPI_COMM_SIZE(communicator(2), nb_procs, ierr)
        m_max_pad = 3 * SIZE(list_mode) * nb_procs / 2
        bloc_size = SIZE(grad1_vel, 1) / nb_procs + 1

        CALL FFT_PAR_DOT_PROD_DCL(communicator(2), grad1_vel, grad1_vel, grad1_vel_sqr, nb_procs, &
                bloc_size, m_max_pad)
        CALL FFT_PAR_DOT_PROD_DCL(communicator(2), grad2_vel, grad2_vel, grad2_vel_sqr, nb_procs, &
                bloc_size, m_max_pad)
        CALL FFT_PAR_DOT_PROD_DCL(communicator(2), grad3_vel, grad3_vel, grad3_vel_sqr, nb_procs, &
                bloc_size, m_max_pad)

        !CALL FFT_PAR_VAR_ETA_PROD_GAUSS_DCL(communicator(2), penal_in_real_space, vv_mesh, &
        !     grad1_vel_sqr, grad1_vel_sqr_penal, nb_procs, bloc_size, m_max_pad, rr_gauss, time)
        !CALL FFT_PAR_VAR_ETA_PROD_GAUSS_DCL(communicator(2), penal_in_real_space, vv_mesh, &
        !     grad2_vel_sqr, grad2_vel_sqr_penal, nb_procs, bloc_size, m_max_pad, rr_gauss, time)
        !CALL FFT_PAR_VAR_ETA_PROD_GAUSS_DCL(communicator(2), penal_in_real_space, vv_mesh, &
        !     grad3_vel_sqr, grad3_vel_sqr_penal, nb_procs, bloc_size, m_max_pad, rr_gauss, time)
        !TEST LC: no ETA_FFT (we compute dissipation in all domain)
        grad1_vel_sqr_penal = grad1_vel_sqr
        grad2_vel_sqr_penal = grad2_vel_sqr
        grad3_vel_sqr_penal = grad3_vel_sqr

        dissipation = 0.d0
        dissipation_tot = 0.d0
        dissipation_loc = 0.d0
        DO i = 1, SIZE(list_mode)
            mode = list_mode(i)
            IF (mode/=0) THEN
                CYCLE
            ELSE
                index = 0
                DO m = 1, vv_mesh%dom_me
                    DO l = 1, vv_mesh%gauss%l_G
                        index = index + 1
                        !===Compute dissipation
                        dissipation_loc = dissipation_loc + &
                                (grad1_vel_sqr_penal(index, 1, i) + &
                                 grad2_vel_sqr_penal(index, 1, i) + &
                                 grad3_vel_sqr_penal(index, 1, i)) &
                                        * rr_gauss(1, index) * vv_mesh%gauss%rj(l, m)
                    END DO
                END DO
            END IF
        END DO
        dissipation_loc = 2 * ACOS(-1.d0) * dissipation_loc
        CALL MPI_ALLREDUCE(dissipation_loc, dissipation_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, communicator(1), ierr)
        CALL MPI_ALLREDUCE(dissipation_tot, dissipation, 1, MPI_DOUBLE_PRECISION, MPI_SUM, communicator(2), ierr)

        IF (rank==0) THEN
            WRITE(*, *) ' Dissipation ', time, dissipation
            WRITE(999, *) time, dissipation
        END IF

    END SUBROUTINE compute_dissipation_ALL

!! CN 02/2024
!  FUNCTION vv_exact_transform(type,rr,m,t) RESULT(vv)
!    USE user_data
!    IMPLICIT NONE
!    INTEGER     ,                        INTENT(IN)   :: type
!    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
!    INTEGER,                             INTENT(IN)   :: m  !mode 
!    REAL(KIND=8),                        INTENT(IN)   :: t
!    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
!
!    vv = 0.d0
!    IF ((m == 0) .and. (type == 3)) THEN
!       vv(:) = user%solid_vel(1)*rr(1,:)
!    ENDIF
!  END FUNCTION vv_exact_transform
!! CN 02/2024

  !DCQ computes  int (f)/ vol , only in mode zero ( which is the only where int(f) is nonzero)
  FUNCTION mean_int_S(communicator, mesh, list_mode, v) RESULT(norm)
    USE def_type_mesh
    USE chaine_caractere
    USE fem_tn_NS_MHD
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage  
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode  
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v
    INTEGER                                     :: mode_idx
    REAL(KIND=8) :: norm_loc, norm_tot, norm
    INTEGER  :: deb, fin, code!, rank
    INTEGER, DIMENSION(1) :: zero_mode
!v5.0!#include "petsc/finclude/petsc.h"
    MPI_Comm, DIMENSION(2)      :: communicator
    norm = 0.d0
    norm_tot = 0.d0
    !CALL MPI_COMM_RANK(PETSC_COMM_WORLD,rank,code)
    !WRITE(*,*) 'rank, listemode', rank, list_mode(:)
    IF (mesh%me==0) THEN
       norm_loc = 0.d0
       !WRITE(*,*) 'rank, listemode pas de points, normloczero', rank, list_mode(:), norm_loc
    ELSE
       IF (MINVAL(list_mode)==0) THEN !Just mode zero
          zero_mode = MINLOC(list_mode)
          norm_loc = mean_int(mesh, list_mode, zero_mode(1), v)
          !WRITE(*,*) 'rank, listemode, normloc', rank, list_mode(:), norm_loc
       ELSE
          norm_loc=0;
       !WRITE(*,*) 'rank, listemode, normloczero', rank, list_mode(:), norm_loc
       END IF
    END IF
    CALL MPI_ALLREDUCE(norm_loc,norm,1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator(1), code)
    !WRITE(*,*) 'rank, listemode, norm', rank, list_mode(:), norm
    norm =(norm) !no 2pi here
  END FUNCTION mean_int_S

  !DCQ, computes int (f) / vol in mode_idx 
  FUNCTION mean_int(mesh, list_mode, mode_idx, v) RESULT(norm)
    USE def_type_mesh
    USE fem_tn_axi
    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage  
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode  
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) ::  v
    INTEGER                                     :: mode_idx
    REAL(KIND=8)                        :: err1, s1, s2,norm
    INTEGER                             :: k, nn
    err1 = 0.d0
    s2 = 0.d0
    IF (SIZE(v,2)==mesh%np) THEN
       DO nn = 1,SIZE(v,1)
          CALL mean(mesh , (v(nn,:,mode_idx)), s1)
          s2 = s2 + s1
       END DO
       err1 = err1 + s2
       ! CN-AR Tue Jan 13 2009
    ELSE
       DO nn = 1,SIZE(v,2)
          CALL mean(mesh , (v(:,nn,mode_idx)), s1)
          s2 = s2 + s1
       END DO
       ! CN-AR Tue Jan 13 2009
       ! JLG/CN correction du bug CN-AR April 7, 2010
       err1 = err1 + s2
       ! CN-AR Tue Jan 13 2009
    END IF
    norm = err1
  END FUNCTION mean_int

  SUBROUTINE mean (mesh, ff,  t)
    !  < f > /vol  ===>   t 
    USE Gauss_points
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
    REAL(KIND=8),                 INTENT(OUT) :: t
    INTEGER ::  m, l ,i ,ni 
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8) ,DIMENSION(:,:), POINTER       :: rr
    REAL(KIND=8)                                :: ray ,vol 
    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    rr=> mesh%rr
    t = 0
    vol=0
    DO m = 1, me
       DO l = 1, l_G
          !===Compute radius of Gauss point
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
          END DO
          t = t + SUM(ff(jj(:,m)) * ww(:,l))* rj(l,m)*ray
          vol = vol + rj(l,m)*ray
       ENDDO
    ENDDO
    t= t / vol
  END SUBROUTINE mean

END PROGRAM mhd_prog
