MODULE user_data_module
  TYPE personalized_data
     !===I declare my own data here==================================================
     LOGICAL                                 :: if_my_stuff
     !.......Continue here ................................
     REAL(KIND = 8)                          :: disk_radius
     REAL(KIND = 8)                          :: omega_Vol

     LOGICAL                                 :: if_impose_vel
     REAL(KIND=8)                            :: solid_vel
     LOGICAL                                 :: if_nonlin
     REAL(KIND=8)                            :: amp_nl
     LOGICAL                                 :: if_nonlin_m
     INTEGER                                 :: mode_ampli
     REAL(KIND=8)                            :: amp_nl_m


     !Time avg info and restart
     LOGICAL                                 :: read_time_average_info
     INTEGER                                 :: time_avg_freq_restart
     !===End I declare my own data here==============================================
  END TYPE personalized_data
END MODULE user_data_module

MODULE user_data
  USE user_data_module
  IMPLICIT NONE
  PUBLIC :: read_user_data
  TYPE(personalized_data), PUBLIC  :: user
  PRIVATE

CONTAINS

  SUBROUTINE read_user_data(data_file)
    USE my_util
    USE chaine_caractere
    IMPLICIT NONE
    CHARACTER(*),       INTENT(IN) :: data_file
    INTEGER                        :: unit_file=22
    LOGICAL                        :: test

    OPEN(UNIT=unit_file, FILE = data_file, FORM = 'formatted', STATUS = 'unknown')

    !===I add lines that the code SFEMaNS reads in the data file=========================
    CALL find_string(unit_file, '===Should I read my stuff? (true/false)', test)
    IF (test) THEN
       READ (unit_file, *) user%if_my_stuff
    ELSE
       user%if_my_stuff = .FALSE.
    END IF
    !.......Continue here ................................

    !===Nonlinear restart================================================
    CALL find_string(unit_file, '===Is it a nonlinear restart on modes/=0 ?(true/false)', test)
    IF (test) THEN
       READ(unit_file,*) user%if_nonlin
       IF (user%if_nonlin) THEN
          CALL find_string(unit_file, '===Amplitude of amp_nl', test)
          READ(unit_file,*) user%amp_nl
       ELSE
          user%if_nonlin = .FALSE.
          user%amp_nl = 1.d0
       END IF
    ELSE
       user%if_nonlin = .FALSE.
       user%amp_nl = 1.d0
    END IF

    CALL find_string(unit_file, '===Is it a nonlinear restart on one mode?(true/false)', test)
    IF (test) THEN
       READ(unit_file,*) user%if_nonlin_m
       IF (user%if_nonlin_m) THEN
          CALL find_string(unit_file, '===Mode amplified', test)
          READ(unit_file,*) user%mode_ampli
          CALL find_string(unit_file, '===Amplitude of amp_nl_m', test)
          READ(unit_file,*) user%amp_nl_m
       ELSE
          user%if_nonlin_m = .FALSE.
          user%mode_ampli = 0
          user%amp_nl_m = 1.d0
       END IF
    ELSE
       user%if_nonlin_m = .FALSE.
       user%mode_ampli = 0
       user%amp_nl_m = 1.d0
    END IF
    WRITE(*,*) 'user%if_nonlin_m = ', user%if_nonlin_m
    WRITE(*,*) 'user%mode_ampli = ', user%mode_ampli
    WRITE(*,*) 'user%amp_nl_m = ', user%amp_nl_m

    !===Imposed Vel========================================================
    CALL find_string(unit_file, '===Imposed Velocity', test)
    IF (test) THEN
       READ(unit_file,*) user%solid_vel
    ELSE
       user%solid_vel = 0.d0
    END IF

    !===Disk's geometry==================================================
! CN 12/02/2024
        CALL find_string(unit_file, '===Disk Geometry: disk_radius, omega_Vol', test)
        IF (test) THEN
            READ(unit_file, *) user%disk_radius, user%omega_Vol
        ELSE
            user%disk_radius = 0.d0
            user%omega_Vol = 0.d0
        END IF

    !===Read time average info ===============================================
    CALL find_string(unit_file, '===Read time average info? (true/false)', test)
    IF (test) THEN
       READ (unit_file,*) user%read_time_average_info
    ELSE
       user%read_time_average_info=.FALSE.
    END IF


    CALL find_string(unit_file, '===Frequency to write time avg restart file', test)
    IF (test) THEN
       READ (unit_file,*) user%time_avg_freq_restart
    ELSE
       user%time_avg_freq_restart=1d5
    END IF

    !===End I add lines that the code SFEMaNS reads in the data file=====================

    CLOSE(unit_file)
  END SUBROUTINE read_user_data

END MODULE user_data
