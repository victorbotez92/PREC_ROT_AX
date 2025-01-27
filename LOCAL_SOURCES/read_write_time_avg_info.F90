MODULE read_write_time_avg_info
  IMPLICIT NONE
  PUBLIC :: time_average_file

CONTAINS
      
  SUBROUTINE time_average_file(time_avg_counter, avg_kin_ener, URMS, avg_delta, avg_delta_sq)

    IMPLICIT NONE
! VB 23/01/2024
    INTEGER,                             INTENT(IN) :: time_avg_counter
    REAL(KIND=8) :: avg_kin_ener, URMS, avg_delta, avg_delta_sq
 !   INTEGER :: time_avg_counter
! VB 23/01/2024

  !  REAL(KIND=8) :: time, avg_kin_ener, URMS, avg_delta, avg_delta_sq, time_avg_counter_from_file
  ! REAL(KIND=8) :: avg_torque, avg_torque_top, avg_torque_bot
  !  REAL(KIND=8) :: real1, real2, real3, real4, real5, real6, real7


  !  CALL system("tail -n 1 fort.97  > toto1")
  !  CALL system("tail -n 1 fort.99  > toto2")

  !  OPEN(unit=12,file='toto1',form='formatted',status='unknown')
  !  OPEN(unit=13,file='toto2',form='formatted',status='unknown')
  ! OPEN(unit=14,file='toto3',form='formatted',status='unknown')

    !===fort.97
  !  READ(12,*) time, time_avg_counter_from_file, real1, real2, real3, URMS, &
  !            real4, avg_kin_ener, real5
  !  CLOSE(12)

    !===fort.99
  !  READ(13,*) time, time_avg_counter_from_file, real1, avg_delta, avg_delta_sq, &
  !            real2, real3, real4, real5, real6, real7
  !  CLOSE(13)

  ! !===fort.667
  ! READ(14,*) time, time_avg_counter_from_file, avg_torque, avg_torque_top, avg_torque_bot
  ! CLOSE(14)

    !======output
    OPEN(UNIT = 10, FILE = 'time_average_info_new', &
            FORM = 'formatted', STATUS = 'replace')
    WRITE(10,'(A)') '===Time average counter'

    WRITE(10,'(I0)') time_avg_counter

      WRITE(10,'(A)') '===Average Kinetic Energy'
      WRITE(10,'(e15.9)') avg_kin_ener 
      WRITE(10,'(A)') '===URMS'
      WRITE(10,'(e15.9)') URMS
      WRITE(10,'(A)') '===Avg Delta'
      WRITE(10,'(e15.9)') avg_delta
      WRITE(10,'(A)') '===Avg Delta Squared'
      WRITE(10,'(e15.9)') avg_delta_sq
    ! WRITE(10,'(A)') '===Avg torque'
    ! WRITE(10,'(e15.9)') avg_torque
    ! WRITE(10,'(A)') '===Avg torque top'
    ! WRITE(10,'(e15.9)') avg_torque_top
    ! WRITE(10,'(A)') '===Avg torque bot'
    ! WRITE(10,'(e15.9)') avg_torque_bot
      CLOSE(10)
      
  !  CALL system("\rm toto1")
  !  CALL system("\rm toto2")

  END SUBROUTINE time_average_file

END MODULE read_write_time_avg_info
