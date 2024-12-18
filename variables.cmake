# Variables to set for Jean-Zay
#set(ENV{CC} "icc")
#set(ENV{CXX} "icpc")
#set(ENV{FC} "ifort")
#set(debug_bounds "-O2 -g -traceback -heap-arrays")
#set(release_bounds "-O3")
#set(SFEMaNS_DIR "/gpfswork/rech/nor/commun/SFEMaNS_GIT/SFEMaNS")
#set(RUN_PRE_PROC "srun")
#set(PROC_CALL "--ntasks=")
#set(RUN_POST_PROC "--hint=nomultithread --job-name=regression_SFEMaNS --time=00:20:00 --partition=visu -A nor@cpu")
##set(RUN_POST_PROC "--hint=nomultithread --job-name=regression_SFEMaNS --time=00:20:00 --qos=qos_cpu-dev -A nor@cpu")

# Variables to set for Whistler
#set(SFEMaNS_DIR "/home/guermond/SFEMaNS/GIT_SFEMaNS/SFEMaNS")
#set(ADDITIONAL_LINKS "-lmetis -lz -L /usr/lib/x86_64-linux-gnu/hdf5/serial")
#set(debug_bounds "-Wall -fimplicit-none -fbounds-check")
#set(release_bounds "-O3")
#set(native_bounds "-march=native -mtune=native -Ofast")
#set(RUN_PRE_PROC "srun")
#set(RUN_PRE_PROC "mpirun")
#set(PROC_CALL "-n ")
#set(RUN_POST_PROC "")

# Variables to set for Ruche
#set(ADDITIONAL_LINKS "-lmetis -lz")
##set(ADDITIONAL_LINKS "-lmetis -lz -L /usr/lib/x86_64-linux-gnu/hdf5/serial")
#set(ENV{CC} "icc")
#set(ENV{CXX} "icpc")
#set(ENV{FC} "ifort")
#set(RUN_PRE_PROC "srun")
#set(PROC_CALL "--ntasks=")
#set(RUN_POST_PROC "--job-name=regression_SFEMaNS --time=00:20:00 --partition=cpu_short")
#set(release_bounds "-O3")
#set(debug_bounds "-O2 -g -traceback -heap-arrays -check bounds -warn all")
#set(SFEMaNS_DIR "/gpfs/users/botezv/SFEMaNS_GIT/SFEMaNS")
#set(SFEMaNS_DIR "/gpfs/users/botezv/SFEMaNS_GIT/SFEMaNS_DEV")

# Variables to set for hemera5
set(ENV{CC} "gcc")
set(ENV{CXX} "gcc")
set(ENV{FC} "gfortran")
set(release_bounds "-O3")
set(debug_bounds "-O2 -g -traceback -heap-arrays")
set(SFEMaNS_DIR "/home/botez18/SFEMaNS")
set(RUN_PRE_PROC "srun")
set(PROC_CALL "--ntasks=")
#set(RUN_PRE_PROC "mpirun")
#set(PROC_CALL "-np ")
set(RUN_POST_PROC "")
#set(ADDITIONAL_LINKS "-L /trinity/shared/pkg/numlib/fftw/3.3.10/gcc/12.2.0/openmpi/4.1.5-ucx/lib/")
set(ADDITIONAL_LINKS "")