set(CMAKE_Fortran_COMPILER "/trinity/shared/pkg/compiler/gcc/12.2.0/bin/gfortran")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "12.2.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_COMPILER_FRONTEND_VARIANT "GNU")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/trinity/shared/pkg/compiler/gcc/12.2.0/bin/gcc-ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "/trinity/shared/pkg/compiler/gcc/12.2.0/bin/gcc-ranlib")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/trinity/shared/pkg/compiler/gcc/12.2.0/lib/gcc/x86_64-pc-linux-gnu/12/finclude;/trinity/shared/pkg/numlib/fftw/3.3.10/gcc/12.2.0/openmpi/4.1.5-ucx/include;/trinity/shared/pkg/devel/metis/5.1.0/include;/trinity/shared/pkg/analysis/petsc/3.22.2/include;/trinity/shared/pkg/devel/valgrind/3.14.0/include;/trinity/shared/pkg/filelib/hdf5/1.14.0/gcc/12.2.0/include;/trinity/shared/pkg/devel/zlib/1.2.11/include;/trinity/shared/pkg/compiler/gcc/12.2.0/include;/trinity/shared/pkg/compiler/gcc/12.2.0/lib/gcc/x86_64-pc-linux-gnu/12/include;/usr/local/include;/trinity/shared/pkg/compiler/gcc/12.2.0/lib/gcc/x86_64-pc-linux-gnu/12/include-fixed;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;m;gcc_s;gcc;quadmath;m;gcc_s;gcc;c;gcc_s;gcc")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/trinity/shared/pkg/compiler/gcc/12.2.0/lib64;/trinity/shared/pkg/compiler/gcc/12.2.0/lib/gcc/x86_64-pc-linux-gnu/12;/lib64;/usr/lib64;/trinity/shared/pkg/devel/metis/5.1.0/lib;/trinity/shared/pkg/numlib/fftw/3.3.10/gcc/12.2.0/openmpi/4.1.5-ucx/lib;/trinity/shared/pkg/devel/valgrind/3.14.0/lib;/trinity/shared/pkg/compiler/gcc/12.2.0/lib;/trinity/shared/pkg/filelib/hdf5/1.14.0/gcc/12.2.0/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
