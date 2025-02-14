MODULE canopy_files_mod

!-------------------------------------------------------------------------------
! Name:     Files
! Purpose:  Contains FORTRAN units and file names.
! Revised:  15 Jul 2022  Original version.  (P. C. Campbell)
!-------------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER,            PARAMETER :: max_mm     = 1
    INTEGER,            PARAMETER :: iutnml     =  8
    CHARACTER(LEN=256)            :: file_vars ( max_mm )
    CHARACTER(LEN=256)            :: file_out  ( max_mm )
    CHARACTER(LEN=*), PARAMETER   :: file_nml   = 'input/namelist.canopy'

END MODULE canopy_files_mod
