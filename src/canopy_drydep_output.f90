program canopy_drydep_output
    use canopy_drydep_constants, only: rk=>dp, npts, ninteg
    use canopy_drydep_main, only: GetDryDepExCoeffs

    implicit none

    integer :: i, j
    real(rk) :: vd(npts, ninteg)

    print *, "vd before call:"  ! (uninitialized)
    print "(*(g0))", ((vd(i,j), " ", j = 1, ninteg), new_line("A"), i = 1, npts)
    call GetDryDepExCoeffs(vd)
    print *, "vd after call:"
    print "(*(g0))", ((vd(i,j), " ", j = 1, ninteg), new_line("A"), i = 1, npts)

end program canopy_drydep_output
