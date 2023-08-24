!=======================================================================================================
!    Program: Simplied ACCESS model
!             Simplied (Dry deposition only) Atmospheric Chemistry and Canopy Exchange Simulation System
!
!    developed by : Dr. Beiming Tang based on Dr. Rick D. Saylor's ACCESS model version 3.1.0
!=======================================================================================================

module canopy_drydep_main
    use canopy_drydep_constants
    use canopy_drydep_modules

    implicit none

    public  GetDryDepExCoeffs, GetSoilDryDepExCoeffs

contains

!============================================================================
!PART 1. CALCULATE LEAF-SCALE DRY DEPOSITION VELOCITIES & COMPENSATION POINTS
!============================================================================

subroutine GetDryDepExCoeffs(vd)
    integer(kind=ik)                                       :: i                       !i is layer
    integer(kind=ik)                                       :: l                       !l is species
    real(kind = rk)                                        :: mdiffl                  !molecular diffusivity of species l in air (cm^2/s)
    real(kind = rk)                                        :: relhumi                 !relative humidity (%)
    real(kind = rk)                                        :: hstarl                  !effective Henry's Law coefficient (M/atm)
    real(kind=rk)                                          :: f01                     !reactivity parameter (0-1) !Wesely's reactivity parameter (dimensionless)
    real(kind = rk)                                        :: srad                    !Solar radiation at canopy top (W/m^2)
    character(len=10)                                      :: rs_select               !Selection of stomatal resistance algorithm
    real(kind=rk)                                          :: hc                      !hc is canopy height (cm)
    real(kind=rk),dimension(npts)                          :: ppfd                    !vertical profile of photosynthetic photon flux (umol/m^2-s) at current simulation time
    real(kind=rk),dimension(npts)                          :: lai                     !layer leaf area index of canopy (cm^2*leaf/cm^2)
    real(kind=rk),dimension(npts)                          :: tk                      !vertical profile of temperature at current simulation time (K)
    real(kind=rk),dimension(npts)                          :: pmb                     !vertical profile of pressure at current simulation time (mb)
    real(kind=rk),dimension(npts)                          :: qh                      !vertical profile of specific humidity at current simulation time (g/kg)
    real(kind=rk),dimension(npts)                          :: ubar                    !vertical profile of mean wind speed at current simulation time (cm/s)
    real(kind = rk),dimension(npts,ninteg)                 :: rb                      !leaf boundary resistance (s/cm)
    real(kind = rk),dimension(npts,ninteg)                 :: rc                      !cuticular resistance (s/cm)
    real(kind = rk),dimension(npts,ninteg)                 :: rm                      !mesophyll resistance (s/cm)
    real(kind = rk),dimension(npts,ninteg)                 :: rs                      !stomatal resistance (s/cm)
    real(kind = rk),dimension(npts)                        :: rs_wgt                  !sun/shade weighted leaf stomatal resistance (s/cm)
    real(kind=rk),dimension(npts,ninteg), intent(inout)    :: vd                      !dry deposition exchange coefficient (cm/s)
    real(kind = rk),dimension(npts,ninteg)                 :: gp                      !dry deposition compensation point (mole/cm^3)
    real(kind = rk)                                        :: rnum
    real(kind = rk)                                        :: rden
    real(kind = rk)                                        :: rlx
    real(kind = rk)                                        :: vdlx
    real(kind=rk),dimension(npts)                          :: fshd                    !shaded fraction of canopy
    real(kind=rk),dimension(npts)                          :: fsun                    !sunlit fraction of canopy
    real(kind=rk),dimension(npts)                          :: gs_sun                  !leaf stomatal conductance in sublit fraction (mol/m^2-s)
    real(kind=rk),dimension(npts)                          :: gs_shd 
    
    !nighttime 1am, get input data from 5m to 40m, every 5m interval
!    pmb    = (/1012.66     , 1012.07     , 1011.48     , 1010.891    , 1010.292    , 1009.702    , 1009.111    , 1008.522    /)              !calc use barometric formula
!    lai    = (/0.444595E+01, 0.395672E+01, 0.331488E+01, 0.211422E+01, 0.148396E+01, 0.392690E+00, 0.000000E+00, 0.000000E+00/)              !data from CANACC
!    tk     = (/0.293843E+03, 0.293680E+03, 0.293552E+03, 0.293149E+03, 0.292633E+03, 0.292133E+03, 0.292153E+03, 0.292183E+03/)              !data from CANACC
!    ppfd   = (/0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00/)              !data from CANACC
!    qh     = (/5.2         , 5.2         , 5.2         , 5.2         , 5.2         , 5.2         , 5.2         , 5.2         /)              !calc use RH=30% T=20C
!    ubar   = (/0.540616E+01, 0.688076E+01, 0.898461E+01, 0.121621E+02, 0.174182E+02, 0.278164E+02, 0.861199E+02, 0.861199E+02/)              !data from CANACC
!    fsun   = (/0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00/)              !data from CANACC
!    fshd   = (/0.100000E+01, 0.100000E+01, 0.100000E+01, 0.100000E+01, 0.100000E+01, 0.100000E+01, 0.100000E+01, 0.100000E+01/)              !data from CANACC
!    gs_sun = (/0.200000E-01, 0.200000E-01, 0.200000E-01, 0.200000E-01, 0.200000E-01, 0.200000E-01, 0.200000E-01, 0.000000E+00/)              !data from CANACC
!    gs_shd = (/0.200000E-01, 0.200000E-01, 0.200000E-01, 0.200000E-01, 0.200000E-01, 0.200000E-01, 0.200000E-01, 0.000000E+00/)              !data from CANACC

    !daytime 11am, get input data from 5m to 40m, every 5m interval
    pmb    = (/1012.665, 1012.09     , 1011.515    , 1010.945    , 1010.37     , 1009.783    , 1009.165    , 1008.567    /)              !calc use barometric formula
    lai    = (/4.45E+00, 0.395672E+01, 0.331488E+01, 0.211422E+01, 0.148396E+01, 0.393690E+00, 0.000000E+00, 0.000000E+00/)              !data from CANACC
    tk     = (/2.96E+02, 0.298247E+03, 0.299576E+03, 0.300502E+03, 0.300629E+03, 0.299304E+03, 0.296860E+03, 0.295961E+03/)              !data from CANACC
    ppfd   = (/1.79E+02, 0.230477E+03, 0.317002E+03, 0.558829E+03, 0.743638E+03, 0.119511E+04, 0.140834E+04, 0.225315E+04/)              !data from CANACC
    qh     = (/5.2     , 5.2         , 5.2         , 5.2         , 5.2         , 5.2         , 5.2         , 5.2         /)              !calc use RH=30% T=20C
    ubar   = (/1.77E+01, 0.225845E+02, 0.294898E+02, 0.399192E+02, 0.571712E+02, 0.913008E+02, 0.282668E+03, 0.282688E+03/)              !data from CANACC
    fsun   = (/8.20E-02, 0.108003E+00, 0.154963E+00, 0.304459E+00, 0.434003E+00, 0.801360E+00, 0.100000E+01, 0.100000E+01/)              !data from CANACC
    fshd   = (/9.18E-01, 0.891997E+00, 0.845037E+00, 0.695541E+00, 0.565997E+00, 0.198640E+00, 0.000000E+00, 0.000000E+00/)              !data from CANACC
    gs_sun = (/1.76E+00, 0.886976E+00, 0.509293E+00, 0.379802E+00, 0.328706E+00, 0.310937E+00, 0.314569E+00, 0.000000E+00/)              !data from CANACC
    gs_shd = (/2.88E-01, 0.192159E+00, 0.156742E+00, 0.182532E+00, 0.188411E+00, 0.198841E+00, 0.195128E+00, 0.000000E+00/)              !data from CANACC

    rs_select = 'zhang_df'                                             !To Do: make this selection via input file
    srad = ppfd(int(hc)+1)/4.57_rk                                      !cannopy top PPFD converted from (umol/m^2-s) to (W/m^2)
    do l = 1,ninteg                                                     !ninteg is integrated # of species in ACCESS model
        do i = 1,npts                                                   !npts is # of vertical layers
            if (lai(i) > 0.0) then                                      !within canopy
                mdiffl  = MolecDiff(l,tk(i),pmb(i))                     !calculate molecular diffusivity (cm^2/s)
                relhumi = RelativeHumidity(tk(i),pmb(i),qh(i))          !calculate relative humidity at layer i
                rb(i,l) = rbl(mdiffl, ubar(i))                          !leaf boundary layer resistance (s/cm)
                
                hstarl  = EffHenrysLawCoeff(l)
                f01     = ReactivityParam(l)
                rc(i,l) = rcl(hstarl, f01)                              !leaf cuticular resistance (s/cm)
                rm(i,l) = rml(hstarl, f01)                              !leaf mesophyll resistance (s/cm)
                select case (rs_select)
                    case('zhang_df')                                    !Zhang et al., (2002,2003) with hardcodes deciduous forest parameters
                        rs(i,l) = rs_zhang_df(mdiffl,tk(i),pmb(i),ppfd(i),srad,relhumi)
                    case('gs_medlyn')                                   !Medlyn et al., (2011) 
                        Call CalcWeightedProfiles(rs_wgt,fsun,fshd,gs_sun,gs_shd,pmb,tk)   !this is to get rs_wgt list
                        rs(i,l) = (mdiffh2o(tk(i),pmb(i))/mdiffl)*rs_wgt(i)* 0.01   !convert rs_wgt from (s/m) to (s/cm)
                    case default                                        !Zhang et la., (2002,2003) wiht hardcodes deciduous forest parameters
                        rs(i,l) = rs_zhang_df(mdiffl,tk(i),pmb(i),ppfd(i),srad,relhumi)
                end select
                rnum = rc(i,l) * (rs(i,l) + rm(i,l))
                rden = rc(i,l) + 2.0 * (rs(i,l) + rm(i,l))
                rlx   = rb(i,l) + (rnum/rden)
                vdlx  = 1.0/rlx
                vd(i,l) = vdlx                                           !calcualte deposition velocity (cm/s)
                gp(i,l) = gpla(i,l)                                      !compensation point concentration (mole/cm^3)
           
            else                                                         !out of canopy
                rb(i,l) = 0.0_rk
                rc(i,l) = 0.0_rk
                rm(i,l) = 0.0_rk
                rs(i,l) = 0.0_rk
                vd(i,l) = 0.0_rk
                gp(i,l) = 0.0_rk

            end if
        end do
    end do

    return
end subroutine GetDryDepExCoeffs

!==================================================================
! PART 2. CALCULATE DRY DEPOSITION VELOCITIES TO THE GROUND SURFACE
!==================================================================

subroutine GetSoilDryDepExCoeffs()
    integer(kind = ik)                :: l                       !l is species
    real(kind = rk)                   :: mdiffl                  !molecular diffusivity (cm^2/s)
    real(kind = rk)                   :: tsoilk                  !soil/litter temperature (K)
    real(kind = rk),dimension(ninteg) :: rsoill                  !resistance to diffusion thru soil pore space for chemical species (s/cm)
    real(kind = rk)                   :: rbg                     !ground boundary layer resistance (s/cm)
    real(kind = rk),dimension(ninteg) :: vs                      !soil exchange coefficients (cm/s)
    real(kind = rk)                   :: ubarg
    real(kind=rk),dimension(npts)     :: pmb
    real(kind=rk),dimension(npts)     :: ubar

    pmb    = (/1.01266     , 1.01207     , 1.01148     , 1.010891    , 1.010292    , 1.009702    , 1.009111    , 1.008522    /)              !calc use barometric formula
    ubar   = (/0.540616E+01, 0.688076E+01, 0.898461E+01, 0.121621E+02, 0.174182E+02, 0.278164E+02, 0.861199E+02, 0.861199E+02/)              !data from CANACC

    do l = 1,ninteg    !ninteg is integrated # of species in ACCESS
        mdiffl = MolecDiff(l, tsoilk, pmb(1))
        rsoill(l) = SoilResist(mdiffl)

        ubarg = ubar(1)                                          !get ground level mean wind speed
        rbg = SoilRbg(ubarg)                                     !Rbg(ground boundary layer resistance, s/cm)
                                                                 !Rbg is invariant to species not layers
        vs(l) = 1.0/(rbg+rsoill(l))                              !deposition velocity to ground surface (cm/s)
    end do

    return
end subroutine GetSoilDryDepExCoeffs

end module canopy_drydep_main
