!===========================================================================================================
!    Program:    Simplified ACCESS Model
!                Simplified (Dry Deposition only) Atmospheric Chemistry and Canopy Exchange Simulation System
!
!    Developed by : Dr. Beiming Tang based on Dr. Rick Saylor's ACCESS model version 3.1.0
!============================================================================================================

module canopy_drydep_modules
    use canopy_drydep_constants
    implicit none
    
    public MolecDiff, EffHenrysLawCoeff,ReactivityParam, &
           SoilResist, rs_zhang_df, SoilRbg,  &
           mdiffh2o, RelativeHumidity,gpla,  &
           rbl, rcl, rml, CalcWeightedProfiles
         
      

contains

!===============================================
!PART 1. Functions directly called by DryDep.f90
!
!
!===============================================

!=========================================================================
!function RelativeHumidity - calculate Relative Humidity from the supplied 
!                            specific humidity, temperature and pressure
!=========================================================================
function RelativeHumidity(tki,pmbi,qhi)
    real(kind = rk), intent(in) :: tki
    real(kind = rk), intent(in) :: pmbi
    real(kind = rk), intent(in) :: qhi
    real(kind = rk)             :: RelativeHumidity
    real(kind = rk)             :: e
    real(kind = rk)             :: es
    real(kind = rk)             :: qhd
    real(kind = rk)             :: pkpa
    real(kind = rk)             :: rhi
    real(kind = rk),parameter   :: rhmin = 0.1
    real(kind = rk),parameter   :: rhmax = 99.0

    qhd  = 0.001*qhi
    pkpa = 0.1*pmbi
    e    = pkpa*qhd/(0.622+qhd)
    es   = esat(tki)
    rhi  = max(rhmin,min(rhmax,100.0*e/es))   !bound RH to (rhmin, rhmax)
    RelativeHumidity = rhi

    return
end function RelativeHumidity


!==========================================================================
!function MolecDiff - calculate molecular diffusivities (cm^2/s) at a given
!                     temperature and pressure
!==========================================================================
function MolecDiff(ispec, tkx, pmbx)
    integer(kind=ik), intent(in) :: ispec     !dummy id for species
    real(kind=rk), intent(in)    :: tkx
    real(kind=rk), intent(in)    :: pmbx
    real(kind=rk)                :: MolecDiff
    real(kind=rk),dimension(ntotal)  ::mdiffstp
    
    call SetMolecDiffSTP(mdiffstp)
    MolecDiff = mdiffstp(ispec)*(1013.25_rk/pmbx)*((tkx/298.15_rk)**1.81)

    return
end function MolecDiff


!=========================================================================
!function mdiffh2o - calculate molecular diffusivity of water vapor in air
!
!source: Tracy (1980)
!=========================================================================
function mdiffh2o(tki,pmbi)
    real(kind = rk), intent(in) :: tki
    real(kind = rk), intent(in) :: pmbi
    real(kind = rk)             :: mdiffh2o

    mdiffh2o = 0.226_rk*((tki/273.15_rk)**1.81_rk)*(1000.0_rk/pmbi)

    return
end function mdiffh2o


!==========================================================================
!function ReactivityParam - calculate reactivity parameters (dimensionless)
!==========================================================================
function ReactivityParam(ispec)
    integer(kind=ik), intent(in)    :: ispec   !dummy id for species
    real(kind=rk)                   :: ReactivityParam
    real(kind=rk),dimension(ntotal) :: f0

    call SetReactivityParams(f0)
    ReactivityParam = f0(ispec)

    return
end function ReactivityParam


!========================================================================
!function EffHenrysLawCoeff - calculate Henry's Law coefficient (M/atm) 
!                             have not include temperature dependence yet
!========================================================================
function EffHenrysLawCoeff(ispec)
    integer(kind=ik), intent(in)    :: ispec  !dummy id for species
    real(kind=rk)                   :: EffHenrysLawCoeff
    real(kind=rk),dimension(ntotal) :: hstar     

    call SetEffHenrysLawCoeffs(hstar)
    EffHenrysLawCoeff = hstar(ispec)

    return
end function EffHenrysLawCoeff


!===============================================================================================================
!function SoilResist - calculate the resistance to diffusion of a species from the free warer surface in the soil
!                      to the soil-atmosphere interface. Rsoil 
!Source - Sakagichi & Zeng (2009)
!================================================================================================================
function SoilResist(mdiffl)
    real(kind=rk), intent(in)  :: mdiffl       !molecular diffusivity of species in air (cm^2/s)
    real(kind=rk)              :: SoilResist   !Soil resistance (s/cm)
    real(kind=rk)              :: xe           !temporary variable
    real(kind=rk)              :: ldry         !diffusion distance through the soil (cm)
    real(kind=rk)              :: mdiffp       !effective diffusivity of species through the soil (cm^2/s)

    xe = (1.0-(stheta/sattheta))**5.0
    ldry = dsoil*(exp(xe)-1.0)/1.7183
    ldry = max(0.0_rk,ldry)

    mdiffp = mdiffl*sattheta*sattheta*(1.0-(rtheta/sattheta))**(2.0+3.0/sbcoef)
    SoilResist = ldry/mdiffp

    return
end function SoilResist


!=====================================================================================
!function SoilRbg - calculate the boundary layer resistance at the ground surface. Rbg
!
!Source - Schuepp (1977) 
!=====================================================================================
function SoilRbg(ubarg)
    real(kind=rk), intent(in)  :: ubarg           !mean wind speed in the 1st model layer (cm/s)
    real(kind=rk)              :: SoilRbg         !Boundary layer resistance at ground surface (s/cm)
    real(kind=rk), parameter   :: rbgmax = 1.67   !maximum ground surface boundary layer resistance (s/cm)
    real(kind=rk)              :: rbg             !temporary variable for boundary layer resistance (s/cm)

    rbg = 11.534/(0.14*ubarg)                     !assume Sc=0.7,del0/zl = 0.02 and ustar = 0.14*ubar (Weber,1999)
    SoilRbg = min(rbgmax,rbg)

    return
end function SoilRbg 


!===================================================================
!function rbl - calculate leaf boundary resistance for trace species
!
!Source - rb formulation from Wu et al., (2003)
!       - ustar = 0.14*ubar from Weber (1999)
!===================================================================
function rbl(mdiffl, ubari)
    real(kind=rk)               :: rbl            !leaf boundary layer resistance (s/cm)
    real(kind=rk), intent(in)   :: mdiffl         !molecular diffusivity of species in air (cm^2/s)
    real(kind=rk), intent(in)   :: ubari          !mean wind speed at layer i (cm/s)
    
    rbl = 10.53_rk/((mdiffl**0.666667)*ubari)

    return
end function rbl


!===============================================================
!function rcl - calculate cuticular resistance for trace species
!
!Source - Wesely (1989)
!===============================================================
function rcl(hstarl,f01)
    real(kind=rk), intent(in)   :: hstarl         ! effective Henry's law coefficient (M/atm)
    real(kind=rk), intent(in)   :: f01            ! reactivity parameter (0-1)
    real(kind=rk)               :: rcl            ! cuticular resistance (s/cm)
    real(kind=rk), parameter    :: rcref=20.0     ! rc for ozone (s/cm) for deciduous forest

    rcl = rcref/((hstarl*1.0D-05)+f01)

    return
end function rcl


!===============================================================
!function rml - calculate mesophyll resistance for trace species
!
!Source - Wesely (1989)
!===============================================================
function rml(hstarl,f01)
    real(kind=rk), intent(in)   :: hstarl          ! effective Henry's law coefficient (M/atm)
    real(kind=rk), intent(in)   :: f01             ! reactivity parameter (0-1)
    real(kind=rk)               :: rml             ! mesophyll resistance (s/cm)

    rml = 1.0_rk/((hstarl/3000.0_rk)+100.0_rk*f01)

    return
end function rml


!==================================================================================
!function gpla - calculate stomatal compensation point for ACCESS-ORG trace species
!
!Source - ???
!==================================================================================
function gpla(ilev,lsp)
    integer(kind=ik), intent(in) :: ilev           !ilev is layer
    integer(kind=ik), intent(in) :: lsp            !lsp is species
    real(kind=rk)                :: gpla

    gpla = 0.0_rk                                  !dummy stub for now

    return
end function gpla


!======================================================================
!function rs_zhang_df - calcualte stomatal resistance for trace species
!
!Source - Zhang et al., (2002 & 2003)
!       - 'rsmin' by Wesely et al (1989)
!       - 'bvpd' by Wolfe and Thornton (2011)
!======================================================================
function rs_zhang_df(mdiffl, tki, pmbi, ppfdi, srad, relhumi)
    real(kind = rk), intent(in) :: mdiffl            !molecular diffusivity of trace species in air (cm^2/s)
    real(kind = rk), intent(in) :: tki               !temperature
    real(kind = rk), intent(in) :: pmbi              !air pressure (mb)
    real(kind = rk), intent(in) :: ppfdi             !photosynthetic photon flux (umol/m^2-s)
    real(kind = rk), intent(in) :: srad              !solar irradiation (W/m^2)
    real(kind = rk), intent(in) :: relhumi           !relative humidity (%)
    real(kind = rk)             :: rs_zhang_df       !stomatal resistance (s/cm)
    real(kind = rk), parameter  :: rsmin =1.0        !minimum leaf stomatal resistance (s/cm) for deciduous forest
    real(kind = rk), parameter  :: rsmax =10000.     !maximum leaf stomatal resistance (s/cm) (stoma are closed)
    real(kind = rk), parameter  :: brsp = 196.5      !empirical constant (umol/m^2-s) for deciduous forest
    real(kind = rk), parameter  :: tmin = 0.0        !temperature correction parameter-deciduous forest
    real(kind = rk), parameter  :: tmax = 45.0       !temperature correction parameter-deciduous forest
    real(kind = rk), parameter  :: topt = 27.0       !temperature correction parameter-deciduous forest
    real(kind = rk), parameter  :: bvpd = 0.10       !empirical constant for VPD correction-deciduous forest
    real(kind = rk), parameter  :: phic1 = -1.9      !empirical constant for water stress correction-deciduous forest
    real(kind = rk), parameter  :: phic2 = -2.5      !empirical constant for water stress correction-deciduous forest
    real(kind = rk)             :: cft
    real(kind = rk)             :: cfvpd
    real(kind = rk)             :: cfphi
    real(kind = rk)             :: tcel
    real(kind = rk)             :: ft1
    real(kind = rk)             :: ft2
    real(kind = rk)             :: et
    real(kind = rk)             :: vpd
    real(kind = rk)             :: phi

    !temperature correction
    tcel = tki - 273.15
    et   = (tmax-topt)/(topt-tmin)
    ft1  =(tcel-tmin)/(topt-tmin)
    ft2  = (tmax-tcel)/(tmax-topt)
    cft  = ft1*(ft2**et)

    !water vapor pressure defit correction
    vpd  = esat(tki)*(1.0_rk - (relhumi/100.0_rk))
    cfvpd= 1.0_rk - bvpd*vpd

    !water stress correction
    phi  = -0.72_rk - 0.0013_rk*srad
    cfphi= (phi-phic2)/(phic1-phic2)

    if (ppfdi >0.0) then
        rs_zhang_df = rsmin*(1.0_rk+brsp/ppfdi)*mdiffh2o(tki,pmbi)/(mdiffl*cft*cfvpd*cfphi)
    else
        rs_zhang_df = rsmax                         !nighttime, stoma are closed
    endif

    return
end function rs_zhang_df

!==================================================================
!CalWeightedProfiles - calculate sun/shade weighted canopy profiles
!
!==================================================================

subroutine CalcWeightedProfiles(rs_wgt,fsun,fshd,gs_sun,gs_shd,pmb,tk)    !Beiming modify this code from CalcWeightedProfiles from CanopyPhysics.f90
    integer(kind=ik)                             :: i                     !i is layer
    real(kind=rk),dimension(npts),intent(out)    :: rs_wgt
    real(kind=rk),dimension(npts),intent(in)     :: fsun
    real(kind=rk),dimension(npts)                :: rs_sun
    real(kind=rk),dimension(npts),intent(in)     :: fshd
    real(kind=rk),dimension(npts)                :: rs_shd
    real(kind=rk),dimension(npts),intent(in)     :: gs_sun
    real(kind=rk),dimension(npts),intent(in)     :: gs_shd    
    real(kind=rk),dimension(npts),intent(in)     :: pmb 
    real(kind=rk),dimension(npts),intent(in)     :: tk    
   
    do i = npts,1,-1
!        ppfd_wgt(i) = fun(i) *ppfd_sun(i) + fshd(i)*ppfd_shd(i)
!        nir_wgt(i) = fun(i) *nir_sun(i) + fshd(i)*nir_shd(i)
!        rt_wgt(i) = fun(i) *rt_sun(i) + fshd(i)*rt_shd(i)
!        rabs_wgt(i) = fun(i) *rabs_sun(i) + fshd(i)*rabs_shd(i)
!        tl_wgt(i) = fun(i) *tl_sun(i) + fshd(i)*tl_shd(i)
!        gs_wgt(i) = fun(i) *gs_sun(i) + fshd(i)*gs_shd(i)
        
        if (i <= ncnpy) then            !ncnpy is # of vertical layers
            rs_sun(i) = gtor(gs_sun(i),pmb(i),tk(i))
            rs_shd(i) = gtor(gs_shd(i),pmb(i),tk(i))
        endif
 !       print *,'rs_sun',rs_sun
 !       print *,'rs_shd',rs_shd
 !       print *,'fsun',fsun
 !       print *, 'fshd',fshd
        rs_wgt(i) = fsun(i)*rs_sun(i) + fshd(i)*rs_shd(i)
!        anet_wft(i) = fsun(i)*anet_sun(i) + fshd(i)*anet_shd(i)
    end do

!    do i = npts,1,-1
!        if (ppfd_sun(ncnpy+1) >0.0) then
!            fj(i) = pppfd_sun(i)/ppfd_sun(ncnpy+1)
!        else
!            fj(i) = 0.0
!        endif
!        fjout(i,nt) = fj(i)
!    enddo

    return
end subroutine CalcWeightedProfiles
 
!=============================================================
!PART 2. 2nd level functions called by functions in DryDep.f90
!
!
!=============================================================

!==============================================================================================
!function esat - calculate the saturation vaport pressure (kPa) of water at a given temperature
!
!source - Rogers et al., (1989)
!notes - valid over -30 C <= T <= 35 C
!==============================================================================================
function esat(tki)
    real(kind=rk),intent(in) :: tki    !temperature,unit = K
    real(kind=rk)            :: esat   !satuaration vapor pressure,unit = kPa
    real(kind=rk)            :: tc     !temperature, unit = C

    tc = tki - 273.15
    esat = 0.6112*exp(17.67*tc/(tc+243.5))

    return
end function esat

!=====================================================================================
!function mdiffstp - set molecular diffusivity data (cm^2/s) for all species at
!                           0 deg C and 1 atm 
!=====================================================================================
subroutine SetMolecDiffSTP(mdiffstp)
    real(kind = rk),dimension(ntotal),intent(out)  :: mdiffstp       !molecular diffusivities of species in air at 0degC and 1 atm [cm^2/s]
!    real(kind = dp), parameter         :: mdiffstp_default = 0.100   !default value of mdiffstp (cm^2/s) with no reliable data
!    integer(kind=i4)                   :: l                          !l is species
    
!    do l=1,ntotal      !ntotal is total species in ACCESS based on chosen chemical mechanim, including transported species
!        mdiffstp(l) = mdiffstp_default
!    end do

!Insert MolecDiffSTP Coefficients for RACM2_plus mechanism
!species in array are:
!         = (/NO,   NO2,    O3,  HONO,  HNO4,   HNO3,   N2O5,     CO,  H2O2,  CH4,
!            MO2,   OP1,   MOH,   NO3,   O3P,    O1D,     HO,    HO2,  ORA1,  HAC,
!            PAA, DHMOB, HPALD,  ISHP, IEPOX, PROPNN, ISOPNB, ISOPND, MACRN, MVKN,
!           ISNP /)  
    mdiffstp =(/0.1802, 0.1361, 0.1444, 0.1349, 0.1041, 0.1041, 0.0808, 0.1807, 0.1300, 0.1952,  &
                0.1297, 0.1200, 0.1297, 0.1153, 0.2773, 0.2773, 0.2543, 0.2000, 0.1340, 0.1060,  &
                0.1040, 0.0892, 0.0845, 0.0837, 0.0837, 0.0834, 0.0750, 0.0750, 0.0745, 0.0745,  &
                0.0712 /)
!species (NO2, O3)
!    mdiffstp =(/0.1361, 0.1444/) 

    return
end subroutine SetMolecDiffSTP

!========================================================
!function f0 - set reactivity parameters for all species
!
!source - Wesely et al., (1989) and Nguyen et al., (2015)
!========================================================
subroutine SetReactivityParams(f0)
!    integer(kind=i4)                              :: l               !l is species
    real(kind = rk),dimension(ntotal),intent(out) :: f0
!    real(kind=dp),parameter                       :: f0_default = 0.0
    
!    do l=1,ntotal     !ntotal is total species in ACCESS based on chosen chemical mechanim, including transported species
!        f0(l) = f0_default
!    end do

!Insert ReactivityParams Coefficients for RACM2_plus mechanism
!species in array are:
!         = (/NO,   NO2,    O3,  HONO,  HNO4,   HNO3,   N2O5,     CO,  H2O2,  CH4,
!            MO2,   OP1,   MOH,   NO3,   O3P,    O1D,     HO,    HO2,  ORA1,  HAC,
!            PAA, DHMOB, HPALD,  ISHP, IEPOX, PROPNN, ISOPNB, ISOPND, MACRN, MVKN,
!           ISNP /)  
    f0 = (/0.0, 0.1, 1.0, 0.1, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,  &
           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.0, 0.0,  &
           0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  &
           0.0 /)

!species (NO2, O3)
!    f0 = (/0.1, 1.0 /)

    return
end subroutine SetReactivityParams 

!====================================================================
!function hstar - set henry's law coefficient for all species (M/atm)
!
!source - Nguyen et al., (2015)
!====================================================================
subroutine SetEffHenrysLawCoeffs(hstar)
!    integer(kind=i4)                            :: l                  !l is species
    real(kind=rk),dimension(ntotal),intent(out) :: hstar
!    real(kind=dp),parameter                     :: hstar_default = 1.000
    
!    do l = 1,ntotal    !ntotal is total species in ACCESS based on chosen chemical mechanim, including transported species
!        hstar(l) = hstar_default
!    End do

!Insert EffHenryLaw Coefficients for RACM2_plus mechanism
!species in array are:
!         = (/NO,   NO2,    O3, HONO,  HNO4,   HNO3,   N2O5,     CO,  H2O2,  CH4,
!            MO2,   OP1,   MOH,  NO3,   O3P,    O1D,     HO,    HO2,  ORA1,  HAC,
!            PAA, DHMOB, HPALD, ISHP, IEPOX, PROPNN, ISOPNB, ISOPND, MACRN, MVKN,
!           ISNP /)    

!    hstar = (/0.0019, 0.0120, 0.0103, 2.6D+05, 1.0D+07, 3.2D+13, 1.0D+14, 9.8D-04, 8.4D+04, 1.4D-03,  &
!             6.9D+02, 3.0D+02, 2.2D+02, 0.0380, 0.0380, 0.0380, 3.9D+01, 6.9D+02, 5.6D+03, 2.0D+03,  & 
!             5.2D+02, 2.0D+03, 4.0D+04, 7.0D+07, 7.0D+07, 1.0D+04, 5.0D+03, 5.0D+03, 6.0D+03, 6.0D+03,  &
!             5.0D+03 /)

    hstar = (/1.9D-03, 1.2D-02,1.03D-02, 2.6D+05, 1.0D+07, 3.2D+13, 1.0D+14, 9.8D-04, 8.4D+04, 1.4D-03,  &
              6.9D+02, 3.0D+02, 2.2D+02, 3.8D-02, 3.8D-02, 3.8D-02, 3.9D+01, 6.9D+02, 5.6D+03, 2.0D+03,  &
              5.2D+02, 2.0D+03, 4.0D+04, 7.0D+07, 7.0D+07, 1.0D+04, 5.0D+03, 5.0D+03, 6.0D+03, 6.0D+03,  &
              5.0D+03 /)

!species (NO2, O3)
!    hstar = (/0.0120, 0.0103/)
    
    return
end subroutine SetEffHenrysLawCoeffs

!==================================================================
!function gtor - convert conductance (mol/m^2-s) to resistance (s/m)
!
!notes - assumes non-zero value for gz!
!==================================================================
function gtor(gz, pmbi,tki)
    real(kind=rk), intent(in)     :: gz                 !conductance, mol/m^2-s
    real(kind=rk), intent(in)     :: pmbi               !air pressure, mb
    real(kind=rk), intent(in)     :: tki                !temperature, K
    real(kind=rk)                 :: gtor               !resistance, s/s
    real(kind=rk)                 :: gms                !conductance, m/s
    real(kind=rk), parameter      :: rgas=8.205D-05     !ideal gas constant, m^3-atm/K-mol

    gms = gz*(rgas*tki)/(pmbi/1013.)
    gtor = 1.0/gms

    return
end function gtor

end module canopy_drydep_modules


