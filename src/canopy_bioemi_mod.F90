module canopy_bioemi_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_BIO( ZK, FCLAI, FCH, LAI, FSUN, PPFD_SUN, &
        PPFD_SHADE, TLEAF_SUN, TLEAF_SHADE, TLEAF_AVE, LU_OPT, &
        VTYPE, MODRES, CCE, VERT, CO2OPT, CO2SET, &
        EMI_IND, EMI_OUT)

!-----------------------------------------------------------------------

! Description:
!     computes parameterized canopy biogenic emissions

! Preconditions:
!     in-canopy FCLAI, model LAI, etc.

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/22 by PCC, based on Clifton et al. (2022) algorithms
! Citation:
! Clifton, O. E. et al. (2022). Large eddy simulation for investigating
! coupled forest canopy and turbulence influences on atmospheric chemistry.
! Journal of Advances in Modeling Earth Systems, 14, e2022MS003078.
! https://doi.org/10.1029/2022MS003078
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Jan 2023 P.C. Campbell: Initial canopy isoprene only version
!     Feb 2023 P.C. Campbell: Modified for multiple biogenic species
!     Jul 2023 P.C. Campbell: Restructured to use FSUN, TLEAF, and PPFD
!                             as inputs
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk,rgasuniv   !constants for canopy models
        use canopy_utils_mod,  ONLY: interp_linear1_internal, GET_GAMMA_CO2
        use canopy_bioparm_mod

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )       :: ZK(:)           ! Model heights (m)
        REAL(RK),    INTENT( IN )       :: FCLAI(:)        ! Fractional (z) shapes of the
        ! plant surface distribution (nondimensional), i.e., a Fractional Culmulative LAI
        REAL(RK),    INTENT( IN )       :: FCH             ! Canopy height (m)
        REAL(RK),    INTENT( IN )       :: LAI             ! Total Leaf Area Index
        REAL(RK),    INTENT( IN )       :: FSUN(:)         ! Sunlit/Shaded fraction from photolysis correction factor
        REAL(RK),    INTENT( IN )       :: PPFD_SUN(:)     ! PPFD for sunlit leaves (umol phot/m2 s)
        REAL(RK),    INTENT( IN )       :: PPFD_SHADE(:)   ! PPFD for shaded leaves (umol phot/m2 s)
        REAL(RK),    INTENT( IN )       :: TLEAF_SUN(:)    ! Leaf temp for sunlit leaves (K)
        REAL(RK),    INTENT( IN )       :: TLEAF_SHADE(:)  ! Leaf temp for shaded leaves (K)
        REAL(RK),    INTENT( IN )       :: TLEAF_AVE(:)    ! Average Leaf temp (K)
        INTEGER,     INTENT( IN )       :: LU_OPT          ! integer for LU type from model mapped to Massman et al. (default = 0/VIIRS)
        INTEGER,     INTENT( IN )       :: VTYPE           ! Grid cell dominant vegetation type
        REAL(RK),    INTENT( IN )       :: MODRES          ! Canopy model input vertical resolution (m)
        REAL(RK),    INTENT( IN )       :: CCE             ! MEGAN Canopy environment coefficient.
        INTEGER,     INTENT( IN )       :: VERT            ! MEGAN vertical integration option (default = 0/no integration)
        INTEGER,     INTENT( IN )       :: CO2OPT          ! Option for co2 inhibition calculation
        REAL(RK),    INTENT( IN )       :: CO2SET          ! User set atmospheric CO2 conc [ppmv]
        INTEGER,     INTENT( IN )       :: EMI_IND         ! Input biogenic emissions index
        REAL(RK),    INTENT( OUT )      :: EMI_OUT(:)      ! Output canopy layer volume emissions (kg m-3 s-1)

! Local Variables
        REAL(RK) :: TLEAF24_AVE(SIZE(ZK))          ! Average Leaf temp over the past 24 hours (K)
        REAL(RK) :: TLEAF240_AVE(SIZE(ZK))         ! Average Leaf temp over the past 240 hours (K)
        REAL(RK) :: GammaTLEAF_SUN_NUM(SIZE(ZK))   ! Numerator in Tleaf sun activity factor
        REAL(RK) :: GammaTLEAF_SHADE_NUM(SIZE(ZK)) ! Numerator in Tleaf shade activity factor
        REAL(RK) :: GammaTLEAF_SUN_DEN(SIZE(ZK))   ! Denominator in Tleaf sun activity factor
        REAL(RK) :: GammaTLEAF_SHADE_DEN(SIZE(ZK)) ! Denominator in Tleaf sun activity factor
        REAL(RK) :: GammaTLEAF_SUN(SIZE(ZK))       ! Tleaf sun activity factor
        REAL(RK) :: GammaTLEAF_SHADE(SIZE(ZK))     ! Tleaf shade activity factor
        REAL(RK) :: GammaTLEAF_AVE(SIZE(ZK))       ! Average Tleaf activity factor
        REAL(RK) :: PPFD24_SUN(SIZE(ZK))           ! Average PPFD sun over the past 24 hours
        REAL(RK) :: PPFD24_SHADE(SIZE(ZK))         ! Average PPFD shade over the past 24 hours
        REAL(RK) :: PPFD240_SUN(SIZE(ZK))          ! Average PPFD sun over the past 240 hours
        REAL(RK) :: PPFD240_SHADE(SIZE(ZK))        ! Average PPFD shade over the past 240 hours
        REAL(RK) :: CP_SUN(SIZE(ZK))               ! Normalized emission capacity sun at PPFD = 1000 umol phot/m2 s
        REAL(RK) :: CP_SHADE(SIZE(ZK))             ! Normalized emission capacity shade at PPFD = 1000 umol phot/m2 s
        REAL(RK) :: ALPHA_P_SUN(SIZE(ZK))          ! Quantum yield of isoprene sunlit (mol/mol)
        REAL(RK) :: ALPHA_P_SHADE(SIZE(ZK))        ! Quantum yield of isoprene shade (mol/mol)
        REAL(RK) :: GammaPPFD_SUN(SIZE(ZK))        ! PPFD activity factor sun (unitless)
        REAL(RK) :: GammaPPFD_SHADE(SIZE(ZK))      ! PPFD activity factor shade (unitless)
        REAL(RK) :: GammaPPFD_AVE(SIZE(ZK))        ! PPFD activity factor ave sun and shade
        REAL(RK) :: E_OPT(SIZE(ZK))                ! maximum normalized emission capacity
        REAL(RK) :: TLEAF_OPT(SIZE(ZK))            ! Tleaf at which E_OPT occurs (K)
        REAL(RK) :: FLAI(SIZE(ZK))                 ! Fractional LAI in layer
        REAL(RK) :: VPGWT(SIZE(ZK))                ! MEGANv3-like in-canopy weighting factor
        REAL(RK) :: GAUSS(SIZE(ZK))                ! MEGANv3-like in-canopy gaussian
        REAL(RK) :: CT1                            ! Activation energy (kJ/mol)
        REAL(RK) :: CEO                            ! Empirical coefficient
        REAL(RK) :: EF                             ! Final Mapped Emission factor (EF) (ug/m2 hr)
        REAL(RK) :: GAMMACO2                       ! CO2 inhibition factor (isoprene only)
        integer i, LAYERS

! Constant Canopy Parameters
        REAL(RK),          PARAMETER     :: PPFD0_SUN       =  200.0      !Constant PPFDo sunlit (umol/m2 s) (Guenther et al.,2012)
        REAL(RK),          PARAMETER     :: PPFD0_SHADE     =  50.0       !Constant PPFDo shaded (umol/m2 s) (Guenther et al.,2012)
        REAL(RK),          PARAMETER     :: CT2             =  230.0_rk   !Deactivation energy (kJ/mol) (Guenther et al., 2012)

! Calculate maximum normalized emission capacity (E_OPT) and Tleaf at E_OPT
        TLEAF240_AVE   = TLEAF_AVE  !Assume instantaneous TLEAF estimate for TLEAF240 and TLEAF24 (could improve...)
        TLEAF24_AVE    = TLEAF_AVE
        TLEAF_OPT = 313.0_rk + (0.6_rk * (TLEAF240_AVE-297.0_rk)) !Guenther et al. (2012)

! Calculate emission species/plant-dependent mapped emission factors
        call canopy_biop(EMI_IND, LU_OPT, VTYPE, EF, CT1, CEO)

        E_OPT = CEO * EXP(0.05_rk * (TLEAF24_AVE-297.0_rk)) * EXP(0.05_rk * (TLEAF240_AVE-297.0_rk))

! Calculate gamma (activity) values for average Tleaf (Clifton et al., 2022)
        GammaTLEAF_SUN_NUM = CT2*exp((CT1/(rgasuniv/1000.0))*((1.0/TLEAF_OPT)-(1.0/TLEAF_SUN)))
        GammaTLEAF_SUN_DEN = (CT2-CT1)*(1.0-exp((CT2/(rgasuniv/1000.0))*((1.0/TLEAF_OPT)-(1.0/TLEAF_SUN))))
        GammaTLEAF_SUN     = E_OPT*(GammaTLEAF_SUN_NUM/GammaTLEAF_SUN_DEN)

        GammaTLEAF_SHADE_NUM = CT2*exp((CT1/(rgasuniv/1000.0))*((1.0/TLEAF_OPT)-(1.0/TLEAF_SHADE)))
        GammaTLEAF_SHADE_DEN = (CT2-CT1)*(1.0-exp((CT2/(rgasuniv/1000.0))*((1.0/TLEAF_OPT)-(1.0/TLEAF_SHADE))))
        GammaTLEAF_SHADE     = E_OPT*(GammaTLEAF_SHADE_NUM/GammaTLEAF_SHADE_DEN)

        GammaTLEAF_AVE = (GammaTLEAF_SUN*FSUN) + (GammaTLEAF_SHADE*(1.0-FSUN)) ! average = sum sun and shade weighted by sunlit fraction

! Calculate gamma (activity) values for average PPFD (Clifton et al., 2022)
        PPFD240_SUN   = PPFD_SUN/2.0  !Clifton et al...halve the instantaneous PPFD estimate to get PPFD240 and PPFD24
        PPFD240_SHADE = PPFD_SHADE/2.0
        PPFD24_SUN    = PPFD_SUN/2.0
        PPFD24_SHADE  = PPFD_SHADE/2.0

        ALPHA_P_SUN = 0.004 - 0.0005*log(PPFD240_SUN)
        ALPHA_P_SHADE = 0.004 - 0.0005*log(PPFD240_SHADE)
        CP_SUN = 0.0468*(PPFD240_SUN**(0.6))*exp(0.005*(PPFD24_SUN-PPFD0_SUN))
        CP_SHADE = 0.0468*(PPFD240_SHADE**(0.6))*exp(0.005*(PPFD24_SHADE-PPFD0_SHADE))
        GammaPPFD_SUN   = CP_SUN*((ALPHA_P_SUN*PPFD_SUN)/SQRT(1.0 + (ALPHA_P_SUN**2.0) * (PPFD_SUN**2.0)))
        GammaPPFD_SHADE = CP_SHADE*((ALPHA_P_SHADE*PPFD_SHADE)/SQRT(1.0 + (ALPHA_P_SHADE**2.0) * (PPFD_SHADE**2.0)))

        GammaPPFD_AVE = (GammaPPFD_SUN*FSUN) + (GammaPPFD_SHADE*(1.0-FSUN)) ! average = sum sun and shade weighted by sunlit fraction

! Get CO2 inhibition factor for isoprene only

        if (EMI_IND .eq. 1) then  !Isoprene
            GAMMACO2 = GET_GAMMA_CO2(CO2OPT,CO2SET)
        else
            GAMMACO2 = 1.0_rk
        end if

! Calculate emissions profile in the canopy
        EMI_OUT = 0.0_rk  ! set initial emissions profile to zero
        FLAI = 0.0_rk  ! set initial fractional FLAI (LAD) profile to zero

        if (VERT .eq. 0) then         !Full 3D leaf-level biogenic emissions (no averaging, summing, or integration)
            do i=1, SIZE(ZK)
                if (ZK(i) .gt. 0.0 .and. ZK(i) .le. FCH) then           ! above ground level and at/below canopy top
                    FLAI(i) = ((FCLAI(i+1) - FCLAI(i)) * LAI)/MODRES    !fractional LAI in each layer converted to LAD (m2 m-3)
                    EMI_OUT(i) = FLAI(i) * EF * GammaTLEAF_AVE(i) * GammaPPFD_AVE(i) * GAMMACO2 * CCE  ! (ug m-3 hr-1)
                    EMI_OUT(i) = EMI_OUT(i) * 2.7777777777778E-13_rk    !convert emissions output to (kg m-3 s-1)
                end if
            end do
        else if (VERT .eq. 1) then       !"MEGANv3-like": Use weighting factors normalized to plant distribution shape (FCLAI)
            !across canopy layers
            LAYERS = floor(FCH/MODRES) + 1
            do i=1,  SIZE(ZK)
                if (ZK(i) .gt. 0.0 .and. ZK(i) .le. FCH) then
                    FLAI(i) = ((FCLAI(i+1) - FCLAI(i)) * LAI)/MODRES    !fractional LAI in each layer converted to LAD (m2 m-3)
                end if
            end do
            do i=1,  SIZE(ZK)
                if (ZK(i) .gt. 0.0 .and. ZK(i) .le. FCH) then
                    VPGWT(i) = (FLAI(i))/sum(FLAI(1:LAYERS))
                end if
            end do
            EMI_OUT(SIZE(ZK)) = LAI * EF * SUM(GammaTLEAF_AVE(1:LAYERS) * GammaPPFD_AVE(1:LAYERS) * &
                VPGWT(1:LAYERS)) * GAMMACO2 * CCE   !put into top model layer (ug m-2 hr-1)
            EMI_OUT = EMI_OUT * 2.7777777777778E-13_rk    !convert emissions output to    (kg m-2 s-1)
        else if (VERT .eq. 2) then       !"MEGANv3-like": Add weighted sum of activity coefficients using normal distribution
            !across canopy layers using 5 layer numbers directly from MEGANv3
            !--warning: weights are not consistent with FCLAI distribution
            !used for biomass distribution used for sunlit/shaded in Gamma TLEAF and GammaPPFD.
            LAYERS = floor(FCH/MODRES) + 1
            do i=1, SIZE(ZK)
                if (ZK(i) .gt. FCH) then
                    GAUSS(i) = 0.0
                else if (ZK(i) .le. FCH .and. ZK(i) .gt. FCH*(4.0_rk/5.0_rk)) then  !Level 1 - 2
                    GAUSS(i)   = interp_linear1_internal((/ FCH*(4.0_rk/5.0_rk),FCH /), &
                        (/ 0.118464_rk,0.0_rk /),ZK(i))
                else if (ZK(i) .le. FCH*(4.0_rk/5.0_rk) .and. ZK(i) .gt. FCH*(3.0_rk/5.0_rk)) then  !Level 2 - 3
                    GAUSS(i)   = interp_linear1_internal((/ FCH*(3.0_rk/5.0_rk),FCH*(4.0_rk/5.0_rk) /), &
                        (/ 0.239314_rk,0.118464_rk /),ZK(i))
                else if (ZK(i) .le. FCH*(3.0_rk/5.0_rk) .and. ZK(i) .gt. FCH*(2.0_rk/5.0_rk)) then  !Level 3 - 4
                    GAUSS(i)   = interp_linear1_internal((/ FCH*(2.0_rk/5.0_rk),FCH*(3.0_rk/5.0_rk) /), &
                        (/ 0.284444_rk,0.239314_rk /),ZK(i))
                else if (ZK(i) .le. FCH*(2.0_rk/5.0_rk) .and. ZK(i) .gt. FCH*(1.0_rk/5.0_rk) ) then  !Level 4 - Bottom
                    GAUSS(i)   = interp_linear1_internal((/ FCH*(1.0_rk/5.0_rk),FCH*(2.0_rk/5.0_rk) /), &
                        (/ 0.239314_rk,0.284444_rk /),ZK(i))
                else if (ZK(i) .le. FCH*(1.0_rk/5.0_rk) ) then  !Level 4 - Bottom
                    GAUSS(i)   = interp_linear1_internal((/ ZK(1),FCH*(1.0_rk/5.0_rk) /), &
                        (/ 0.118464_rk,0.239314_rk /),ZK(i))
                end if
            end do

            do i=1, SIZE(ZK)
                VPGWT(i) = GAUSS(i)/sum(GAUSS(1:LAYERS))
            end do
            EMI_OUT(SIZE(ZK)) = LAI * EF * SUM(GammaTLEAF_AVE(1:LAYERS) * GammaPPFD_AVE(1:LAYERS) * &
                VPGWT(1:LAYERS)) * GAMMACO2 * CCE   !put into top model layer (ug m-2 hr-1)
            EMI_OUT = EMI_OUT * 2.7777777777778E-13_rk    !convert emissions output to    (kg m-2 s-1)
        else if (VERT .eq. 3) then       !"MEGANv3-like": Add weighted sum of activity coefficients equally
            !across canopy layers
            !--warning: weights are not consistent with FCLAI distribution
            !used for biomass distribution used for sunlit/shaded in Gamma TLEAF and GammaPPFD.
            LAYERS = floor(FCH/MODRES) + 1
            do i=1,  SIZE(ZK)
                VPGWT(i) = 1.0_rk/LAYERS
            end do
            EMI_OUT(SIZE(ZK)) = LAI * EF * SUM(GammaTLEAF_AVE(1:LAYERS) * GammaPPFD_AVE(1:LAYERS) * &
                VPGWT(1:LAYERS)) * GAMMACO2 * CCE   !put into top model layer (ug m-2 hr-1)
            EMI_OUT = EMI_OUT * 2.7777777777778E-13_rk    !convert emissions output to    (kg m-2 s-1)
        else
            write(*,*)  'Wrong BIOVERT_OPT choice of ', VERT, ' in namelist...exiting'
            call exit(2)
        end if

    END SUBROUTINE CANOPY_BIO

end module canopy_bioemi_mod
