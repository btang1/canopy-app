!============================================================================================================
!    Program:    Simplified ACCESS Model
!                Simplified (Dry Deposition only) Atmospheric Chemistry and Canopy Exchange Simulation System
!
!    Developed by : Dr. Beiming Tang based on Dr. Rick Saylor's ACCESS model version 3.1.0
!============================================================================================================

module canopy_drydep_constants
    
    implicit none
    integer, parameter :: rk = selected_real_kind(15,307)
    integer, parameter :: ik = selected_real_kind(9)
!===========================================
!part 1.define constants used in calculation
!===========================================

    integer(kind=ik), parameter             :: ninteg = 31        !# of integrated species, depend on chemical mechanism used in ACCESS
    integer(kind=ik), parameter             :: ntotal = 31        !# of total species in ACCESS based on chosen chemical mechanim, including transported species
    integer(kind=ik), parameter             :: npts   = 8        !# of vertical layers
    integer(kind=ik), parameter             :: ncnpy  = 8        !# of vertical layers


!==================================================
!part 2.define global variables used in calculation
!==================================================

    real(kind=rk), parameter                :: tsoilk = 298      !soil/litter temperature(K)  
    real(kind=rk), parameter                :: stheta = 0.01     !volumetric soil water content (m^3/m^3)
    real(kind=rk), parameter                :: sattheta = 0.1    !saturation volumetric soil water content (m^3/m^3)
    real(kind=rk), parameter                :: dsoil = 4         !depth of topsoil (cm)
    real(kind=rk), parameter                :: rtheta = 0.05     !residual volumetic soil water content (m^3/m^3)
    real(kind=rk), parameter                :: sbcoef = 0.2      !clapp and hornberger exponent

end module canopy_drydep_constants
