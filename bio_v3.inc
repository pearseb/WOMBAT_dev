!

subroutine bio_v3(isc, iec, jsc, jec, isd, ied, jsd, jed, T_prog, Grid, Time, dtts, Thickness, Dens, swflx, sw_frac_zt, &
                  days_in_this_year)


! Based on the pt_npzd.aos version  - which is the optimized version
!
!
!     (\___/)
!     / o o \
!    (   "   )
!     \__ __/
!     
!
! World Ocean Model of Biogeochemistry And Trophic-dynamics (WOMBAT)
!
!  Authors:
!   Richard Matear          (Richard.Matear@csiro.au)
!   Matthew Chamberlain     (Matt.Chamberlain@csiro.au)
!   Pearse Buchanan         (Pearse.Buchanan@csiro.au)
!
!
!   Sequence
!   --------
!       1.    Set local variables and constants
!       2.    Retrieve key biotic parameters from input file 
!       3.    Get basic information (time, indices, 
!       4.    Initialise diagnostic arrays
!       5.    Retrieve indices of the active tracers
!       6.    Retrieve active tracer values in temporary arrays
!       7.    Create tracer mask 
!
!                [ BEGIN do loop over k and i indices ]
!
!       8.    Do light attenuation through water column 
!       9.    Calculate depth of the euphotic zone
!       10.   Aggregate light terms within mixed layer and calculate phytoplankton maximum growth rate
!
!                [ END do loop over k and i indices ]
!
!                [ BEGIN do loop over t, k and i indices ]
!
!       11.   Begin nested time-stepping for ecosystem model
!       11.1  Update the tracer values within the nested timestep
!       11.2  Get important ratios
!       11.3  Dissolved Fe chemistry 
!       11.4  Nutrient limitation terms for phytoplankton growth  
!       11.5  Light limitation terms for phytoplankton growth  
!       11.6  Realised growth rate given resource availability
!       11.7  Growth in chlorophyll (mg Chl / m3)
!       11.8  Phytoplankton uptake of Fe (nmol Fe / m3)
!       11.9  Temperature dependence of heterotrophy / respiration
!       11.10 Grazing functions 
!       11.11 Mortality scalers
!       11.12 Collect major source and sink terms for tendency equations
!       11.13 Apply sources and sinks to tendency equations 
!       11.14 Collect key diagnostics
!       11.15 Check for conservation of mass by ecosystem component
!
!               [ END do loop over t, k and i indices ]
!
!
!               [ BEGIN do loop over k and i indices ]
!
!       12.   Add biotically-induced tendency to biotracers
!       13.   Collect other diagnostic terms outside of the nested ecosystem time-step
!
!               [ END do loop over k and i indices ]
!
!
!               [ BEGIN do loop over k and i indices ]
!
!       14.   Sink detritus, detrital iron, detrital silicate and CaCO3
!
!               [ END do loop over k and i indices ]
!
!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                                             :: isc, iec
integer, intent(in)                                             :: jsc, jec
integer, intent(in)                                             :: isd, ied
integer, intent(in)                                             :: jsd, jed
type(ocean_prog_tracer_type), dimension(:), intent(inout)       :: T_prog
type(ocean_grid_type), intent(in)                               :: Grid
type(ocean_time_type), intent(in)                               :: Time
real, intent(in)                                                :: dtts
type(ocean_thickness_type), intent(in)                          :: Thickness
type(ocean_density_type), intent(in)                            :: Dens
real, intent(in), dimension(isd:ied,jsd:jed)                    :: swflx        ! short wave radiation flux (W/m^2)
real, intent(in), dimension(isd:,jsd:,:)                        :: sw_frac_zt   ! short wave fraction on T grid (none)
real, intent(in)                                                :: days_in_this_year


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [1] - Set local variables and constants                                                 !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

integer :: i
integer :: j
integer :: k
integer :: n

logical :: used
integer :: index_temp, index_salt
integer :: ts_npzd                      ! number of time steps within NPZD model
integer :: tn, trn
integer :: ichl                         ! index for chlorophyll-dependent light attenuation coefficient table

! ---------------------- !
!!! ------------------ !!!
!!! Time and constants !!!
!!! ------------------ !!!
! ---------------------- !
real :: pi = 3.14159265358979           ! pi
real :: epsi = 1e-16                    ! small number for ensuring division by zero doesn't happen
real :: rdtts                           ! 1/dtts
real :: dtsb                            ! number of seconds in the time-split of the ecosystem component
real :: day_angle                       ! Fractional year in radians as a function of days_in_this_year
real :: declination                     ! Angle between the rays of the Sun and the plane of the Earth's equator
real :: cos_hour_angle                  ! A measure of time since solar noon (radians)

! ---------------------------- !
!!! ------------------------ !!!
!!! Temperature and salinity !!!
!!! ------------------------ !!!
! ---------------------------- !
real :: zval                            ! Generic value. Used many times in different parts of this code.
real :: ztem                            ! Temperature (deg C)
real :: zsal                            ! Salinity (psu)
real :: ztemk                           ! Temperature (Kelvin)
real :: mld(isd:ied,jsd:jed)            ! Mixed layer depth

!-------------------!
!!! ------------- !!!
!!! Tracer values !!!
!!! ------------- !!!
!-------------------!
real :: biodic                          ! MAX(0.0, dissolved inorganic carbon)
real :: biopo4                          ! MAX(0.0, phosphate)
real :: biono3                          ! MAX(0.0, nitrate)
real :: biofe                           ! MAX(0.0, dissolved iron)
real :: biooxy                          ! MAX(0.0, oxygen)
real :: biophy                          ! MAX(0.0, nanophytoplankton)
real :: biozoo                          ! MAX(0.0, microzooplankton)
real :: biodet                          ! MAX(0.0, slow-sinking detritus)
real :: biocaco3                        ! MAX(0.0, calcium carbonate)
real :: biopchl                         ! MAX(0.0, nanophytoplankton chlorophyll)
real :: biophyfe                        ! MAX(0.0, nanophytoplankton iron)
real :: biozoofe                        ! MAX(0.0, microzooplankton iron)
real :: biodetfe                        ! MAX(0.0, slow-sinking detrital iron)
real :: biophy_min                      ! MAX(some_threshold, nanophytoplankton)     || beneath some_threshold no mortality occurs
real :: biodet_min                      ! MAX(some_threshold, detritus)              || beneath some_threshold no mortality occurs
real :: biomask(isc:iec,grid%nk)        ! Land-sea masking
real :: biotr(isc:iec,grid%nk,ntr_bmax) ! Biogeochemical tracer array

! --------------------------------------------- !
!!! ----------------------------------------- !!!
!!! Limitation terms for phytoplankton growth !!!
!!! ----------------------------------------- !!!
! --------------------------------------------- !
real :: phy_limnut                      ! Nanophytoplankton limitation by all nutrients   || 0 == full limitation
real :: phy_conc1                       ! Nanophytoplankton biomass exceeding bloom threshold
real :: phy_conc2                       ! Nanophytoplankton biomass (but maximum == bloom threshold)
real :: phy_k_nit                       ! Nanophytoplankton half-saturation coefficient for dissolved inorganic nitrogen uptake
real :: phy_k_fer                       ! Nanophytoplankton half-saturation coefficient for dissolved iron uptake
real :: phy_k_pho                       ! Nanophytoplankton half-saturation coefficient for phosphate uptake
real :: phy_limnit(isc:iec,grid%nk)     ! Nanophytoplankton limitation by nitrogen
real :: phy_limpho(isc:iec,grid%nk)     ! Nanophytoplankton limitation by phosphorus
real :: phy_limfe(isc:iec,grid%nk)      ! Nanophytoplankton limitation by iron

! ------------------------------------- !
!!! --------------------------------- !!!
!!! Mortality terms for phytoplankton !!!
!!! --------------------------------- !!!
! ------------------------------------- !
real :: phy_slmor                       ! Scaler [0-1] to depress linear nanophytoplankton mortality in nutrient-limited regimes

! ---------------------------- !
!!! ------------------------ !!!
!!! Growth and grazing rates !!!
!!! ------------------------ !!!
! ---------------------------- !
real :: u_phy                           ! Realised growth rate of nanophytoplankton (nutrients & light limitations)
real :: g_zoo                           ! Realised grazing rate of microzooplankton
real :: zooprey                         ! Total prey biomass for microzooplankton
real :: fbc                             ! Heterotrophic temperature-dependence function
real :: mumax_phy(isc:iec,grid%nk)      ! Maximum potential (T-dependent) growth rate of nanophytoplankton
real :: mu_phy(isc:iec,grid%nk)         ! Realised growth rate of nanophytoplankton (light limitation)

! --------------------------------------- !
!!! ----------------------------------- !!!
!!! Tendency terms for tracer equations !!!
!!! ----------------------------------- !!!
! --------------------------------------- !
real :: f11phy                          ! Nanophytoplankton biomass growth
real :: z21phy                          ! Microzooplankton grazing of nanophytoplankton
real :: z21det                          ! Microzooplankton grazing of slow-sinking detritus
real :: f22phy                          ! Nanophytoplankton biomass loss via linear mortality (respiration)
real :: f23phy                          ! Nanophytoplankton biomass loss via quadratic mortality (viral lysis, etc.)
real :: z31                             ! Microzooplankton biomass loss via linear mortality (excretion)
real :: z32                             ! Mesozooplankton biomass loss via quadratic mortality (predation, infection, etc.)
real :: f41det                          ! Remineralisation of slow-sinking detritus to inorganic constituents
real :: f51                             ! Remineralisation of CaCO3 to DIC and Alkalinity
real :: caco3_bgc_change                ! Total change in CaCO3 following time-splitting of ecosystem component
real :: no3_bgc_change                  ! Total change in nitrate following time-splitting of ecosystem component
real :: c_pools(isc:iec,grid%nk,2)      ! Total biological carbon budget (PHY + ZOO + DET + DIC + CaCO3)
real :: n_pools(isc:iec,grid%nk,2)      ! Total biological nitrogen budget (PHY + ZOO + DET + NO3)

! ------------------------- !
!!! --------------------- !!!
!!! Light and chlorophyll !!!
!!! --------------------- !!!
! ------------------------- !
real :: zchl                            ! Total chlorophyll in mg/m3  || used to find light attenuation coefficients
real :: phy_pisl2                       ! Photosynthetic-Irradiance slope of nanophytoplankton (~ Chl:C, respiration, daylength, zeuph/zmld)
real :: phy_chl_pisl                    ! Photosynthetic-Irradiance slope of nanophytoplankton chlorophyll (~ Chl:C, µmax, daylength, zeuph/zmld)
real :: mumin_pchl                      ! Rate of chlorophyll growth in nanophytoplankton to achieve minimum Chl:C ratio
real :: muopt_pchl                      ! Rate of chlorophyll growth in nanophytoplankton to achieve optimal Chl:C ratio
real :: keuphot(isc:iec)                ! Estimated depth of the euphotic zone
real :: daylight(isc:iec)               ! Number of daylight hours in a day
real :: par_phy_mldsum(isc:iec)         ! Sum of photosynthetic active radiation seen by nanophytoplankton in the mixed layer
real :: par_tot_mldsum(isc:iec)         ! Sum of photosynthetic active radiation in the mixed layer
real :: par_z_mldsum(isc:iec)           ! Sum of depths in the mixed layer
real :: par_tot(isc:iec,grid%nk)        ! Total photosynthetic active radiation (red + green + blue)
real :: par_phy(isc:iec,grid%nk)        ! Photosynthetic active radiation seen by nanophytoplankton (red + green + blue)
real :: par_phymld(isc:iec,grid%nk)     ! Average photosynthetic active radiation seen by nanophytoplankton (red + green + blue) in mixed layer
real :: par_eup(isc:iec,grid%nk)        ! Photosynthetic active radiation for calculation of euphotic zone depth
real :: phy_pisl(isc:iec,grid%nk)       ! Initial guess at photosynthesis-irradiance slope (~ Chl:C) for nanophytoplankton
real :: phy_lday(isc:iec,grid%nk)       ! Daylength and euphotic < MLD limitation term for nanophytoplankton
real :: chl_lday(isc:iec,grid%nk)       ! Daylength and euphotic < MLD limitation term for chlorophyll production ????
real :: phy_lpar(isc:iec,grid%nk)       ! Light limitation term for nanophytoplankton
real :: phy_chl_lpar(isc:iec,grid%nk)   ! Light limitation term for nanophytoplankton chlorophyll
real :: phy_chlc(isc:iec,grid%nk)       ! Nanophytoplankton Chl:C ratio
real :: mu_pchl(isc:iec,grid%nk)        ! Growth rate of nanophytoplankton chlorophyll
real :: ek_rgb(isc:iec,grid%nk,3)       ! Attenutation coefficients for red, green and blue light
real :: par_rgb(isc:iec,grid%nk,3)      ! Photosynthetic active radiation split into red, green and blue light wavelength bands
real :: zrgb(4,61)                      ! Table of attenutation coefficients for red, green and blue light dependent on chlorophyll concentration


! -------------- !
!!! ---------- !!!
!!! Iron cycle !!!
!!! ---------- !!!
! -------------- !
real :: phy_minFeC                      ! Minimum Fe:C quota of nanophytoplankton given their Chlorophyll content and N limitations
real :: phy_maxQFe                      ! Maximum Fe:C quota of nanophytoplankton
real :: phy_Feupt_upreg                 ! Upregulation of Fe uptake by nanophytoplankton in response to low Fe cellular quota
real :: phy_Feupt_downreg               ! Downregulation of Fe uptake by nanophytoplankton in response to achieving high cellular Fe quotas
real :: ligand                          ! Concentration of ligand in the ocean
real :: fe_keq                          ! Equilibrium coefficient determining fraction of dFe in Fe' or L-Fe (ligand bound) (~ temperature)
real :: fe_III                          ! Estimated concentration of Fe' (ferric iron form, Fe(III))
real :: fe_lig                          ! Amount of dFe bound to organic ligands
real :: fe_col                          ! Amount of dFe within the colloidal form
real :: fesol1                          ! Fe solubility constant 1 (~ Salinity and Temperature)
real :: fesol2                          ! Fe solubility constant 2 (~ Salinity and Temperature)
real :: fesol3                          ! Fe solubility constant 3 (~ Salinity and Temperature)
real :: fesol4                          ! Fe solubility constant 4 (~ Salinity and Temperature)
real :: fesol5                          ! Fe solubility constant 5 (~ Salinity and Temperature)
real :: hp                              ! concentration of hydrogen ions (H+) derived from pH
real :: fe3sol                          ! Solubility of Fe(III) as Fe' (free iron)
real :: precip                          ! Rate of Fe' precipitation to Fe hydroxides and oxyhydroxides (permanent loss of dFe)
real :: partic                          ! Concentration of particles (Carbon detritus + CaCO3 + BSi)
real :: scaven                          ! Rate of scavenging of Fe' (Fe(III)) to all particles
real :: scadet                          ! Rate of scavenging of Fe' (Fe(III)) to detritus
real :: feloss                          ! Rate of ambiguous dFe losses in near-coastal environments (i.e., amplified scavenging, etc.)
real :: fe2det                          ! Rate of coagulation of colloidal dFe to detritus
real :: phy_FeC(isc:iec,grid%nk)        ! Nanophytoplankton Fe:C ratio
real :: zoo_FeC(isc:iec,grid%nk)        ! Microzooplankton Fe:C ratio
real :: det_FeC(isc:iec,grid%nk)        ! Slow-sinking detritus Fe:C ratio
real :: phy_Feupt(isc:iec,grid%nk)      ! Nanophytoplankton Fe uptake rate

! -------------------- !
!!! ---------------- !!!
!!! Detrital sinking !!!
!!! ---------------- !!!
! -------------------- !
real :: biodetsed(isc:iec)              ! Temporary array for rate of accumulation in sediments 
real :: biodetfesed(isc:iec)            ! Temporary array for rate of accumulation in sediments
real :: adv_fb1(isc:iec,1:grid%nk+1)    ! Sinking of detritus (carbon)
real :: adv_fb2(isc:iec,1:grid%nk+1)    ! Sinking of detritus (iron)


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read in attenuation coefficients for blue, green and red light depending on chlorophyll !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Chlorophyll      ! Blue attenuation    ! Green attenuation   ! Red attenuation
      zrgb(1, 1) =  0.010; zrgb(2, 1) = 0.01618; zrgb(3, 1) = 0.07464; zrgb(4, 1) = 0.3780
      zrgb(1, 2) =  0.011; zrgb(2, 2) = 0.01654; zrgb(3, 2) = 0.07480; zrgb(4, 2) = 0.37823
      zrgb(1, 3) =  0.013; zrgb(2, 3) = 0.01693; zrgb(3, 3) = 0.07499; zrgb(4, 3) = 0.37840
      zrgb(1, 4) =  0.014; zrgb(2, 4) = 0.01736; zrgb(3, 4) = 0.07518; zrgb(4, 4) = 0.37859
      zrgb(1, 5) =  0.016; zrgb(2, 5) = 0.01782; zrgb(3, 5) = 0.07539; zrgb(4, 5) = 0.37879
      zrgb(1, 6) =  0.018; zrgb(2, 6) = 0.01831; zrgb(3, 6) = 0.07562; zrgb(4, 6) = 0.37900
      zrgb(1, 7) =  0.020; zrgb(2, 7) = 0.01885; zrgb(3, 7) = 0.07586; zrgb(4, 7) = 0.37923
      zrgb(1, 8) =  0.022; zrgb(2, 8) = 0.01943; zrgb(3, 8) = 0.07613; zrgb(4, 8) = 0.37948
      zrgb(1, 9) =  0.025; zrgb(2, 9) = 0.02005; zrgb(3, 9) = 0.07641; zrgb(4, 9) = 0.37976
      zrgb(1,10) =  0.028; zrgb(2,10) = 0.02073; zrgb(3,10) = 0.07672; zrgb(4,10) = 0.38005
      zrgb(1,11) =  0.032; zrgb(2,11) = 0.02146; zrgb(3,11) = 0.07705; zrgb(4,11) = 0.38036
      zrgb(1,12) =  0.035; zrgb(2,12) = 0.02224; zrgb(3,12) = 0.07741; zrgb(4,12) = 0.38070
      zrgb(1,13) =  0.040; zrgb(2,13) = 0.02310; zrgb(3,13) = 0.07780; zrgb(4,13) = 0.38107
      zrgb(1,14) =  0.045; zrgb(2,14) = 0.02402; zrgb(3,14) = 0.07821; zrgb(4,14) = 0.38146
      zrgb(1,15) =  0.050; zrgb(2,15) = 0.02501; zrgb(3,15) = 0.07866; zrgb(4,15) = 0.38189
      zrgb(1,16) =  0.056; zrgb(2,16) = 0.02608; zrgb(3,16) = 0.07914; zrgb(4,16) = 0.38235
      zrgb(1,17) =  0.063; zrgb(2,17) = 0.02724; zrgb(3,17) = 0.07967; zrgb(4,17) = 0.38285
      zrgb(1,18) =  0.071; zrgb(2,18) = 0.02849; zrgb(3,18) = 0.08023; zrgb(4,18) = 0.38338
      zrgb(1,19) =  0.079; zrgb(2,19) = 0.02984; zrgb(3,19) = 0.08083; zrgb(4,19) = 0.38396
      zrgb(1,20) =  0.089; zrgb(2,20) = 0.03131; zrgb(3,20) = 0.08149; zrgb(4,20) = 0.38458
      zrgb(1,21) =  0.100; zrgb(2,21) = 0.03288; zrgb(3,21) = 0.08219; zrgb(4,21) = 0.38526
      zrgb(1,22) =  0.112; zrgb(2,22) = 0.03459; zrgb(3,22) = 0.08295; zrgb(4,22) = 0.38598
      zrgb(1,23) =  0.126; zrgb(2,23) = 0.03643; zrgb(3,23) = 0.08377; zrgb(4,23) = 0.38676
      zrgb(1,24) =  0.141; zrgb(2,24) = 0.03842; zrgb(3,24) = 0.08466; zrgb(4,24) = 0.38761
      zrgb(1,25) =  0.158; zrgb(2,25) = 0.04057; zrgb(3,25) = 0.08561; zrgb(4,25) = 0.38852
      zrgb(1,26) =  0.178; zrgb(2,26) = 0.04289; zrgb(3,26) = 0.08664; zrgb(4,26) = 0.38950
      zrgb(1,27) =  0.200; zrgb(2,27) = 0.04540; zrgb(3,27) = 0.08775; zrgb(4,27) = 0.39056
      zrgb(1,28) =  0.224; zrgb(2,28) = 0.04811; zrgb(3,28) = 0.08894; zrgb(4,28) = 0.39171
      zrgb(1,29) =  0.251; zrgb(2,29) = 0.05103; zrgb(3,29) = 0.09023; zrgb(4,29) = 0.39294
      zrgb(1,30) =  0.282; zrgb(2,30) = 0.05420; zrgb(3,30) = 0.09162; zrgb(4,30) = 0.39428
      zrgb(1,31) =  0.316; zrgb(2,31) = 0.05761; zrgb(3,31) = 0.09312; zrgb(4,31) = 0.39572
      zrgb(1,32) =  0.355; zrgb(2,32) = 0.06130; zrgb(3,32) = 0.09474; zrgb(4,32) = 0.39727
      zrgb(1,33) =  0.398; zrgb(2,33) = 0.06529; zrgb(3,33) = 0.09649; zrgb(4,33) = 0.39894
      zrgb(1,34) =  0.447; zrgb(2,34) = 0.06959; zrgb(3,34) = 0.09837; zrgb(4,34) = 0.40075
      zrgb(1,35) =  0.501; zrgb(2,35) = 0.07424; zrgb(3,35) = 0.10040; zrgb(4,35) = 0.40270
      zrgb(1,36) =  0.562; zrgb(2,36) = 0.07927; zrgb(3,36) = 0.10259; zrgb(4,36) = 0.40480
      zrgb(1,37) =  0.631; zrgb(2,37) = 0.08470; zrgb(3,37) = 0.10495; zrgb(4,37) = 0.40707
      zrgb(1,38) =  0.708; zrgb(2,38) = 0.09056; zrgb(3,38) = 0.10749; zrgb(4,38) = 0.40952
      zrgb(1,39) =  0.794; zrgb(2,39) = 0.09690; zrgb(3,39) = 0.11024; zrgb(4,39) = 0.41216
      zrgb(1,40) =  0.891; zrgb(2,40) = 0.10374; zrgb(3,40) = 0.11320; zrgb(4,40) = 0.41502
      zrgb(1,41) =  1.000; zrgb(2,41) = 0.11114; zrgb(3,41) = 0.11639; zrgb(4,41) = 0.41809
      zrgb(1,42) =  1.122; zrgb(2,42) = 0.11912; zrgb(3,42) = 0.11984; zrgb(4,42) = 0.42142
      zrgb(1,43) =  1.259; zrgb(2,43) = 0.12775; zrgb(3,43) = 0.12356; zrgb(4,43) = 0.42500
      zrgb(1,44) =  1.413; zrgb(2,44) = 0.13707; zrgb(3,44) = 0.12757; zrgb(4,44) = 0.42887
      zrgb(1,45) =  1.585; zrgb(2,45) = 0.14715; zrgb(3,45) = 0.13189; zrgb(4,45) = 0.43304
      zrgb(1,46) =  1.778; zrgb(2,46) = 0.15803; zrgb(3,46) = 0.13655; zrgb(4,46) = 0.43754
      zrgb(1,47) =  1.995; zrgb(2,47) = 0.16978; zrgb(3,47) = 0.14158; zrgb(4,47) = 0.44240
      zrgb(1,48) =  2.239; zrgb(2,48) = 0.18248; zrgb(3,48) = 0.14701; zrgb(4,48) = 0.44765
      zrgb(1,49) =  2.512; zrgb(2,49) = 0.19620; zrgb(3,49) = 0.15286; zrgb(4,49) = 0.45331
      zrgb(1,50) =  2.818; zrgb(2,50) = 0.21102; zrgb(3,50) = 0.15918; zrgb(4,50) = 0.45942
      zrgb(1,51) =  3.162; zrgb(2,51) = 0.22703; zrgb(3,51) = 0.16599; zrgb(4,51) = 0.46601
      zrgb(1,52) =  3.548; zrgb(2,52) = 0.24433; zrgb(3,52) = 0.17334; zrgb(4,52) = 0.47313
      zrgb(1,53) =  3.981; zrgb(2,53) = 0.26301; zrgb(3,53) = 0.18126; zrgb(4,53) = 0.48080
      zrgb(1,54) =  4.467; zrgb(2,54) = 0.28320; zrgb(3,54) = 0.18981; zrgb(4,54) = 0.48909
      zrgb(1,55) =  5.012; zrgb(2,55) = 0.30502; zrgb(3,55) = 0.19903; zrgb(4,55) = 0.49803
      zrgb(1,56) =  5.623; zrgb(2,56) = 0.32858; zrgb(3,56) = 0.20898; zrgb(4,56) = 0.50768
      zrgb(1,57) =  6.310; zrgb(2,57) = 0.35404; zrgb(3,57) = 0.21971; zrgb(4,57) = 0.51810
      zrgb(1,58) =  7.079; zrgb(2,58) = 0.38154; zrgb(3,58) = 0.23129; zrgb(4,58) = 0.52934
      zrgb(1,59) =  7.943; zrgb(2,59) = 0.41125; zrgb(3,59) = 0.24378; zrgb(4,59) = 0.54147
      zrgb(1,60) =  8.912; zrgb(2,60) = 0.44336; zrgb(3,60) = 0.25725; zrgb(4,60) = 0.55457
      zrgb(1,61) = 10.000; zrgb(2,61) = 0.47804; zrgb(3,61) = 0.27178; zrgb(4,61) = 0.56870


!
! =====================================================================
!     begin executable code
! =====================================================================

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [2] - Retrieve key biotic parameters from input file                                    !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

  call time_interp_external(alphabio_id, time%model_time, alphabio) !Initial slope of PI curve
  call time_interp_external(parbio_id, time%model_time, parbio)     !Proportion of SW rad that is PAR
  call time_interp_external(phybiot_id, time%model_time, phybiot)   !Phytoplankton biomass threshold
  call time_interp_external(phyminqc_id, time%model_time, phyminqc) !Phytoplankton minimum Chl:C quota
  call time_interp_external(phyoptqc_id, time%model_time, phyoptqc) !Phytoplankton optimal Chl:C quota
  call time_interp_external(phymaxqf_id, time%model_time, phymaxqf) !Phytoplankton maximum Fe:C quota
  call time_interp_external(phyoptqf_id, time%model_time, phyoptqf) !Phytoplankton optimal Fe:C quota
  call time_interp_external(abioa_phy_id, time%model_time, abioa_phy) !Autotrophic T-dependent scalar
  call time_interp_external(bbioa_id, time%model_time, bbioa)       !Autotrophic T-dependent base
  call time_interp_external(cbioa_id, time%model_time, cbioa)       !Autotrophic T-dependent exponent
  call time_interp_external(abioh_id, time%model_time, abioh)       !Heterotrophic T-dependent scalar
  call time_interp_external(bbioh_id, time%model_time, bbioh)       !Heterotrophic T-dependent base
  call time_interp_external(cbioh_id, time%model_time, cbioh)       !Heterotrophic T-dependent exponent
  call time_interp_external(phykn_id, time%model_time, phykn)       !Phytoplankton N half-saturation
  call time_interp_external(phykf_id, time%model_time, phykf)       !Phytoplankton Fe half-saturation
  call time_interp_external(phykp_id, time%model_time, phykp)       !Phytoplankton P half-saturation
  call time_interp_external(phylmor_id, time%model_time, phylmor)   !Phytoplankton linear mortality
  call time_interp_external(phyqmor_id, time%model_time, phyqmor)   !Phytoplankton quadratic mortality
  call time_interp_external(zooassi_id, time%model_time, zooassi)   !Zooplankton assimilation efficiency
  call time_interp_external(zoogmax_id, time%model_time, zoogmax)   !Zooplankton max grazing rate
  call time_interp_external(epszoo_id, time%model_time, epszoo)     !Zooplankton prey capture rate
  call time_interp_external(zprefphy_id, time%model_time, zprefphy) !Zooplankton preference for phytoplankton
  call time_interp_external(zprefdet_id, time%model_time, zprefdet) !Zooplankton preference for Det
  call time_interp_external(zooqmor_id, time%model_time, zooqmor)   !Zooplankton quadratic mortality
  call time_interp_external(zooexcr_id, time%model_time, zooexcr)   !Zooplankton excretion rate
  call time_interp_external(detlrem_id, time%model_time, detlrem)   !Detritus remineralisation rate
  call time_interp_external(wdetbio_id, time%model_time, wdetbio)   !Sinking rate of slow detritus
  call time_interp_external(caco3rem_id, time%model_time, caco3rem) !CaCO3 remineralisation rate
  call time_interp_external(wcaco3_id, time%model_time, wcaco3)     !CaCO3 sinking rate
  call time_interp_external(tscav_fe_id, time%model_time, tscav_fe) !Iron scavenging rate
  call time_interp_external(fe_bkgnd_id, time%model_time, fe_bkgnd) !Background dFe concentration
  call time_interp_external(f_inorg_id, time%model_time, f_inorg)   !CaCO3 inorganic fraction
  call time_interp_external(knano_dfe_id, time%model_time, knano_dfe) !Rate of dFe nanoparticle formation (dFe-->)
  call time_interp_external(kscav_dfe_id, time%model_time, kscav_dfe) !Rate of dFe scavenging (dFe-->particles)
  call time_interp_external(kcoag_dfe_id, time%model_time, kcoag_dfe) !Rate of dFe colloidal coagulation (dFe-->detritus)


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [3] - Get basic information (time, physical values)                                     !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

  ! Get indices of temperature and salinity
  index_temp = fm_get_index('/ocean_mod/prog_tracers/temp')
  index_salt = fm_get_index('/ocean_mod/prog_tracers/salt')

  ! Get the number of timesteps for the ecosystem model per model timestep
  ts_npzd = max(1, nint(dtts / 900.))

  ! Calculate inverse of total number of seconds per timestep
  rdtts  = 1/dtts

  ! Number of seconds per nested ecosystem timestep
  dtsb=dtts/float(ts_npzd)

  ! Calculate the mixed layer depth.
  call calc_mixed_layer_depth(Thickness,                                                           &
                              T_prog(index_salt)%field(isd:ied,jsd:jed,:,Time%tau),                &
                              T_prog(index_temp)%field(isd:ied,jsd:jed,:,Time%tau),                &
                              Dens%rho(isd:ied,jsd:jed,:,Time%tau),                                &
                              Dens%pressure_at_depth(isd:ied,jsd:jed,:), mld)

!  write (stdout(),*) ' AO-NPZD model will do ',ts_npzd,' time steps'
!  write (stdout(),*) ' time step in NPZD model will be ',  dtts/ts_npzd,'sec.'  

!
!
!-----------------------------------------------------------------------
!     calculate the source terms for BIOTICs
!-----------------------------------------------------------------------
!

!
!       Loop over multiple instances
!



do n = 1, instances  !{

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [4] - Initialise diagnostic arrays                                                      !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

  pprod_gross(:,:,:) = 0.0
  zprod_gross(:,:,:) = 0.0
  light_limit(:,:) = 0.0
  radbio3d(:,:,:) = 0.0
  phyxsize(:,:,:) = 0.0
  npp3d(:,:,:) = 0.0
  nsp3d(:,:,:) = 0.0
  adic_intmld(:,:) = 0.0
  dic_intmld(:,:) = 0.0
  o2_intmld(:,:) = 0.0
  no3_intmld(:,:) = 0.0
  fe_intmld(:,:) = 0.0
  phy_intmld(:,:) = 0.0
  det_intmld(:,:) = 0.0
  pprod_gross_intmld(:,:) = 0.0
  npp_intmld(:,:) = 0.0
  radbio_intmld(:,:) = 0.0
  adic_int100(:,:) = 0.0
  dic_int100(:,:) = 0.0
  o2_int100(:,:) = 0.0
  no3_int100(:,:) = 0.0
  fe_int100(:,:) = 0.0
  phy_int100(:,:) = 0.0
  det_int100(:,:) = 0.0
  pprod_gross_int100(:,:) = 0.0
  npp_int100(:,:) = 0.0
  radbio_int100(:,:) = 0.0
  zeuphot(:,:) = 10.0       ! minimum euphotic zone depth is 20 metres
  chlorophyll(:,:) = 0.0    ! total surface chlorophyll
  wdet(:,:,:) = 0.0         ! sinking rate of detritus
  phy_parlimit(:,:,:) = 0.0
  phy_Felimit(:,:,:) = 0.0
  phy_Nlimit(:,:,:) = 0.0
  phy_Plimit(:,:,:) = 0.0
  zoo_grazpres(:,:,:) = 0.0


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [5] - Retrieve indices of the active tracers                                            !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

  ind_po4  = biotic(n)%ind_bgc(id_po4)      ! mmol P / m3
  ind_no3  = biotic(n)%ind_bgc(id_no3)      ! mmol N / m3
  ind_fe   = biotic(n)%ind_bgc(id_fe)       ! µmol Fe / m3
  ind_o2   = biotic(n)%ind_bgc(id_o2)       ! mmol O2 / m3
  ind_dic  = biotic(n)%ind_bgc(id_dic)      ! mmol DIC / m3
  ind_alk  = biotic(n)%ind_bgc(id_alk)      ! mmol Eq / m3
  ind_caco3= biotic(n)%ind_bgc(id_caco3)    ! mmol C / m3
  ind_adic = biotic(n)%ind_bgc(id_adic)     ! mmol DIC / m3
  ind_phy  = biotic(n)%ind_bgc(id_phy)      ! mmol C / m3
  ind_zoo  = biotic(n)%ind_bgc(id_zoo)      ! mmol C / m3
  ind_det  = biotic(n)%ind_bgc(id_det)      ! mmol C / m3
  ind_pchl = biotic(n)%ind_bgc(id_pchl)     ! mg / m3
  ind_phyfe = biotic(n)%ind_bgc(id_phyfe)   ! mmol Fe / m3
  ind_zoofe = biotic(n)%ind_bgc(id_zoofe)   ! mmol Fe / m3
  ind_detfe = biotic(n)%ind_bgc(id_detfe)   ! mmol Fe / m3


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! biotic source-sink terms using Euler forward timestepping !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do j = jsc, jec  !{ 

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [6] - Retrieve active tracer values in temporary arrays (tau-1)                         !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

       do k=1,grid%nk
        do i=isc,iec

          ! ESSENTIAL TRACERS
          biotr(i,k,id_dic) = max(0.0,t_prog(ind_dic)%field(i,j,k,Time%taum1)) 
          biotr(i,k,id_no3) = max(0.0,t_prog(ind_no3)%field(i,j,k,Time%taum1)) 
          biotr(i,k,id_phy) = max(0.0,t_prog(ind_phy)%field(i,j,k,Time%taum1))
          biotr(i,k,id_zoo) = max(0.0,t_prog(ind_zoo)%field(i,j,k,Time%taum1))
          biotr(i,k,id_det) = max(0.0,t_prog(ind_det)%field(i,j,k,Time%taum1))
          biotr(i,k,id_o2) = max(0.0,t_prog(ind_o2)%field(i,j,k,Time%taum1)) 
          ! OPTIONAL TRACERS
          if (id_caco3.ne.0) biotr(i,k,id_caco3) = max(0.0,t_prog(ind_caco3)%field(i,j,k,Time%taum1)) 
          if (id_po4.ne.0)   biotr(i,k,id_po4) = max(0.0,t_prog(ind_po4)%field(i,j,k,Time%taum1)) 
          if (id_fe.ne.0)    biotr(i,k,id_fe) = max(0.0,t_prog(ind_fe)%field(i,j,k,Time%taum1)) 
          if (id_pchl.ne.0)  biotr(i,k,id_pchl) = max(0.0,t_prog(ind_pchl)%field(i,j,k,Time%taum1))
          if (id_phyfe.ne.0) biotr(i,k,id_phyfe) = max(0.0,t_prog(ind_phyfe)%field(i,j,k,Time%taum1))
          if (id_zoofe.ne.0) biotr(i,k,id_zoofe) = max(0.0,t_prog(ind_zoofe)%field(i,j,k,Time%taum1))
          if (id_detfe.ne.0) biotr(i,k,id_detfe) = max(0.0,t_prog(ind_detfe)%field(i,j,k,Time%taum1))
          !! IF OPTIONAL TRACER IS OFF, SET TO A CONSTANT VALUE 
          !if (id_caco3.eq.0) biotr(i,k,id_caco3) = 0.0 
          !if (id_po4.eq.0)   biotr(i,k,id_po4) = 1000.0    ! Make PO4 unlimiting to growth
          !if (id_fe.eq.0)    biotr(i,k,id_fe) = 1000.0     ! Make dFe unlimiting to growth
          !if (id_pchl.eq.0)  biotr(i,k,id_pchl) = 0.10
          !if (id_phyfe.eq.0) biotr(i,k,id_phyfe) = 10e-6
          !if (id_zoofe.eq.0) biotr(i,k,id_zoofe) = 10e-6
          !if (id_detfe.eq.0) biotr(i,k,id_detfe) = 10e-6

         enddo ! i
       enddo ! k


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [7] - Create tracer mask (currently based on where DIC == 0)                            !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

       biomask(:,:) = 0.0
       do k=1,grid%nk
         do i=isc,iec

           if (biotr(i,k,id_dic).gt.epsi) biomask(i,k) = 1.0

         enddo
       enddo


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [8] - Do light attenuation through water column & euphotic depth                        !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

       ! Reset some arrays
       ek_rgb(:,:,:) = 0.0
       par_rgb(:,:,:) = 0.0
       par_tot(:,:) = 0.0
       par_eup(:,:) = 0.0
       par_phy(:,:) = 0.0
       par_phymld(:,:) = 0.0
       par_tot_mldsum = 0.0
       par_phy_mldsum = 0.0
       par_z_mldsum = 0.0
       phy_lday(:,:) = 1e-2
       chl_lday(:,:) = 1e-2
       keuphot(:) = 4  ! set k level of euphotic zone to 10 metres
       daylight(:) = 0.0


       do i=isc,iec

         ! Calculate the daylength (hours per day) given the latitude and day of year
         day_angle = 2.0*pi*(days_in_this_year-1)/365.0
         declination = 0.409 * sin(2*pi*days_in_this_year/365 - 1.39)
         cos_hour_angle = max(min(-tan(declination) * tan(grid%yt(i,j) * pi/180.0), 1.0), -1.0)
         daylight(i) = 2 * (acos(cos_hour_angle) * 180 / pi) / 15.0

         ! Save the surface chlorophyll concentrations
         chlorophyll(i,j) = biotr(i,1,id_pchl)

       enddo

       do k=1,grid%nk
         do i=isc,iec

           ! Attenuation coefficients for Red, Green and Blue light given chlorophyll concentration
           zchl = MAX(0.05, MIN(10.0, biotr(i,k,id_pchl) ))
           ichl = NINT( 41 + 20.0*LOG10(zchl) + epsi )
           ek_rgb(i,k,1) = zrgb(2,ichl) * thickness%dzt(i,j,k)
           ek_rgb(i,k,2) = zrgb(3,ichl) * thickness%dzt(i,j,k)
           ek_rgb(i,k,3) = zrgb(4,ichl) * thickness%dzt(i,j,k)

           ! RGB light available in the water column given RGB attenuation
           if (k.eq.1) then
            if (swflx(i,j).gt.0.0) then
             par_rgb(i,k,1) = swflx(i,j) * parbio(i,j) * 1./3. * EXP(-0.5 * ek_rgb(i,k,1))
             par_rgb(i,k,2) = swflx(i,j) * parbio(i,j) * 1./3. * EXP(-0.5 * ek_rgb(i,k,2))
             par_rgb(i,k,3) = swflx(i,j) * parbio(i,j) * 1./3. * EXP(-0.5 * ek_rgb(i,k,3))
            else
             par_rgb(i,k,1) = 1./3. * EXP(-0.5 * ek_rgb(i,k,1))
             par_rgb(i,k,2) = 1./3. * EXP(-0.5 * ek_rgb(i,k,2))
             par_rgb(i,k,3) = 1./3. * EXP(-0.5 * ek_rgb(i,k,3))
            endif
           else
             par_rgb(i,k,1) = par_rgb(i,k-1,1) * EXP(-0.5 * (ek_rgb(i,k-1,1)+ek_rgb(i,k,1)))
             par_rgb(i,k,2) = par_rgb(i,k-1,2) * EXP(-0.5 * (ek_rgb(i,k-1,2)+ek_rgb(i,k,2)))
             par_rgb(i,k,3) = par_rgb(i,k-1,3) * EXP(-0.5 * (ek_rgb(i,k-1,3)+ek_rgb(i,k,3)))
           endif

           ! Light available for phytoplankton
           if (swflx(i,j).gt.0.0) then
             par_tot(i,k) =        par_rgb(i,k,1) +        par_rgb(i,k,2) +        par_rgb(i,k,3)
             par_phy(i,k) = 1.85 * par_rgb(i,k,1) + 0.68 * par_rgb(i,k,2) + 0.46 * par_rgb(i,k,3)
           else
             par_tot(i,k) = 0.0
             par_phy(i,k) = 0.0
           endif
           par_eup(i,k) = par_rgb(i,k,1) + par_rgb(i,k,2) + par_rgb(i,k,3)

           ! Integrated light level in the mixed layer
           if (grid%zw(k) .le. mld(i,j)) then
             par_tot_mldsum(i) = par_tot_mldsum(i) + par_tot(i,k) * thickness%dzt(i,j,k)
             par_phy_mldsum(i) = par_phy_mldsum(i) + par_phy(i,k) * thickness%dzt(i,j,k)
             par_z_mldsum(i) = par_z_mldsum(i) + thickness%dzt(i,j,k)
           endif

         enddo  ! i
       enddo  ! k


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [9] - Calculate depth of the euphotic zone                                              !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

       ! Calculate the euphotic zone depth, making sure it's not greater than the depth of the grid
       do i=isc,iec
         do k=1,grid%kmt(i,j)

           if (swflx(i,j).gt.0.0) then
             if (par_eup(i,k) .gt. (swflx(i,j)*parbio(i,j)*0.01)) then
               zeuphot(i,j) = grid%zw(k)
               keuphot(i) = k
             endif
           else
             if (par_eup(i,k) .gt. 0.01) then
               zeuphot(i,j) = grid%zw(k)
               keuphot(i) = k
             endif
           endif

         enddo  ! k
       enddo  ! i


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [10] - Aggregate light in mixed layer and calculate phytoplankton maximum growth rate   !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

       mumax_phy(:,:) = 0.0     ! phytoplankton maximum growth rate (temperature-dependent)

       do k=1,grid%nk
         do i=isc,iec

           ! Calculate average light level in the mixed layer
           if (grid%zw(k) .le. mld(i,j)) then
             zval = 1.0/(par_z_mldsum(i)+epsi)
             par_phymld(i,k) = par_phy_mldsum(i) * zval
           else
             par_phymld(i,k) = par_phy(i,k)
           endif

           ! Calculate impact of daylength on phytoplankton and chlorophyll production
           zval = MAX(1.0, daylight(i))
           if (grid%zw(k) .le. mld(i,j)) then
             zval = zval * MIN(1.0, (zeuphot(i,j)+epsi)/(mld(i,j)+epsi))
           endif
           chl_lday(i,k) = zval / 24.0
           phy_lday(i,k) = 1.5 * zval / (12.0 + zval)

           ! Save total PAR to radbio array for diagnostic output
           radbio3d(i,j,k) = par_tot(i,k)

           ! Temperature-dependent maximum growth rate (Eppley curve)
           ztem = t_prog(index_temp)%field(i,j,k,time%tau)
           mumax_phy(i,k) = abioa_phy(i,j) * bbioa(i,j)**(cbioa(i,j) * ztem) * 1.0/86400.0

         enddo  ! i
       enddo  ! k


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11] - Begin nested time-stepping for ecosystem model                                   !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

       ! Arrays for checking that ecosystem component equations conserve mass
       n_pools(:,:,:) = 0.0
       c_pools(:,:,:) = 0.0

       do tn = 1,ts_npzd

         ! Initialise some arrays to zero
         phy_chlc(:,:) = 0.0      ! chlorophyll to carbon ratio of phytoplankton
         phy_FeC(:,:) = 0.0       ! iron to carbon ratio of phytoplankton
         zoo_FeC(:,:) = 0.0       ! iron to carbon ratio of zooplankton
         det_FeC(:,:) = 0.0       ! iron to carbon ratio of slow sinking detritus
         phy_Feupt(:,:) = 0.0     ! iron uptake rate of phytoplankton
         phy_lpar(:,:) = 1.0      ! light limitation of phytoplankton growth
         mu_phy(:,:) = 0.0        ! phytoplankton realised growth rate (light limited)
         phy_chl_lpar(:,:) = 0.0  ! light limitation of chlorophyll production by phytoplankton
         mu_pchl(:,:) = 0.0       ! Growth of chlorophyll in phytoplankton
         phy_limnit(:,:) = 1.0    ! Phytoplankton nitrogen limitation
         phy_limfe(:,:) = 1.0     ! Phytoplankton iron limitation
         phy_limpho(:,:) = 1.0    ! Phytoplankton phosphorus limitation

         n_pools(:,:,2) = n_pools(:,:,1)  ! Save values from the last timestep
         c_pools(:,:,2) = c_pools(:,:,1)  ! Save values from the last timestep

         do k=1,grid%nk
           do i=isc,iec
         
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.1] - Update the tracer values within the nested timestep                            !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             biodic = max(0.0,biotr(i,k,id_dic))
             biono3 = max(0.0,biotr(i,k,id_no3))
             biophy = max(0.0,biotr(i,k,id_phy))
             biozoo = max(0.0,biotr(i,k,id_zoo))
             biodet = max(0.0,biotr(i,k,id_det))
             biooxy = max(0.0,biotr(i,k,id_o2))
             ! NOTE: for some reason the if statements here are essential to prevent model crash
             if (id_po4.ne.0)   biopo4 = max(0.0,biotr(i,k,id_po4))
             if (id_caco3.ne.0) biocaco3 = max(0.0,biotr(i,k,id_caco3))
             if (id_fe.ne.0)    biofe = max(0.0,biotr(i,k,id_fe))
             if (id_pchl.ne.0)  biopchl = max(0.0,biotr(i,k,id_pchl))
             if (id_phyfe.ne.0) biophyfe = max(0.0,biotr(i,k,id_phyfe))
             if (id_zoofe.ne.0) biozoofe = max(0.0,biotr(i,k,id_zoofe))
             if (id_detfe.ne.0) biodetfe = max(0.0,biotr(i,k,id_detfe))

             !if (k.lt.10 .and. i.eq.170 .and. j.eq.120) then
             !  print*, "DIC", biodic
             !  print*, "Nitrate", biono3
             !  print*, "Nanophy", biophy
             !  print*, "Microzoo", biozoo
             !  print*, "Detritus", biodet
             !  print*, "Oxygen", biooxy
             !  print*, "PO4", biopo4
             !  print*, "CaCO3", biocaco3
             !  print*, "iron (nM)", biofe
             !  print*, "Nanophy Chl", biopchl
             !  print*, "Nanophy Fe", biophyfe
             !  print*, "Microzoo Fe", biozoofe
             !  print*, "Detrital Fe", biodetfe
             !endif
   
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.2] - Get important ratios and values                                                !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             if (biophy*biopchl.gt.0.0)  phy_chlc(i,k) = biopchl / (biophy * 12)
             if (biophy*biophyfe.gt.0.0) phy_FeC(i,k)  = biophyfe/biophy
             if (biozoo*biozoofe.gt.0.0) zoo_FeC(i,k)  = biozoofe/biozoo
             if (biodet*biodetfe.gt.0.0) det_FeC(i,k)  = biodetfe/biodet
   
             ! Temperature and salinity
             ztem = t_prog(index_temp)%field(i,j,k,time%tau)
             ztemk = ztem + 273.15
             zsal = t_prog(index_salt)%field(i,j,k,time%tau)

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.3] - Dissolved Iron chemistry                                                       !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
             !!! THINGS NEEDING TO BE DONE:
             !!!   - Variable ligand concentration rather than constant?
             !!!   - Variable pH value (currently set to constant 7.9)
             !!!

             ! 1. Determine equilibrium fractionation of total dFe into Fe` and L-Fe
             ligand = 0.7  ! nM ligand concentration
             fe_keq = 10**( 17.27 - 1565.7 / ztemk )
             fe_III = ( -( 1. + fe_keq * ligand*1e-9 - fe_keq * biofe*1e-9 )                       &
                        + SQRT( ( 1. + fe_keq * ligand*1e-9 - fe_keq * biofe*1e-9 )**2             &
                                + 4. * biofe*1e-9 * fe_keq) ) / ( 2. * fe_keq + epsi ) * 1e9
             fe_lig = max(0.0, biofe - fe_III)

             ! 2. Precipitation of Fe` (creation of nanoparticles)
             zval = 19.924 * zsal / ( 1000. - 1.005 * zsal)
             fesol1 = 10**(-13.486 - 0.1856*zval**0.5 + 0.3073*zval + 5254.0/max(ztemk, 278.15) )
             fesol2 = 10**(2.517 - 0.8885*zval**0.5 + 0.2139*zval - 1320.0/max(ztemk, 278.15) )
             fesol3 = 10**(0.4511 - 0.3305*zval**0.5 - 1996.0/max(ztemk, 278.15) )
             fesol4 = 10**(-0.2965 - 0.7881*zval**0.5 - 4086.0/max(ztemk, 278.15) )
             fesol5 = 10**(4.4466 - 0.8505*zval**0.5 - 7980.0/max(ztemk, 278.15) )
             hp = 10**(-7.9)
             fe3sol = fesol1 * ( hp**3 + fesol2 * hp**2 + fesol3 * hp + fesol4 + fesol5 / hp ) *1e9
             precip = max(0.0, ( fe_III - fe3sol ) ) * knano_dfe(i,j)/86400.0

             ! 3. Scavenging of Fe` (Fe` --> DET)
             partic = (biodet + biocaco3)
             scaven = fe_III * (3e-8 + kscav_dfe(i,j) * partic) / 86400.0
             scadet = scaven * (biodet+epsi) / (partic+epsi)

             ! 4. Increased loss of Fe in near-coastal environments
             zval = min(1.0, 1000.0 / (grid%zt(grid%kmt(i,j))+1) )
             feloss = 1e-4 / 86400.0 * zval * biofe

             ! 5. Coagulation of colloidal Fe (nM) into small and large particles (uM)
             fe_col = fe_lig * 0.5
             if (grid%zw(k) .le. mld(i,j)) then
               zval =        ( 0.001 + biodet*kcoag_dfe(i,j) )
             else
               zval = 0.01 * ( 0.001 + biodet*kcoag_dfe(i,j) )
             endif
             fe2det = fe_col * zval / 86400.0


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.4] - Nutrient limitation terms for phytoplankton growth                             !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             ! Calculate the variable half-saturation Ks (change with population biomass)
             phy_conc1 = max(0.0, biophy - phybiot(i,j))  ! phybiot is threshold of phyto biomass in µM
             phy_conc2 = biophy - phy_conc1
             phy_k_nit = max( phykn(i,j),                                                          &
                              ( phy_conc2 * phykn(i,j) + phy_conc1 * 3.0 * phykn(i,j) )            &
                               / ( biophy + epsi ) )
             phy_k_fer = max( phykf(i,j),                                                          &
                              ( phy_conc2 * phykf(i,j) + phy_conc1 * 3.0 * phykf(i,j) )            &
                              / ( biophy + epsi ) )
             phy_k_pho = max( phykp(i,j),                                                          &
                              ( phy_conc2 * phykp(i,j) + phy_conc1 * 3.0 * phykp(i,j) )            &
                              / ( biophy + epsi ) )


             ! Nitrogen (NO3)
             phy_limnit(i,k) = biono3 / ( biono3 + phy_k_nit + epsi )


             ! Iron
             if (biopchl*biophy.gt.0.0) then
               phy_minFeC = 0.0016 / 55.85 * biopchl/biophy +                                      &
                            1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * phy_limnit(i,k) +         &
                            1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * phy_limnit(i,k)
             else
               phy_minFeC = 0.0016 / 55.85 * phyminqc(i,j)*12 +                                    &
                            1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * phy_limnit(i,k) +         &
                            1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * phy_limnit(i,k)
             endif
             phy_limfe(i,k) = min(1.0, max(0.0, (phy_FeC(i,k) - phy_minFeC)/phyoptqf(i,j) ))


             ! Phosphorus (PO4)
             if (id_po4.ne.0) then
               phy_limpho(i,k) = biopo4 / ( biopo4 + phy_k_pho + epsi )
             else
               phy_limpho(i,k) = 1.0
             endif


             !if (k.lt.10 .and. i.eq.170 .and. j.eq.120) then
             !  print*, "Nutrient limitation at lon =", grid%xt(i,j)
             !  print*, "Nutrient limitation at lat =", grid%yt(i,j)
             !  print*, "Nutrient limitation at depth =", grid%zt(k)
             !  print*, "phy_k_nit", phy_k_nit
             !  print*, "phy_k_fer", phy_k_fer
             !  print*, "phy_k_pho", phy_k_pho
             !  print*, "phy_limnit(i,k)", phy_limnit(i,k)
             !  print*, "phy_limfe(i,k)", phy_limfe(i,k)
             !  print*, "phy_limpho(i,k)", phy_limpho(i,k)
             !  print*, "phy_FeC(i,k)", phy_FeC(i,k)
             !  print*, "phy_minFeC", phy_minFeC
             !  print*, "biopchl/biophy", biopchl/biophy
             !  print*, "phyoptqf(i,j)", phyoptqf(i,j)
             !  print*, " "
             !endif

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.5] - Light limitation terms for phytoplankton growth                                !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             ! 1. initial slope of Photosynthesis-Irradiance curve
             ! 2. Alter the slope to account for respiration and daylength limitation
             ! 3. Light limitation
             ! 4. Apply to maximum growth rate

             phy_pisl(i,k) = MAX(alphabio(i,j) * phy_chlc(i,k), alphabio(i,j) * phyminqc(i,j))
             phy_pisl2= phy_pisl(i,k) / ( (1.0 + phylmor(i,j)*86400.0) * max(daylight(i),1.0)/24.0 + epsi )
             phy_lpar(i,k) = (1.0 - EXP(-phy_pisl2 * par_phy(i,k))) * phy_lday(i,k)
             mu_phy(i,k) = mumax_phy(i,k) * phy_lpar(i,k)

             !if (k.lt.10 .and. i.eq.170 .and. j.eq.120) then
             !  print*, "Light limitation at depth =", grid%zt(k)
             !  print*, "par_phy(i,k)", par_phy(i,k)
             !  print*, "phy_chlc(i,k)", phy_chlc(i,k)
             !  print*, "denominator", ( (1.0 + phylmor(i,j)*86400.0) * max(daylight(i),1.0)/24.0 + epsi )
             !  print*, "phy_pisl(i,k)", phy_pisl(i,k)
             !  print*, "phy_pisl2", phy_pisl2
             !  print*, "phy_lpar(i,k)", phy_lpar(i,k)
             !  print*, "mumax_phy(i,k)", mumax_phy(i,k)
             !  print*, "mu_phy(i,k)", mu_phy(i,k)
             !  print*, " "
             !endif


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.6] - Realised growth rate given resource availability                               !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             phy_limnut = min( phy_limnit(i,k), phy_limpho(i,k), phy_limfe(i,k) )
             u_phy = phy_limnut * mu_phy(i,k)


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.7] - Growth in chlorophyll (mg Chl / m3)                                            !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             !     1. Light-limitation of chlorophyll production        (phy_chl_lpar)
             !     2. Minimum and optimal rates of chlorophyll growth   (mumin_pchl & muopt_pchl)
             !     3. Calculate mg Chl m-3 s-1                          (mu_pchl)

             phy_chl_pisl = phy_pisl(i,k) / ( mumax_phy(i,k) * 86400.0                             &
                                              * max(daylight(i),1.0)/24.0                          &   
                                              * (1.0 - phy_limnut) + epsi)
             phy_chl_lpar(i,k) = (1.0 - EXP( -phy_chl_pisl * par_phymld(i,k))) * phy_lday(i,k)
             mumin_pchl = phyminqc(i,j) * u_phy * biophy * 12.
             muopt_pchl = phyoptqc(i,j) * u_phy * biophy * 12.
             mu_pchl(i,k) = (muopt_pchl - mumin_pchl) * phy_chl_lpar(i,k) * phy_limnut
             if (phy_pisl(i,k)*par_phymld(i,k).gt.0.0) then
               mu_pchl(i,k) = mumin_pchl + mu_pchl(i,k) / (phy_pisl(i,k) * par_phymld(i,k))
             endif


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.8] - Phytoplankton uptake of Fe (nmol Fe / m3)                                     !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             ! 1. Maximum iron content of phytoplankton biomass
             ! 2. Ensure that Fe uptake increases or decreases in response to standing quotas
             ! 3. Iron uptake by phytoplankton

             if (phy_FeC(i,k).gt.0.0 .and. biophy.gt.0.0) then
               phy_maxQFe = biophy * phymaxqf(i,j)
               phy_Feupt_upreg   = (4.0 - 4.5 * phy_limfe(i,k) / (phy_limfe(i,k) + 0.5))
               phy_Feupt_downreg = max(0.0, (1.0 - biophyfe/phy_maxQFe) / abs(1.05 - biophyfe/phy_maxQFe) )
               phy_Feupt(i,k) = biophy * mumax_phy(i,k) * phymaxqf(i,j) *                           &
                                biofe / (biofe + phy_k_fer + epsi) * phy_Feupt_downreg * phy_Feupt_upreg
             endif

             !if (k.lt.10 .and. i.eq.170 .and. j.eq.120) then
             !  print*, "Phyto Fe uptake at depth =", grid%zt(k) 
             !  print*, "phy_maxQFe", phy_maxQFe
             !  print*, "phy_Feupt_upreg", phy_Feupt_upreg
             !  print*, "phy_Feupt_downreg", phy_Feupt_downreg
             !  print*, "phy_Feupt(i,k)", phy_Feupt(i,k)
             !  print*, " "
             !endif


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.9] - Temperature dependence of heterotrophy / respiration                          !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             fbc = abioh(i,j) * bbioh(i,j)**(cbioh(i,j) * ztem)


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.10] - Grazing functions                                                             !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             biophy_min = MAX(0.0, biophy - 0.001)  ! Minimum Phy beneath which no grazing
             biodet_min = MAX(0.0, biodet - 0.001)  ! Minimum Det beneath which no grazing

             ! preferences and prey availability
             zooprey = zprefphy(i,j)*biophy_min + zprefdet(i,j)*biodet_min
             g_zoo = zoogmax(i,j)*fbc * epszoo(i,j) * zooprey*zooprey /                            & 
                     ( zoogmax(i,j)*fbc + epszoo(i,j) * zooprey*zooprey + epsi )


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.11] - Mortality scalers                                                             !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             ! Depress mortality in nutrient limited regime (i.e., investment in cell maintainence)
             phy_slmor = MIN(1.0, MAX(0.0,                                                         &
                         min(phy_limnit(i,k), phy_limpho(i,k), phy_limfe(i,k))/0.3))


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.12] - Collect major source and sink terms for tendency equations                    !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             ! PRIMARY PRODUCTION
             f11phy = u_phy * biophy                                        * biomask(i,k)         ! DIC --> PHY
             

             ! GRAZING
             if (zooprey.gt.0.0) then
               z21phy = g_zoo * biozoo * (zprefphy(i,j)*biophy_min/zooprey) * biomask(i,k)         ! PHY --> ZOO & DET & POC
               z21det = g_zoo * biozoo * (zprefdet(i,j)*biodet_min/zooprey) * biomask(i,k)         ! DET --> ZOO & DET
             else
               z21phy = 0.0
               z21det = 0.0
             endif

             ! MORTALITY
             f22phy = phylmor(i,j) * fbc * phy_slmor * biophy               * biomask(i,k)         ! PHY --> DIC
             f23phy = phyqmor(i,j) * biophy * biophy                        * biomask(i,k)         ! PHY --> DET
             z31 = zooexcr(i,j) * biozoo * fbc                              * biomask(i,k)         ! ZOO --> DIC
             z32 = zooqmor(i,j) * biozoo*biozoo                             * biomask(i,k)         ! ZOO --> DET

             ! REMINERALISATION
             f41det = detlrem(i,j) * biodet * fbc                           * biomask(i,k)         ! DET --> DIC
             if (grid%zw(k) .ge. 180) then
               f41det = f41det * 0.5   ! reduce decay below 180m
             endif
             f51 = caco3rem(i,j) * biocaco3                                 * biomask(i,k)


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.13] - Apply sources and sinks to tendency equations                                 !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             !!!!!!!!!!!!!!!!!!!
             ! oxygen equation !
             !!!!!!!!!!!!!!!!!!!
             biotr(i,k,id_o2)   = biotr(i,k,id_o2)  - biomask(i,k) * 172./122. * dtsb * (          &
                                                      f41det + z31 + f22phy - f11phy )

             !!!!!!!!!!!!!!!!!!!!!!
             ! phosphate equation !
             !!!!!!!!!!!!!!!!!!!!!!
             if (id_po4.ne.0)                                                                      &
             biotr(i,k,id_po4)  = biotr(i,k,id_po4) + dtsb * 1./122. * (                           &
                                                      f41det + z31 + f22phy - f11phy )

             !!!!!!!!!!!!!!!!!!!!
             ! nitrate equation !
             !!!!!!!!!!!!!!!!!!!!
             biotr(i,k,id_no3)  = biotr(i,k,id_no3) + dtsb * 16./122. * (                          &
                                                             f41det + z31 + f22phy - f11phy ) 


             !!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! dissolved iron equation !
             !!!!!!!!!!!!!!!!!!!!!!!!!!!
             if (id_fe.ne.0)                                                                       &
             biotr(i,k,id_fe)   = biotr(i,k,id_fe)  + dtsb * (f41det * det_FeC(i,k) +              &
                                                              z31 * zoo_FeC(i,k) +                 &
                                                              f22phy * phy_FeC(i,k) -              &
                                                              phy_Feupt(i,k)) * 1e3                &
                                                    - dtsb * (precip + scaven + feloss + fe2det)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! nanophytoplankton equation !
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             biotr(i,k,id_phy)  = biotr(i,k,id_phy) + dtsb * (f11phy - z21phy - f22phy - f23phy)
             if (id_phyfe.ne.0)                                                                    &
             biotr(i,k,id_phyfe)= biotr(i,k,id_phyfe)+dtsb * (phy_Feupt(i,k) -                     &
                                                              (z21phy+f22phy+f23phy)*phy_FeC(i,k))

             !if (k.lt.10 .and. i.eq.170 .and. j.eq.120) then
             !  print*, "Phyto sources and sinks terms =", grid%zt(k) 
             !  print*, "f11phy", f11phy
             !  print*, "z21phy", z21phy 
             !  print*, "f22phy", f22phy
             !  print*, "f23phy", f23phy
             !  print*, "phy_Feupt(i,k)", phy_Feupt(i,k)
             !  print*, "biomask(i,k)", biomask(i,k)
             !  print*, "biotr(i,k,id_dic)", biotr(i,k,id_dic)
             !  print*, " "
             !endif

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! microzooplankton equation !
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             biotr(i,k,id_zoo)  = biotr(i,k,id_zoo) + dtsb * (zooassi(i,j) *                       &
                                                             (z21phy + z21det) - z31 - z32)
             if (id_zoofe.ne.0)                                                                    &
             biotr(i,k,id_zoofe)= biotr(i,k,id_zoofe)+dtsb * (zooassi(i,j)*(z21phy*phy_FeC(i,k) +  &
                                                                            z21det*det_FeC(i,k)) - &
                                                              (z31+z32)*zoo_FeC(i,k) )


             phyxsize(i,j,k) = ( min(phybiot(i,j),biophy) + 3*max(0.0,biophy-phybiot(i,j)) )       &
                                / (biophy + epsi)

             !!!!!!!!!!!!
             ! Detritus !
             !!!!!!!!!!!!
             biotr(i,k,id_det)  = biotr(i,k,id_det) + dtsb * ( (1.-zooassi(i,j)) *(z21det+z21phy)  &
                                                              + f23phy + z32 - z21det - f41det ) 
             if (id_detfe.ne.0)                                                                    &
             biotr(i,k,id_detfe)= biotr(i,k,id_detfe)+dtsb * ((1.-zooassi(i,j))*z21det*det_FeC(i,k)&
                                                             +(1.-zooassi(i,j))*z21phy*phy_FeC(i,k)&
                                                             + f23phy * phy_FeC(i,k)               &
                                                             + z32 * zoo_FeC(i,k)                  &
                                                             - (z21det + f41det)*det_FeC(i,k)      &
                                                             + (scadet + fe2det)*1e-3 )

             !!!!!!!!!!!!!!!!!!
             ! CaCO3 equation !
             !!!!!!!!!!!!!!!!!!
             if (id_caco3.ne.0)                                                                    &
             biotr(i,k,id_caco3) = biotr(i,k,id_caco3) + dtsb * ( (                                &
                                   (1.-zooassi(i,j))*z21phy + f23phy + z32 ) * f_inorg(i,j) - f51)
  

             !!!!!!!!!!!!!!!!
             ! DIC equation !
             !!!!!!!!!!!!!!!!
             biotr(i,k,id_dic)  = biotr(i,k,id_dic) + dtsb * ( f41det + f51                        &
                                                              + z31 + f22phy                       & 
                                                              - f11phy                             &
                                      - ( (1.-zooassi(i,j))*z21phy + f23phy + z32 )*f_inorg(i,j) )

             !!!!!!!!!!!!!!!!!!!!!!!!!
             ! Chlorophyll equations !
             !!!!!!!!!!!!!!!!!!!!!!!!!
             if (id_pchl.ne.0)                                                                     &
             biotr(i,k,id_pchl) = biotr(i,k,id_pchl) + dtsb * ( mu_pchl(i,k) - phy_chlc(i,k) *     &
                                                                (z21phy + f22phy + f23phy) * 12.0 )
  

            !print*, "Oxygen concentration", biotr(i,k,id_o2)
            !print*, "Phosphate concentration", biotr(i,k,id_po4)
            !print*, "Nitrate concentration", biotr(i,k,id_no3)
            !print*, "dFe concentration", biotr(i,k,id_fe)
            !print*, "Nano concentration", biotr(i,k,id_phy)
            !print*, "Nano Fe concentration", biotr(i,k,id_phyfe)
            !print*, "Zoo concentration", biotr(i,k,id_zoo)
            !print*, "Zoo Fe concentration", biotr(i,k,id_zoofe)
            !print*, "Det concentration", biotr(i,k,id_det)
            !print*, "Det Fe concentration", biotr(i,k,id_detfe)
            !print*, "CaCO3 concentration", biotr(i,k,id_caco3)
            !print*, "pChl concentration", biotr(i,k,id_pchl)

   
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.14] - Collect key diagnostics                                                       !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             ! Gross NPP
             pprod_gross(i,j,k) = pprod_gross(i,j,k) + dtsb*(f11phy)
   
             ! Net NPP (gross minus linear mortality (respiration))
             npp3d(i,j,k) = npp3d(i,j,k) + dtsb*(f11phy - f22phy)
   
             ! Secondary productivity
             zprod_gross(i,j,k) = zprod_gross(i,j,k) + dtsb*zooassi(i,j)*(z21phy+z21det)
   
             ! Net Secondary production (gross minus respiration losses)
             nsp3d(i,j,k) = nsp3d(i,j,k) + dtsb*( (z21phy + z21det)*zooassi(i,j) - z31 )

             ! Light limitation of phytoplankton
             phy_parlimit(i,j,k) = phy_parlimit(i,j,k) + phy_lpar(i,k)*dtsb/dtts
             ! Iron limitation of phytoplankton
             phy_Felimit(i,j,k) = phy_Felimit(i,j,k) + phy_limfe(i,k)*dtsb/dtts
             ! Nitrogen limitation of phytoplankton
             phy_Nlimit(i,j,k) = phy_Nlimit(i,j,k) + phy_limnit(i,k)*dtsb/dtts
             ! Phosphorus limitation of phytoplankton
             phy_Plimit(i,j,k) = phy_Plimit(i,j,k) + phy_limpho(i,k)*dtsb/dtts
   
             ! Specific grazing pressure of zooplankton (µM Z per µM P per second)
             if (zooprey.gt.0.0) then
               zoo_grazpres(i,j,k) = zoo_grazpres(i,j,k) + dtsb * (z21phy+z21det) / (biophy+biodet)
             endif


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.15] - Check for conservation of mass by ecosystem component                         !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             n_pools(i,k,1) = biotr(i,k,id_no3) + (biotr(i,k,id_phy) + biotr(i,k,id_det)           &
                              + biotr(i,k,id_zoo))*16/122.
             c_pools(i,k,1) = biotr(i,k,id_dic) + biotr(i,k,id_phy) + biotr(i,k,id_det)            &
                              + biotr(i,k,id_zoo) + biotr(i,k,id_caco3)
             
             if (tn.gt.1) then
               if (abs(n_pools(i,k,2) - n_pools(i,k,1)) .gt. 1e-6 ) then
                 print*, "Error: Ecosystem model is not conserving nitrogen"
                 print*, "       Longitude = ", grid%xt(i,j)
                 print*, "       Latitude  = ", grid%yt(i,j)
                 print*, "       Depth     = ", grid%zt(k)
                 print*, " "
                 print*, "       Biological N budget (mmol N m-3) at two timesteps ", n_pools(i,k,:)
                 stop
               endif
               if (abs(c_pools(i,k,2) - c_pools(i,k,1)) .gt. 1e-6 ) then
                 print*, "Error: Ecosystem model is not conserving carbon"
                 print*, "       Longitude = ", grid%xt(i,j)
                 print*, "       Latitude  = ", grid%yt(i,j)
                 print*, "       Depth     = ", grid%zt(k)
                 print*, " "
                 print*, "       Biological C budget (mmol C m-3) at two timesteps ", c_pools(i,k,:)
                 stop
               endif
             endif
             

             !if (i.eq.170 .and. j.eq.120 .and. k.lt.10) then
             !  print*, "Longitude =", grid%xt(i,j)
             !  print*, "Latitude =", grid%yt(i,j)
             !  print*, "Depth =", grid%zt(k)
             !  print*, "timestep", tn, "of",ts_npzd
             !  print*, "biotr(i,k,id_phy)*16/122.", biotr(i,k,id_phy)*16/122.
             !  print*, "biotr(i,k,id_zoo)*16/122.", biotr(i,k,id_zoo)*16/122.
             !  print*, "biotr(i,k,id_det)*16/122.", biotr(i,k,id_det)*16/122.
             !  print*, "biotr(i,k,id_no3)", biotr(i,k,id_no3)
             !  print*, "Bio Nitrogen budget", n_pools(i,k,:)
             !  print*, "Bio Carbon budget", c_pools(i,k,:)
             !endif

           enddo  ! i
         enddo  ! k
       enddo  ! tn


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [12] - Add biotically-induced tendency to biotracers                                    !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

       do k = 1,grid%nk
         do i = isc,iec

           ! Collect total change in NO3, but add back NO3 if negative at tau-1
           no3_bgc_change = grid%tmask(i,j,k) *                                                    &
                            (biotr(i,k,id_no3) - max(0.0,t_prog(ind_no3)%field(i,j,k,time%taum1)))

           ! Collect total change in CaCO3, but add back CaCO3 if negative at tau-1
           caco3_bgc_change = grid%tmask(i,j,k) *                                                  &
                              (biotr(i,k,id_caco3) - max(0.0,t_prog(ind_caco3)%field(i,j,k,time%taum1)))


           ! Add sources and sinks of BGC tracers to tracer arrays
           t_prog(ind_o2)%field(i,j,k,time%taum1)  = biotr(i,k,id_o2)
           t_prog(ind_no3)%field(i,j,k,time%taum1) = biotr(i,k,id_no3)
           t_prog(ind_phy)%field(i,j,k,time%taum1) = biotr(i,k,id_phy)
           t_prog(ind_zoo)%field(i,j,k,time%taum1) = biotr(i,k,id_zoo)
           t_prog(ind_det)%field(i,j,k,time%taum1) = biotr(i,k,id_det)
           t_prog(ind_dic)%field(i,j,k,time%taum1) = biotr(i,k,id_dic) 
           if ( id_po4.ne.0) t_prog(ind_po4)%field(i,j,k,time%taum1) = biotr(i,k,id_po4)
           if ( id_caco3.ne.0) t_prog(ind_caco3)%field(i,j,k,time%taum1) = biotr(i,k,id_caco3)
           if ( id_fe.ne.0) t_prog(ind_fe)%field(i,j,k,time%taum1) = biotr(i,k,id_fe)
           if (id_adic.ne.0) &
           t_prog(ind_adic)%field(i,j,k,time%taum1) = t_prog(ind_adic)%field(i,j,k,time%taum1) +   &
                                                      122./16. * no3_bgc_change - caco3_bgc_change
           if (id_alk.ne.0) &
           t_prog(ind_alk)%field(i,j,k,time%taum1) = t_prog(ind_alk)%field(i,j,k,time%taum1) +     &
                                                     ( -2.0 * caco3_bgc_change - no3_bgc_change )
           if (id_pchl.ne.0)  t_prog(ind_pchl)%field(i,j,k,time%taum1) = biotr(i,k,id_pchl)
           if (id_phyfe.ne.0) t_prog(ind_phyfe)%field(i,j,k,time%taum1) = biotr(i,k,id_phyfe)
           if (id_zoofe.ne.0) t_prog(ind_zoofe)%field(i,j,k,time%taum1) = biotr(i,k,id_zoofe)
           if (id_detfe.ne.0) t_prog(ind_detfe)%field(i,j,k,time%taum1) = biotr(i,k,id_detfe)



!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [13] - Collect other diagnostic terms outside of the nested ecosystem time-step         !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

           ! Collect primary and secondary production diagnostic terms
           pprod_gross(i,j,k) = rdtts * pprod_gross(i,j,k) * grid%tmask(i,j,k)
           zprod_gross(i,j,k) = rdtts * zprod_gross(i,j,k) * grid%tmask(i,j,k)
           npp3d(i,j,k)       = rdtts * npp3d(i,j,k)       * grid%tmask(i,j,k)
           nsp3d(i,j,k)       = rdtts * nsp3d(i,j,k)       * grid%tmask(i,j,k)
           zoo_grazpres(i,j,k)= rdtts * zoo_grazpres(i,j,k)* grid%tmask(i,j,k)
  
           ! Collect other diagnostic terms
           if (Grid%zw(k) .le. mld(i,j)) then
              adic_intmld(i,j)= adic_intmld(i,j)+ t_prog(ind_adic)%field(i,j,k,time%taum1)* thickness%dzt(i,j,k)
              dic_intmld(i,j) = dic_intmld(i,j) + t_prog(ind_dic)%field(i,j,k,time%taum1) * thickness%dzt(i,j,k)
              o2_intmld(i,j)  = o2_intmld(i,j)  + t_prog(ind_o2)%field(i,j,k,time%taum1)  * thickness%dzt(i,j,k)
              no3_intmld(i,j) = no3_intmld(i,j) + t_prog(ind_no3)%field(i,j,k,time%taum1) * thickness%dzt(i,j,k)
              fe_intmld(i,j)  = fe_intmld(i,j)  + t_prog(ind_fe)%field(i,j,k,time%taum1)  * thickness%dzt(i,j,k)
              phy_intmld(i,j) = phy_intmld(i,j) + t_prog(ind_phy)%field(i,j,k,time%taum1) * thickness%dzt(i,j,k)
              det_intmld(i,j) = det_intmld(i,j) + t_prog(ind_det)%field(i,j,k,time%taum1) * thickness%dzt(i,j,k)
              pprod_gross_intmld(i,j) = pprod_gross_intmld(i,j) + pprod_gross(i,j,k) * thickness%dzt(i,j,k)
              npp_intmld(i,j) = npp_intmld(i,j) + npp3d(i,j,k) * thickness%dzt(i,j,k)
              radbio_intmld(i,j) = radbio_intmld(i,j) + radbio3d(i,j,k) * thickness%dzt(i,j,k)
              light_limit(i,j) = light_limit(i,j) + phy_lpar(i,k) * thickness%dzt(i,j,k)
           endif
           if (Grid%zw(k) .le. 100) then
              adic_int100(i,j)= adic_int100(i,j)+ t_prog(ind_adic)%field(i,j,k,time%taum1)* thickness%dzt(i,j,k)
              dic_int100(i,j) = dic_int100(i,j) + t_prog(ind_dic)%field(i,j,k,time%taum1) * thickness%dzt(i,j,k)
              o2_int100(i,j)  = o2_int100(i,j)  + t_prog(ind_o2)%field(i,j,k,time%taum1)  * thickness%dzt(i,j,k)
              no3_int100(i,j) = no3_int100(i,j) + t_prog(ind_no3)%field(i,j,k,time%taum1) * thickness%dzt(i,j,k)
              fe_int100(i,j)  = fe_int100(i,j)  + t_prog(ind_fe)%field(i,j,k,time%taum1)  * thickness%dzt(i,j,k)
              phy_int100(i,j) = phy_int100(i,j) + t_prog(ind_phy)%field(i,j,k,time%taum1) * thickness%dzt(i,j,k)
              det_int100(i,j) = det_int100(i,j) + t_prog(ind_det)%field(i,j,k,time%taum1) * thickness%dzt(i,j,k)
              pprod_gross_int100(i,j) = pprod_gross_int100(i,j) + pprod_gross(i,j,k) * thickness%dzt(i,j,k)
              npp_int100(i,j) = npp_int100(i,j) + npp3d(i,j,k) * thickness%dzt(i,j,k)
              radbio_int100(i,j) = radbio_int100(i,j) + radbio3d(i,j,k) * thickness%dzt(i,j,k)
           endif

         enddo  ! i
       enddo  ! k


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [14] - Sink detritus, detrital iron, detrital silicate and CaCO3                        !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

       ! collect tracer values before sinking and find sinking rates
       do i=isc,iec

         ! Find sinking rate as a function of average phytoplankton community size in euphotic zone
         wdet(i,j,:) = wdetbio(i,j) * max(1.0, min(3.0,                                            &
                       SUM(phyxsize(i,j,1:keuphot(i)) * Thickness%dzt(i,j,1:keuphot(i)), dim=1)    &
                       / SUM(Thickness%dzt(i,j,1:keuphot(i)), dim=1) ))

         do k=1,grid%nk

           ! Collect tracer values before we sink them 
           biotr(i,k,id_det) = max(0.0,t_prog(ind_det)%field(i,j,k,Time%taum1))
           if (id_detfe.ne.0) biotr(i,k,id_detfe) = max(0.0,t_prog(ind_detfe)%field(i,j,k,Time%taum1))

           ! Sinking rate increases logarithmically with pressure (also decreased very shallow)
           wdet(i,j,k) = wdet(i,j,k) * max(0.1, min(2.0, 0.5*LOG10(grid%zt(k)/4 + 1.0) + 0.1))

           !if (i.eq.170 .and. j.eq.120) then
           !  print*, " "
           !  print*, "Longitude =", grid%xt(i,j)
           !  print*, "Latitude =", grid%yt(i,j)
           !  print*, "Depth =", grid%zt(k)
           !  print*, "Sinking rate (m/s)", wdet(i,j,k) * 86400.0
           !endif

         enddo
       enddo
      
       ! Initialise deposition arrays
       biodetsed(:) = 0.0
       biodetfesed(:) = 0.0

       ! start timesplitting
       do tn = 1,ts_npzd !{

         ! initialise advective sinking arrays
         adv_fb1(:,:) = 0.0
         adv_fb2(:,:) = 0.0


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! RASF upstream sinking of detritus and particulates !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do i=isc,iec
           do k=2,grid%nk+1
             adv_fb1(i,k) = wdet(i,j,k-1) * max(0.0,biotr(i,k-1,id_det))
             if (id_detfe.ne.0) adv_fb2(i,k) = wdet(i,j,k-1) * max(0.0,biotr(i,k-1,id_detfe))
           enddo  ! i
         enddo  ! k

         ! Deposit tracer to sediment as tracer sinks through base of column,
         !    accounting for burial and denitrification in csiro_bgc.F90
         do i = isc,iec
           k = grid%kmt(i,j)
           if (k.gt.0) then
             biodetsed(i) = biodetsed(i) + adv_fb1(i,k+1)
             if (id_detfe.ne.0) biodetfesed(i) = biodetfesed(i) + adv_fb2(i,k+1)
           endif  ! k.gt.0
         enddo  ! i

         ! Update the tracer arrays at each timesplit
         do k =1,grid%nk  !{
           do i =isc,iec  !{
             biotr(i,k,id_det) = biotr(i,k,id_det) + grid%tmask(i,j,k) * dtsb *                    &
                                 (-adv_fb1(i,k+1) + adv_fb1(i,k))/Thickness%dzt(i,j,k)
             if (id_detfe.ne.0)                                                                    &
             biotr(i,k,id_detfe) = biotr(i,k,id_detfe) + grid%tmask(i,j,k) * dtsb *                &
                                 (-adv_fb2(i,k+1) + adv_fb2(i,k))/Thickness%dzt(i,j,k)
           enddo  !} i
         enddo  !} k

       enddo !} tn

       ! Update the actual tracer values and deposit material to sediments after the timesplitting
       do i = isc,iec !{
         do k = 1,grid%nk !{
           t_prog(ind_det)%field(i,j,k,time%taum1) = biotr(i,k,id_det)
           if (id_detfe.ne.0) t_prog(ind_detfe)%field(i,j,k,time%taum1) = biotr(i,k,id_detfe)
         enddo !}
         biotic(n)%det_sed_depst(i,j) = biodetsed(i)
         biotic(n)%detfe_sed_depst(i,j) = biodetfesed(i)
       enddo !}

       ! Make sure we pass this zero array of sediment silicate deposition to subroutine "csiro_bgc_bbc"
       biotic(n)%detsi_sed_depst(:,:) = 0.0

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! RASF upstream sinking of caco3 !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       adv_fb1(:,:) = 0.0

       if (id_caco3.ne.0) then

         do k=2,grid%nk+1
           do i=isc,iec
             adv_fb1(i,k) = wcaco3(i,j) * biotr(i,k-1,id_caco3)
           enddo
         enddo

         ! Deposit tracer to sediment as tracer sinks through base of column.  mac, nov12
         do i = isc, iec
           k = grid%kmt(i,j)
           if (k .gt. 0) then
             biotic(n)%caco3_sed_depst(i,j) = adv_fb1(i,k+1)
           endif ! k .gt. 0
         enddo ! i

         do k =1,grid%nk  !{
           do i =isc,iec  !{
             t_prog(ind_caco3)%field(i,j,k,time%taum1)=t_prog(ind_caco3)%field(i,j,k,time%taum1) + &
                 grid%tmask(i,j,k) * dtts * (-adv_fb1(i,k+1) + adv_fb1(i,k))/Thickness%dzt(i,j,k)
           enddo  ! i
         enddo  ! k

       endif  ! end loop for caco3

 enddo  !} j

enddo  !} n

return
end subroutine bio_v3
