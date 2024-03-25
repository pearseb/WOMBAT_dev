!

subroutine bio_v4(isc, iec, jsc, jec, isd, ied, jsd, jed, T_prog, Grid, Time, dtts, Thickness, Dens, swflx, sw_frac_zt, &
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
!       11.7  Find rate of nitrogen fixation (facultative)
!       11.8  Si:C ratio for diatom productivity 
!       11.9  Silicate dissolution rate 
!       11.10 Growth in chlorophyll (mg Chl / m3)
!       11.11 Phytoplankton uptake of Fe (nmol Fe / m3)
!       11.12 Temperature dependence of heterotrophy / respiration
!       11.13 Grazing functions 
!       11.14 Mortality scalers
!       11.15 Bacterial processes
!       11.16 Collect major source and sink terms for tendency equations
!       11.17 Apply sources and sinks to tendency equations 
!       11.18 Collect key diagnostics
!       11.19 Check for conservation of mass by ecosystem component
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
real :: bionh4                          ! MAX(0.0, ammonium)
real :: biopo4                          ! MAX(0.0, phosphate) 
real :: biono3                          ! MAX(0.0, nitrate)
real :: biosil                          ! MAX(0.0, silicic acid)
real :: biofe                           ! MAX(0.0, dissolved iron)
real :: biooxy                          ! MAX(0.0, oxygen)
real :: biophy                          ! MAX(0.0, nanophytoplankton)
real :: biodia                          ! MAX(0.0, microphytoplankton)
real :: biodiz                          ! MAX(0.0, diazotrophs)
real :: biozoo                          ! MAX(0.0, microzooplankton)
real :: biomes                          ! MAX(0.0, mesozooplankton)
real :: biodet                          ! MAX(0.0, slow-sinking detritus)
real :: biopoc                          ! MAX(0.0, fast-sinking detrius)
real :: biocaco3                        ! MAX(0.0, calcium carbonate)
real :: biopchl                         ! MAX(0.0, nanophytoplankton chlorophyll) 
real :: biodchl                         ! MAX(0.0, microphytoplankton chlorophyll) 
real :: biozchl                         ! MAX(0.0, diazotroph chlorophyll) 
real :: biophyfe                        ! MAX(0.0, nanophytoplankton iron)
real :: biodiafe                        ! MAX(0.0, microphytoplankton iron)
real :: biodizfe                        ! MAX(0.0, diazotroph iron)
real :: biozoofe                        ! MAX(0.0, microzooplankton iron)
real :: biomesfe                        ! MAX(0.0, mesozooplankton iron)
real :: biodiasi                        ! MAX(0.0, microphytoplankton silicate)
real :: biodetfe                        ! MAX(0.0, slow-sinking detrital iron)
real :: biopocfe                        ! MAX(0.0, fast-sinking detrital iron)
real :: biopocsi                        ! MAX(0.0, fast-sinking detrital silicate)
real :: biophy_min                      ! MAX(some_threshold, nanophytoplankton)     || beneath some_threshold no mortality occurs
real :: biodia_min                      ! MAX(some_threshold, microphytoplankton     || beneath some_threshold no mortality occurs
real :: biodiz_min                      ! MAX(some_threshold, diazotroph             || beneath some_threshold no mortality occurs
real :: biozoo_min                      ! MAX(some_threshold, mesozooplankton)       || beneath some_threshold no mortality occurs
real :: biodet_min                      ! MAX(some_threshold, slow-sinking detritus) || beneath some_threshold no mortality occurs
real :: biopoc_min                      ! MAX(some_threshold, fast-sinking detritus) || beneath some_threshold no mortality occurs
real :: biomask(isc:iec,grid%nk)        ! Land-sea masking
real :: biotr(isc:iec,grid%nk,ntr_bmax) ! Biogeochemical tracer array

! --------------------------------------------- !
!!! ----------------------------------------- !!!
!!! Limitation terms for phytoplankton growth !!!
!!! ----------------------------------------- !!!
! --------------------------------------------- !
real :: phy_limnut                      ! Nanophytoplankton limitation by all nutrients   || 0 == full limitation
real :: phy_limnh4                      ! Nanophytoplankton limitation by ammonium        || 0 == full limitation
real :: phy_limno3                      ! Nanophytoplankton limitation by nitrate         || 0 == full limitation
real :: dia_limnut                      ! Microphytoplankton limitation by all nutrients  || 0 == full limitation
real :: dia_limnh4                      ! Microphytoplankton limitation by ammonium       || 0 == full limitation
real :: dia_limno3                      ! Mircophytoplankton limitation by nitrate        || 0 == full limitation
real :: diz_limnut                      ! Diazotroph limitation by all nutrients  || 0 == full limitation
real :: diz_limnh4                      ! Diazotroph limitation by ammonium       || 0 == full limitation
real :: diz_limno3                      ! Diazotroph limitation by nitrate        || 0 == full limitation
real :: phy_zlimnit                     ! Nanophytoplankton initial guess at nitrogen limitation 
real :: phy_zlimnh4                     ! Nanophytoplankton initial guess at ammonium limitation
real :: phy_zlimno3                     ! Nanophytoplankton initial guess at nitrate limitation
real :: dia_zlimnit                     ! Microphytoplankton initial guess at nitrogen limitation
real :: dia_zlimnh4                     ! Microphytoplankton initial guess at ammonium limitation
real :: dia_zlimno3                     ! Microphytoplankton initial guess at nitrate limitation
real :: diz_zlimnit                     ! Diazotroph initial guess at nitrogen limitation
real :: diz_zlimnh4                     ! Diazotroph initial guess at ammonium limitation
real :: diz_zlimno3                     ! Diazotroph initial guess at nitrate limitation
real :: phy_conc1                       ! Nanophytoplankton biomass exceeding bloom threshold
real :: phy_conc2                       ! Nanophytoplankton biomass (but maximum == bloom threshold)
real :: dia_conc1                       ! Microphytoplankton biomass exceeding bloom threshold
real :: dia_conc2                       ! Microphytoplankton biomass (but maximum == bloom threshold) 
real :: diz_conc1                       ! Diazotroph biomass exceeding bloom threshold
real :: diz_conc2                       ! Diazotroph biomass (but maximum == bloom threshold) 
real :: phy_k_nit                       ! Nanophytoplankton half-saturation coefficient for dissolved inorganic nitrogen uptake
real :: phy_k_fer                       ! Nanophytoplankton half-saturation coefficient for dissolved iron uptake
real :: phy_k_pho                       ! Nanophytoplankton half-saturation coefficient for phosphate uptake
real :: dia_k_nit                       ! Microphytoplankton half-saturation coefficient for dissolved inorganic nitrogen uptake
real :: dia_k_fer                       ! Microphytoplankton half-saturation coefficient for dissolved iron uptake
real :: dia_k_pho                       ! Microphytoplankton half-saturation coefficient for phosphate uptake
real :: dia_k_sil                       ! Microphytoplankton half-saturation coefficient for slicic acid uptake 
real :: diz_k_nit                       ! Diazotroph half-saturation coefficient for dissolved inorganic nitrogen uptake
real :: diz_k_fer                       ! Diazotroph half-saturation coefficient for dissolved iron uptake
real :: diz_k_pho                       ! Diazotroph half-saturation coefficient for phosphate uptake
real :: phy_limnit(isc:iec,grid%nk)     ! Nanophytoplankton limitation by nitrogen
real :: phy_limpho(isc:iec,grid%nk)     ! Nanophytoplankton limitation by phosphorus
real :: phy_limfe(isc:iec,grid%nk)      ! Nanophytoplankton limitation by iron
real :: dia_limnit(isc:iec,grid%nk)     ! Microphytoplankton limitation by nitrogen
real :: dia_limpho(isc:iec,grid%nk)     ! Microphytoplankton limitation by phosphorus
real :: dia_limsil(isc:iec,grid%nk)     ! Microphytoplankton limitation by silicic acid
real :: dia_limfe(isc:iec,grid%nk)      ! Microphytoplankton limitation by iron
real :: diz_limnit(isc:iec,grid%nk)     ! Diazotroph limitation by nitrogen
real :: diz_limpho(isc:iec,grid%nk)     ! Diazotroph limitation by phosphorus
real :: diz_limfe(isc:iec,grid%nk)      ! Diazotroph limitation by iron

! ------------------------------------- !
!!! --------------------------------- !!!
!!! Mortality terms for phytoplankton !!!
!!! --------------------------------- !!!
! ------------------------------------- !
real :: phy_slmor                       ! Scaler [0-1] to depress linear nanophytoplankton mortality in nutrient-limited regimes
real :: dia_slmor                       ! Scaler [0-1] to depress linear microphytoplankton mortality in nutrient-limited regimes
real :: diz_slmor                       ! Scaler [0-1] to depress linear diazotroph mortality in nutrient-limited regimes

! ---------------------------- !
!!! ------------------------ !!!
!!! Growth and grazing rates !!!
!!! ------------------------ !!!
! ---------------------------- !
real :: u_phy                           ! Realised growth rate of nanophytoplankton (nutrients & light limitations)
real :: u_dia                           ! Realised growth rate of microphytoplankton (nutrients & light limitations)
real :: u_diz                           ! Realised growth rate of microphytoplankton (nutrients & light limitations)
real :: g_zoo                           ! Realised grazing rate of microzooplankton
real :: g_mes                           ! Realised grazing rate of mesozooplankton
real :: zooprey                         ! Total prey biomass for microzooplankton
real :: mesprey                         ! Total prey biomass for mesozooplankton
real :: fbc                             ! Heterotrophic temperature-dependence function
real :: mumax_phy(isc:iec,grid%nk)      ! Maximum potential (T-dependent) growth rate of nanophytoplankton
real :: mumax_dia(isc:iec,grid%nk)      ! Maximum potential (T-dependent) growth rate of microphytoplankton
real :: mumax_diz(isc:iec,grid%nk)      ! Maximum potential (T-dependent) growth rate of microphytoplankton
real :: mu_phy(isc:iec,grid%nk)         ! Realised growth rate of nanophytoplankton (light limitation)
real :: mu_dia(isc:iec,grid%nk)         ! Realised growth rate of microphytoplankton (light limitation)
real :: mu_diz(isc:iec,grid%nk)         ! Realised growth rate of microphytoplankton (light limitation)

! --------------------------------------- !
!!! ----------------------------------- !!!
!!! Tendency terms for tracer equations !!!
!!! ----------------------------------- !!!
! --------------------------------------- !
real :: f11phy                          ! Nanophytoplankton biomass growth 
real :: f11dia                          ! Microphytoplankton biomass growth
real :: f11diz                          ! Diazotroph biomass growth
real :: z21phy                          ! Microzooplankton grazing of nanophytoplankton
real :: z21dia                          ! Microzooplankton grazing of microphytoplankton
real :: z21diz                          ! Microzooplankton grazing of diazotroph
real :: z21det                          ! Microzooplankton grazing of slow-sinking detritus 
real :: z21poc                          ! Microzooplankton grazing of fast-sinking detritus 
real :: m21phy                          ! Mesozooplankton grazing of nanophytoplankton 
real :: m21dia                          ! Mesozooplankton grazing of microphytoplankton 
real :: m21diz                          ! Mesozooplankton grazing of diazotroph 
real :: m21det                          ! Mesozooplankton grazing of slow-sinking detritus 
real :: m21poc                          ! Mesozooplankton grazing of fast-sinking detritus 
real :: m21zoo                          ! Mesozooplankton grazing of microzooplankton 
real :: f22phy                          ! Nanophytoplankton biomass loss via linear mortality (respiration)
real :: f23phy                          ! Nanophytoplankton biomass loss via quadratic mortality (viral lysis, etc.)
real :: f22dia                          ! Microphytoplankton biomass loss via linear mortality (respiration)
real :: f23dia                          ! Microphytoplankton biomass loss via quadratic mortality (viral lysis, etc.)
real :: f22diz                          ! Diazotroph biomass loss via linear mortality (respiration)
real :: f23diz                          ! Diazotroph biomass loss via quadratic mortality (viral lysis, etc.)
real :: z31                             ! Microzooplankton biomass loss via linear mortality (excretion)
real :: z32                             ! Mesozooplankton biomass loss via quadratic mortality (predation, infection, etc.)
real :: m31                             ! Mesozooplankton biomass loss via linear mortality (excretion)
real :: m32                             ! Mesozooplankton biomass loss via quadratic mortality (predation, infection, etc.)
real :: f41det                          ! Remineralisation of slow-sinking detritus to inorganic constituents   
real :: f41poc                          ! Remineralisation of fast-sinking detritus to inorganic constituents   
real :: f41pocsi                        ! Remineralisation of fast-sinking detrital silicate to silicic acid
real :: f51                             ! Remineralisation of CaCO3 to DIC and Alkalinity
real :: caco3_bgc_change                ! Total change in CaCO3 following time-splitting of ecosystem component 
real :: no3_bgc_change                  ! Total change in nitrate following time-splitting of ecosystem component 
real :: nh4_bgc_change                  ! Total change in ammonium following time-splitting of ecosystem component
real :: c_pools(isc:iec,grid%nk,2)      ! Total biological carbon budget (PHY + DIA + ZOO + MES + DET + POC + DIC + CaCO3)
real :: n_pools(isc:iec,grid%nk,2)      ! Total biological nitrogen budget (PHY + DIA + ZOO + MES + DET + POC + NO3 + NH4) 

! ------------------------- !
!!! --------------------- !!!
!!! Light and chlorophyll !!!
!!! --------------------- !!!
! ------------------------- !
real :: zchl                            ! Total chlorophyll in mg/m3  || used to find light attenuation coefficients
real :: phy_pisl2                       ! Photosynthetic-Irradiance slope of nanophytoplankton (~ Chl:C, respiration, daylength, zeuph/zmld)
real :: dia_pisl2                       ! Photosynthetic-Irradiance slope of microphytoplankton (~ Chl:C, respiration, daylength, zeuph/zmld)
real :: diz_pisl2                       ! Photosynthetic-Irradiznce slope of diazotroph (~ Chl:C, respiration, daylength, zeuph/zmld)
real :: phy_chl_pisl                    ! Photosynthetic-Irradiance slope of nanophytoplankton chlorophyll (~ Chl:C, µmax, daylength, zeuph/zmld) 
real :: dia_chl_pisl                    ! Photosynthetic-Irradiance slope of microphytoplankton chlorophyll (~ Chl:C, µmax, daylength, zeuph/zmld) 
real :: diz_chl_pisl                    ! Photosynthetic-Irradiznce slope of diazotroph chlorophyll (~ Chl:C, µmax, daylength, zeuph/zmld) 
real :: mumin_pchl                      ! Rate of chlorophyll growth in nanophytoplankton to achieve minimum Chl:C ratio
real :: muopt_pchl                      ! Rate of chlorophyll growth in nanophytoplankton to achieve optimal Chl:C ratio
real :: mumin_dchl                      ! Rate of chlorophyll growth in microphytoplankton to achieve minimum Chl:C ratiol 
real :: muopt_dchl                      ! Rate of chlorophyll growth in microphytoplankton to achieve optimal Chl:C ratiol 
real :: mumin_zchl                      ! Rate of chlorophyll growth in diazotroph to achieve minimum Chl:C ratiol 
real :: muopt_zchl                      ! Rate of chlorophyll growth in diazotroph to achieve optimal Chl:C ratiol 
real :: keuphot(isc:iec)                ! Estimated depth of the euphotic zone
real :: daylight(isc:iec)               ! Number of daylight hours in a day
real :: par_phy_mldsum(isc:iec)         ! Sum of photosynthetic active radiation seen by nanophytoplankton in the mixed layer
real :: par_dia_mldsum(isc:iec)         ! Sum of photosynthetic active radiation seen by microphytoplankton in the mixed layer
real :: par_diz_mldsum(isc:iec)         ! Sum of photosynthetic active radiztion seen by diazotroph in the mixed layer
real :: par_tot_mldsum(isc:iec)         ! Sum of photosynthetic active radiation in the mixed layer
real :: par_z_mldsum(isc:iec)           ! Sum of depths in the mixed layer
real :: par_tot(isc:iec,grid%nk)        ! Total photosynthetic active radiation (red + green + blue)
real :: par_phy(isc:iec,grid%nk)        ! Photosynthetic active radiation seen by nanophytoplankton (red + green + blue)
real :: par_dia(isc:iec,grid%nk)        ! Photosynthetic active radiation seen by microphytoplankton (red + green + blue)
real :: par_diz(isc:iec,grid%nk)        ! Photosynthetic active radiztion seen by diazotroph (red + green + blue)
real :: par_phymld(isc:iec,grid%nk)     ! Average photosynthetic active radiation seen by nanophytoplankton (red + green + blue) in mixed layer 
real :: par_diamld(isc:iec,grid%nk)     ! Average photosynthetic active radiation seen by microphytoplankton (red + green + blue) in mixed layer
real :: par_dizmld(isc:iec,grid%nk)     ! Average photosynthetic active radiztion seen by diazotroph (red + green + blue) in mixed layer
real :: par_eup(isc:iec,grid%nk)        ! Photosynthetic active radiation for calculation of euphotic zone depth
real :: phy_pisl(isc:iec,grid%nk)       ! Initial guess at photosynthesis-irradiance slope (~ Chl:C) for nanophytoplankton
real :: dia_pisl(isc:iec,grid%nk)       ! Initial guess at photosynthesis-irradiance slope (~ Chl:C) for microphytoplankton
real :: diz_pisl(isc:iec,grid%nk)       ! Initial guess at photosynthesis-irradiznce slope (~ Chl:C) for diazotroph
real :: phy_lday(isc:iec,grid%nk)       ! Daylength and euphotic < MLD limitation term for nanophytoplankton 
real :: dia_lday(isc:iec,grid%nk)       ! Daylength and euphotic < MLD limitation term for microphytoplankton
real :: diz_lday(isc:iec,grid%nk)       ! Daylength and euphotic < MLD limitation term for diazotroph
real :: chl_lday(isc:iec,grid%nk)       ! Daylength and euphotic < MLD limitation term for chlorophyll production
real :: phy_lpar(isc:iec,grid%nk)       ! Light limitation term for nanophytoplankton
real :: dia_lpar(isc:iec,grid%nk)       ! Light limitation term for microphytoplankton
real :: diz_lpar(isc:iec,grid%nk)       ! Light limitation term for diazotroph
real :: phy_chl_lpar(isc:iec,grid%nk)   ! Light limitation term for nanophytoplankton chlorophyll
real :: dia_chl_lpar(isc:iec,grid%nk)   ! Light limitation term for microphytoplankton chlorophyll
real :: diz_chl_lpar(isc:iec,grid%nk)   ! Light limitation term for diazotroph chlorophyll
real :: phy_chlc(isc:iec,grid%nk)       ! Nanophytoplankton Chl:C ratio
real :: dia_chlc(isc:iec,grid%nk)       ! Microphytoplankton Chl:C ratio
real :: diz_chlc(isc:iec,grid%nk)       ! Diazotroph Chl:C ratio
real :: mu_pchl(isc:iec,grid%nk)        ! Growth rate of nanophytoplankton chlorophyll
real :: mu_dchl(isc:iec,grid%nk)        ! Growth rate of microphytoplankton chlorophyll
real :: mu_zchl(isc:iec,grid%nk)        ! Growth rate of diazotroph chlorophyll
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
real :: dia_minFeC                      ! Minimum Fe:C quota of microphytoplankton given their Chlorophyll content and N limitations         
real :: dia_maxQFe                      ! Maximum Fe:C quota of microphytoplankton
real :: dia_Feupt_upreg                 ! Upregulation of Fe uptake by microphytoplankton in response to low Fe cellular quota
real :: dia_Feupt_downreg               ! Downregulation of Fe uptake by microphytoplankton in response to achieving high cellular Fe quotas
real :: diz_minFeC                      ! Minimum Fe:C quota of diazotroph given their Chlorophyll content and N limitations         
real :: diz_maxQFe                      ! Maximum Fe:C quota of diazotroph
real :: diz_Feupt_upreg                 ! Upregulation of Fe uptake by diazotroph in response to low Fe cellular quota
real :: diz_Feupt_downreg               ! Downregulation of Fe uptake by diazotroph in response to achieving high cellular Fe quotas
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
real :: scaven                          ! Rate of scavenging of Fe' (Fe(III)) to total detritus
real :: scadet                          ! Rate of scavenging of Fe' (Fe(III)) to slow-sinking fraction of detritus
real :: scapoc                          ! Rate of scavenging of Fe' (Fe(III)) to fast-sinking fraction of detritus
real :: feloss                          ! Rate of ambiguous dFe losses in near-coastal environments (i.e., amplified scavenging, etc.)
real :: fe2det                          ! Rate of coagulation of colloidal dFe to slow-sinking detritus
real :: fe2poc                          ! Rate of coagulation of colloidal dFe to fast-sinking detritus
real :: phy_FeC(isc:iec,grid%nk)        ! Nanophytoplankton Fe:C ratio
real :: dia_FeC(isc:iec,grid%nk)        ! Microphytoplankton Fe:C ratio
real :: diz_FeC(isc:iec,grid%nk)        ! Diazotroph Fe:C ratio
real :: zoo_FeC(isc:iec,grid%nk)        ! Microzooplankton Fe:C ratio
real :: mes_FeC(isc:iec,grid%nk)        ! Mesozooplankton Fe:C ratio
real :: det_FeC(isc:iec,grid%nk)        ! Slow-sinking detritus Fe:C ratio
real :: poc_FeC(isc:iec,grid%nk)        ! Fast-sinking detritus Fe:C ratio
real :: phy_Feupt(isc:iec,grid%nk)      ! Nanophytoplankton Fe uptake rate
real :: dia_Feupt(isc:iec,grid%nk)      ! Microphytoplankton Fe uptake rate
real :: diz_Feupt(isc:iec,grid%nk)      ! Diazotroph Fe uptake rate

! ----------------- !
!!! ------------- !!!
!!! Silicon cycle !!!
!!! ------------- !!!
! ----------------- !
real :: sicfac1                         ! Factor determining Si:C ratio of diatoms --> varies from [1,5.4] as silicate increases and as nutrient/daylight limitation increases 
real :: sicfac2                         ! Factor determining Si:C ratio of diatoms --> varies from [1,3] as silicate increases
real :: sio3eq                          ! Solubility of silicic acid in seawater. Increases as temperature increases
real :: sisat1                          ! Silicic acid saturation factor --> goes to 1 when silicic acid concentration is much less than the solubility
real :: sisat2                          ! Silicic acid saturation factor --> 0.8 at 0ºC and 14.5 at 30ºC
real :: sidisr                          ! Overall temperature-dependent scaler for the dissolution rate of biogenic silicate 
real :: siremin                         ! Rate of biogenic silicate dissolution into silicic acid
real :: dia_SiC(isc:iec,grid%nk)        ! Microphytoplankton Si:C ratio
real :: poc_SiC(isc:iec,grid%nk)        ! Fast-sinking detritus Si:C ratio
real :: si_xlabile(isc:iec,grid%nk)     ! Fraction of biogenic silicate that is labile
real :: si_lab2non(isc:iec,grid%nk)     ! Ratio of labile to non-labile biogenic silicate
real :: dia_SiCupt(isc:iec,grid%nk)     ! Microphytoplankton Si:C uptake rate

! -------------------- !
!!! ---------------- !!!
!!! Detrital routing !!!
!!! ---------------- !!!
! -------------------- !
real :: phy_det2poc                     ! Fraction of nanophytoplankton biomass that is routed to slow-sinking detritus
real :: dia_det2poc                     ! Fraction of microphytoplankton biomass that is routed to slow-sinking detritus
real :: diz_det2poc                     ! Fraction of diazotroph biomass that is routed to slow-sinking detritus
real :: zoo_det2poc                     ! Fraction of microzooplankton biomass that is routed to slow-sinking detritus
real :: mes_det2poc                     ! Fraction of mesozooplankton biomass that is routed to slow-sinking detritus

! -------------------- !
!!! ---------------- !!!
!!! Detrital sinking !!!
!!! ---------------- !!!
! -------------------- !
real :: adv_fb1(isc:iec,1:grid%nk+1)    ! Sinking of slow-sinking detritus (carbon)
real :: adv_fb2(isc:iec,1:grid%nk+1)    ! Sinking of slow-sinking detritus (iron)
real :: adv_fb3(isc:iec,1:grid%nk+1)    ! Sinking of fast-sinking detritus (carbon)
real :: adv_fb4(isc:iec,1:grid%nk+1)    ! Sinking of fast-sinking detritus (iron)
real :: adv_fb5(isc:iec,1:grid%nk+1)    ! Sinking of fast-sinking detritus (silicate)

! ----------------------- !
!!! ------------------- !!!
!!! Bacterial processes !!!
!!! ------------------- !!!
! ----------------------- !
real :: nitrif(isc:iec,grid%nk)         ! Nitrification rate (NH4 --> NO3)
real :: denitr(isc:iec,grid%nk)         ! Denitrification rate (NO3 --> N2)



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
!

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
  call time_interp_external(diabiot_id, time%model_time, diabiot)   !Microphytoplankton biomass threshold
  call time_interp_external(dizbiot_id, time%model_time, dizbiot)   !Diazotroph biomass threshold
  call time_interp_external(phyminqc_id, time%model_time, phyminqc) !Phytoplankton minimum Chl:C quota
  call time_interp_external(diaminqc_id, time%model_time, diaminqc) !Microphytoplankton minimum Chl:C quota 
  call time_interp_external(dizminqc_id, time%model_time, dizminqc) !Diazotroph minimum Chl:C quota 
  call time_interp_external(phyoptqc_id, time%model_time, phyoptqc) !Phytoplankton maximum Chl:C quota
  call time_interp_external(diaoptqc_id, time%model_time, diaoptqc) !Microphytoplankton maximum Chl:C quota 
  call time_interp_external(dizoptqc_id, time%model_time, dizoptqc) !Diazotroph maximum Chl:C quota 
  call time_interp_external(phymaxqf_id, time%model_time, phymaxqf) !Phytoplankton maximum Fe:C quota
  call time_interp_external(diamaxqf_id, time%model_time, diamaxqf) !Microphytoplankton maximum Fe:C quota 
  call time_interp_external(dizmaxqf_id, time%model_time, dizmaxqf) !Diazotroph maximum Fe:C quota 
  call time_interp_external(phyoptqf_id, time%model_time, phyoptqf) !Phytoplankton optimal Fe:C quota
  call time_interp_external(diaoptqf_id, time%model_time, diaoptqf) !Microphytoplankton optimal Fe:C quota 
  call time_interp_external(dizoptqf_id, time%model_time, dizoptqf) !Diazotroph optimal Fe:C quota 
  call time_interp_external(abioa_phy_id, time%model_time, abioa_phy) !Autotrophic T-dependent scalar
  call time_interp_external(abioa_dia_id, time%model_time, abioa_dia) !Autotrophic T-dependent scalar
  call time_interp_external(abioa_diz_id, time%model_time, abioa_diz) !Autotrophic T-dependent scalar (maximum growth rate)
  call time_interp_external(bbioa_id, time%model_time, bbioa)       !Autotrophic T-dependent base
  call time_interp_external(beta_diz_id, time%model_time, beta_diz) !Temperature at which maximum growth occurs for diazotrophs
  call time_interp_external(cbioa_id, time%model_time, cbioa)       !Autotrophic T-dependent exponent
  call time_interp_external(abioh_id, time%model_time, abioh)       !Heterotrophic T-dependent scalar
  call time_interp_external(bbioh_id, time%model_time, bbioh)       !Heterotrophic T-dependent base
  call time_interp_external(cbioh_id, time%model_time, cbioh)       !Heterotrophic T-dependent exponent
  call time_interp_external(phykn_id, time%model_time, phykn)       !Phytoplankton N half-saturation
  call time_interp_external(phykf_id, time%model_time, phykf)       !Phytoplankton Fe half-saturation
  call time_interp_external(phykp_id, time%model_time, phykp)       !Phytoplankton P half-saturation
  call time_interp_external(diakn_id, time%model_time, diakn)       !Microphytoplankton N half-saturation
  call time_interp_external(diakf_id, time%model_time, diakf)       !Microphytoplankton Fe half-saturation
  call time_interp_external(diakp_id, time%model_time, diakp)       !Microphytoplankton P half-saturation
  call time_interp_external(dizkn_id, time%model_time, dizkn)       !Diazotroph N half-saturation
  call time_interp_external(dizkf_id, time%model_time, dizkf)       !Diazotroph Fe half-saturation
  call time_interp_external(dizkp_id, time%model_time, dizkp)       !Diazotroph P half-saturation
  call time_interp_external(phylmor_id, time%model_time, phylmor)   !Phytoplankton linear mortality
  call time_interp_external(phyqmor_id, time%model_time, phyqmor)   !Phytoplankton quadratic mortality
  call time_interp_external(zooassi_id, time%model_time, zooassi)   !Zooplankton assimilation efficiency
  call time_interp_external(zoogmax_id, time%model_time, zoogmax)   !Zooplankton max grazing rate
  call time_interp_external(mesgmax_id, time%model_time, mesgmax)   !Mesozooplankton max grazing rate
  call time_interp_external(epszoo_id, time%model_time, epszoo)     !Zooplankton prey capture rate
  call time_interp_external(epsmes_id, time%model_time, epsmes)     !Mesozooplankton prey capture rate
  call time_interp_external(zprefphy_id, time%model_time, zprefphy) !Zooplankton preference for phytoplankton
  call time_interp_external(zprefdia_id, time%model_time, zprefdia) !Zooplankton preference for microphytoplankton
  call time_interp_external(zprefdiz_id, time%model_time, zprefdiz) !Zooplankton preference for diazotroph
  call time_interp_external(zprefdet_id, time%model_time, zprefdet) !Zooplankton preference for Det
  call time_interp_external(zprefpoc_id, time%model_time, zprefpoc) !Zooplankton preference for POC
  call time_interp_external(mprefphy_id, time%model_time, mprefphy) !Mesozooplankton preference for phytoplankton
  call time_interp_external(mprefdia_id, time%model_time, mprefdia) !Mesozooplankton preference for microphytoplankton
  call time_interp_external(mprefdiz_id, time%model_time, mprefdiz) !Mesozooplankton preference for diazotroph
  call time_interp_external(mprefdet_id, time%model_time, mprefdet) !Mesozooplankton preference for Det
  call time_interp_external(mprefpoc_id, time%model_time, mprefpoc) !Mesozooplankton preference for POC
  call time_interp_external(mprefzoo_id, time%model_time, mprefzoo) !Mesozooplankton preference for Zooplankton
  call time_interp_external(zooqmor_id, time%model_time, zooqmor)   !Zooplankton quadratic mortality
  call time_interp_external(mesqmor_id, time%model_time, mesqmor)   !Mesozooplankton quadratic mortality
  call time_interp_external(zooexcr_id, time%model_time, zooexcr)   !Zooplankton excretion rate
  call time_interp_external(detlrem_id, time%model_time, detlrem)   !Detritus remineralisation rate
  call time_interp_external(wdetbio_id, time%model_time, wdetbio)   !Sinking rate of slow detritus
  call time_interp_external(caco3rem_id, time%model_time, caco3rem) !CaCO3 remineralisation rate
  call time_interp_external(wcaco3_id, time%model_time, wcaco3)     !CaCO3 sinking rate
  call time_interp_external(tscav_fe_id, time%model_time, tscav_fe) !Iron scavenging rate
  call time_interp_external(fe_bkgnd_id, time%model_time, fe_bkgnd) !Background dFe concentration
  call time_interp_external(f_inorg_id, time%model_time, f_inorg)   !CaCO3 inorganic fraction


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
  mprod_gross(:,:,:) = 0.0
  light_limit(:,:) = 0.0
  radbio3d(:,:,:) = 0.0
  phyxsize(:,:,:) = 0.0
  diaxsize(:,:,:) = 0.0
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
  phy_parlimit(:,:,:) = 0.0
  dia_parlimit(:,:,:) = 0.0
  diz_parlimit(:,:,:) = 0.0
  phy_Felimit(:,:,:) = 0.0
  dia_Felimit(:,:,:) = 0.0
  diz_Felimit(:,:,:) = 0.0
  phy_Nlimit(:,:,:) = 0.0
  dia_Nlimit(:,:,:) = 0.0
  diz_Nlimit(:,:,:) = 0.0
  phy_Plimit(:,:,:) = 0.0
  dia_Plimit(:,:,:) = 0.0
  diz_Plimit(:,:,:) = 0.0
  dia_Silimit(:,:,:) = 0.0
  dia_SiCupta(:,:,:) = 0.0
  zoo_grazpres(:,:,:) = 0.0
  mes_grazpres(:,:,:) = 0.0
  nitrific(:,:,:) = 0.0
  denitrif(:,:,:) = 0.0
  diazofix(:,:,:) = 0.0


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
  ind_nh4  = biotic(n)%ind_bgc(id_nh4)      ! mmol N / m3
  ind_no3  = biotic(n)%ind_bgc(id_no3)      ! mmol N / m3
  ind_sil  = biotic(n)%ind_bgc(id_sil)      ! mmol Si / m3
  ind_fe   = biotic(n)%ind_bgc(id_fe)       ! µmol Fe / m3
  ind_o2   = biotic(n)%ind_bgc(id_o2)       ! mmol O2 / m3
  ind_dic  = biotic(n)%ind_bgc(id_dic)      ! mmol DIC / m3
  ind_alk  = biotic(n)%ind_bgc(id_alk)      ! mmol Eq / m3
  ind_caco3= biotic(n)%ind_bgc(id_caco3)    ! mmol C / m3
  ind_adic = biotic(n)%ind_bgc(id_adic)     ! mmol DIC / m3
  ind_phy  = biotic(n)%ind_bgc(id_phy)      ! mmol C / m3
  ind_dia  = biotic(n)%ind_bgc(id_dia)      ! mmol C / m3
  ind_diz= biotic(n)%ind_bgc(id_diz)        ! mmol C / m3
  ind_zoo  = biotic(n)%ind_bgc(id_zoo)      ! mmol C / m3
  ind_mes  = biotic(n)%ind_bgc(id_mes)      ! mmol C / m3
  ind_det  = biotic(n)%ind_bgc(id_det)      ! mmol C / m3
  ind_poc  = biotic(n)%ind_bgc(id_poc)      ! mmol C / m3
  ind_pchl = biotic(n)%ind_bgc(id_pchl)     ! mg / m3
  ind_dchl = biotic(n)%ind_bgc(id_dchl)     ! mg / m3
  ind_zchl = biotic(n)%ind_bgc(id_zchl)     ! mg / m3
  ind_phyfe = biotic(n)%ind_bgc(id_phyfe)   ! mmol Fe / m3
  ind_diafe = biotic(n)%ind_bgc(id_diafe)   ! mmol Fe / m3
  ind_dizfe = biotic(n)%ind_bgc(id_dizfe)   ! mmol Fe / m3
  ind_diasi = biotic(n)%ind_bgc(id_diasi)   ! mmol Si / m3
  ind_zoofe = biotic(n)%ind_bgc(id_zoofe)   ! mmol Fe / m3
  ind_mesfe = biotic(n)%ind_bgc(id_mesfe)   ! mmol Fe / m3
  ind_detfe = biotic(n)%ind_bgc(id_detfe)   ! mmol Fe / m3
  ind_pocfe = biotic(n)%ind_bgc(id_pocfe)   ! mmol Fe / m3
  ind_pocsi = biotic(n)%ind_bgc(id_pocsi)   ! mmol Si / m3


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

           biotr(i,k,id_dic) = max(0.0,t_prog(ind_dic)%field(i,j,k,Time%taum1)) 
           biotr(i,k,id_no3) = max(0.0,t_prog(ind_no3)%field(i,j,k,Time%taum1)) 
           biotr(i,k,id_phy) = max(0.0,t_prog(ind_phy)%field(i,j,k,Time%taum1))
           biotr(i,k,id_zoo) = max(0.0,t_prog(ind_zoo)%field(i,j,k,Time%taum1))
           biotr(i,k,id_det) = max(0.0,t_prog(ind_det)%field(i,j,k,Time%taum1))
           biotr(i,k,id_o2)  = max(0.0,t_prog(ind_o2)%field(i,j,k,Time%taum1)) 
           if (id_nh4.ne.0)   biotr(i,k,id_nh4) = max(0.0,t_prog(ind_nh4)%field(i,j,k,Time%taum1))
           if (id_po4.ne.0)   biotr(i,k,id_po4) = max(0.0,t_prog(ind_po4)%field(i,j,k,Time%taum1))
           if (id_sil.ne.0)   biotr(i,k,id_sil) = max(0.0,t_prog(ind_sil)%field(i,j,k,Time%taum1))
           if (id_dia.ne.0)   biotr(i,k,id_dia) = max(0.0,t_prog(ind_dia)%field(i,j,k,Time%taum1))
           if (id_diz.ne.0) biotr(i,k,id_diz) = max(0.0,t_prog(ind_diz)%field(i,j,k,Time%taum1))
           if (id_mes.ne.0)   biotr(i,k,id_mes) = max(0.0,t_prog(ind_mes)%field(i,j,k,Time%taum1))
           if (id_poc.ne.0)   biotr(i,k,id_poc) = max(0.0,t_prog(ind_poc)%field(i,j,k,Time%taum1))
           if (id_caco3.ne.0) biotr(i,k,id_caco3) = max(0.0,t_prog(ind_caco3)%field(i,j,k,Time%taum1)) 
           if (id_fe.ne.0)    biotr(i,k,id_fe) = max(0.0,t_prog(ind_fe)%field(i,j,k,Time%taum1)) 
           if (id_pchl.ne.0)  biotr(i,k,id_pchl) = max(0.0,t_prog(ind_pchl)%field(i,j,k,Time%taum1))
           if (id_dchl.ne.0)  biotr(i,k,id_dchl) = max(0.0,t_prog(ind_dchl)%field(i,j,k,Time%taum1))
           if (id_zchl.ne.0)  biotr(i,k,id_zchl) = max(0.0,t_prog(ind_zchl)%field(i,j,k,Time%taum1))
           if (id_phyfe.ne.0) biotr(i,k,id_phyfe) = max(0.0,t_prog(ind_phyfe)%field(i,j,k,Time%taum1))
           if (id_diafe.ne.0) biotr(i,k,id_diafe) = max(0.0,t_prog(ind_diafe)%field(i,j,k,Time%taum1))
           if (id_dizfe.ne.0) biotr(i,k,id_dizfe) = max(0.0,t_prog(ind_dizfe)%field(i,j,k,Time%taum1))
           if (id_diasi.ne.0) biotr(i,k,id_diasi) = max(0.0,t_prog(ind_diasi)%field(i,j,k,Time%taum1))
           if (id_zoofe.ne.0) biotr(i,k,id_zoofe) = max(0.0,t_prog(ind_zoofe)%field(i,j,k,Time%taum1))
           if (id_mesfe.ne.0) biotr(i,k,id_mesfe) = max(0.0,t_prog(ind_mesfe)%field(i,j,k,Time%taum1))
           if (id_detfe.ne.0) biotr(i,k,id_detfe) = max(0.0,t_prog(ind_detfe)%field(i,j,k,Time%taum1))
           if (id_pocfe.ne.0) biotr(i,k,id_pocfe) = max(0.0,t_prog(ind_pocfe)%field(i,j,k,Time%taum1))
           if (id_pocsi.ne.0) biotr(i,k,id_pocsi) = max(0.0,t_prog(ind_pocsi)%field(i,j,k,Time%taum1))

         enddo  ! i
       enddo  ! k


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
       par_dia(:,:) = 0.0
       par_diz(:,:) = 0.0
       par_phymld(:,:) = 0.0
       par_diamld(:,:) = 0.0
       par_dizmld(:,:) = 0.0
       par_tot_mldsum = 0.0
       par_phy_mldsum = 0.0
       par_dia_mldsum = 0.0
       par_diz_mldsum = 0.0
       par_z_mldsum = 0.0
       phy_lday(:,:) = 1e-2
       dia_lday(:,:) = 1e-2
       diz_lday(:,:) = 1e-2
       chl_lday(:,:) = 1e-2
       keuphot(:) = 4  ! set k level of euphotic zone to 10 metres
       daylight(:) = 0.0

       do i=isc,iec
         ! Calculate the daylength (hours per day) given the latitude and day of year
         day_angle = 2.0*pi*(days_in_this_year-1)/365.0
         declination = 0.409 * sin(2*pi*days_in_this_year/365 - 1.39)
         cos_hour_angle = max(min(-tan(declination) * tan(grid%yt(i,j) * pi/180.0), 1.0), -1.0)
         daylight(i) = 2 * (acos(cos_hour_angle) * 180 / pi) / 15.0
       enddo

       do k=1,grid%nk
         do i=isc,iec

           ! Attenuation coefficients for Red, Green and Blue light given chlorophyll concentration
           zchl = MAX(0.05, MIN(10.0, biotr(i,k,id_pchl)+biotr(i,k,id_dchl) ))
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
             par_dia(i,k) = 1.62 * par_rgb(i,k,1) + 0.74 * par_rgb(i,k,2) + 0.63 * par_rgb(i,k,3)
             par_diz(i,k) = 1.42 * par_rgb(i,k,1) + 0.84 * par_rgb(i,k,2) + 0.83 * par_rgb(i,k,3) !!! TMP
           else
             par_tot(i,k) = 0.0
             par_phy(i,k) = 0.0
             par_dia(i,k) = 0.0
             par_diz(i,k) = 0.0
           endif
           par_eup(i,k) = par_rgb(i,k,1) + par_rgb(i,k,2) + par_rgb(i,k,3)

           ! Integrated light level in the mixed layer
           if (grid%zw(k) .le. mld(i,j)) then
             par_tot_mldsum(i) = par_tot_mldsum(i) + par_tot(i,k) * thickness%dzt(i,j,k)
             par_phy_mldsum(i) = par_phy_mldsum(i) + par_phy(i,k) * thickness%dzt(i,j,k)
             par_dia_mldsum(i) = par_dia_mldsum(i) + par_dia(i,k) * thickness%dzt(i,j,k)
             par_diz_mldsum(i) = par_diz_mldsum(i) + par_diz(i,k) * thickness%dzt(i,j,k)
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

       ! Euphotic zone depth, making sure it's not greater than the depth of the grid
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
       mumax_dia(:,:) = 0.0     ! microphytoplankton maximum growth rate (temperature-dependent)
       mumax_diz(:,:) = 0.0     ! diazotroph maximum growth rate (temperature-dependent)

       do k=1,grid%nk
         do i=isc,iec

! pjb - Calculate average light level in the mixed layer
           if (grid%zw(k) .le. mld(i,j)) then
             zval = 1.0/par_z_mldsum(i)
             par_phymld(i,k) = par_phy_mldsum(i) * zval
             par_diamld(i,k) = par_dia_mldsum(i) * zval
             par_dizmld(i,k) = par_diz_mldsum(i) * zval
           else
             par_phymld(i,k) = par_phy(i,k)
             par_diamld(i,k) = par_dia(i,k)
             par_dizmld(i,k) = par_diz(i,k)
           endif

! pjb - Calculate impact of daylength on phytoplankton and chlorophyll production
           zval = MAX(1.0, daylight(i))
           if (grid%zw(k) .le. mld(i,j)) then
             zval = zval * MIN(1.0, (zeuphot(i,j)+epsi)/(mld(i,j)+epsi))
           endif
           chl_lday(i,k) = zval / 24.0
           phy_lday(i,k) = 1.5 * zval / (12.0 + zval)
           dia_lday(i,k) = 1.5 * zval / (12.0 + zval)
           diz_lday(i,k) = 1.5 * zval / (12.0 + zval)

           ! Save total PAR to radbio array for diagnostic output
           radbio3d(i,j,k) = par_tot(i,k)

           ! Temperature-dependent maximum growth rate (Eppley curve)
           ztem = t_prog(index_temp)%field(i,j,k,time%tau)
           mumax_phy(i,k) = abioa_phy(i,j) * bbioa(i,j)**(cbioa(i,j) * ztem) / 86400.0
           mumax_dia(i,k) = abioa_dia(i,j) * bbioa(i,j)**(cbioa(i,j) * ztem) / 86400.0
           mumax_diz(i,k) = abioa_diz(i,j) * exp(- (ztem - beta_diz(i,j)) /                        &
                                                   (8.0 - 0.5*(ztem - beta_diz(i,j)))) / 86400.0

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
         dia_chlc(:,:) = 0.0      ! chlorophyll to carbon ratio of microphytoplankton
         diz_chlc(:,:) = 0.0      ! chlorophyll to carbon ratio of diazotroph
         phy_FeC(:,:) = 0.0       ! iron to carbon ratio of phytoplankton
         dia_FeC(:,:) = 0.0       ! iron to carbon ratio of microphytoplankton
         diz_FeC(:,:) = 0.0       ! iron to carbon ratio of diazotroph
         dia_SiC(:,:) = 0.0       ! Silicon to carbon ratio of microphytoplankton
         dia_SiCupt(:,:) = 0.0    ! Silicon to carbon ratio of microphytoplankton uptake
         zoo_FeC(:,:) = 0.0       ! iron to carbon ratio of zooplankton
         mes_FeC(:,:) = 0.0       ! iron to carbon ratio of mesozooplankton
         det_FeC(:,:) = 0.0       ! iron to carbon ratio of slow sinking detritus
         poc_FeC(:,:) = 0.0       ! iron to carbon ratio of fast sinking detritus
         poc_SiC(:,:) = 0.0       ! Silicon to carbon ratio of fast sinking detritus
         phy_Feupt(:,:) = 0.0     ! iron uptake rate of phytoplankton
         dia_Feupt(:,:) = 0.0     ! iron uptake rate of microphytoplankton
         diz_Feupt(:,:) = 0.0     ! iron uptake rate of diazotroph
         phy_lpar(:,:) = 1.0      ! light limitation of phytoplankton growth
         dia_lpar(:,:) = 1.0      ! light limitation of microphytoplankton growth
         diz_lpar(:,:) = 1.0      ! light limitation of diazotroph growth
         mu_phy(:,:) = 0.0        ! phytoplankton realised growth rate (light limited)
         mu_dia(:,:) = 0.0        ! microphytoplankton realised growth rate (light limited)
         mu_diz(:,:) = 0.0        ! diazotroph realised growth rate (light limited)
         phy_chl_lpar(:,:) = 0.0  ! light limitation of chlorophyll production by phytoplankton
         dia_chl_lpar(:,:) = 0.0  ! light limitation of chlorodiall production by microphytoplankton
         diz_chl_lpar(:,:) = 0.0  ! light limitation of chlorodiall production by diazotroph
         mu_pchl(:,:) = 0.0       ! Growth of chlorophyll in phytoplankton
         mu_dchl(:,:) = 0.0       ! Growth of chlorophyll in microphytoplanktons
         mu_zchl(:,:) = 0.0       ! Growth of chlorophyll in diazotrophs
         phy_limnit(:,:) = 1.0    ! Phytoplankton nitrogen limitation
         dia_limnit(:,:) = 1.0    ! Microphytoplankton nitrogen limitation
         diz_limnit(:,:) = 1.0    ! diazotroph nitrogen limitation
         phy_limfe(:,:) = 1.0     ! Phytoplankton iron limitation
         dia_limfe(:,:) = 1.0     ! Microphytoplankton iron limitation
         diz_limfe(:,:) = 1.0     ! diazotroph iron limitation
         phy_limpho(:,:) = 1.0    ! Phytoplankton phosphorus limitation
         dia_limpho(:,:) = 1.0    ! Microphytoplankton phosphorus limitation
         diz_limpho(:,:) = 1.0    ! diazotroph phosphorus limitation
         dia_limsil(:,:) = 1.0    ! Microphytoplankton silicic acid limitation
         si_xlabile(:,:) = 0.5    ! Initialise that half of all particulate silicate is labile
         nitrif(:,:) = 0.0        ! nitrification rate (NH4 --> NO3)
         denitr(:,:) = 0.0        ! denitrification rate (NO3 --> N2)

         ! Calculate the ratio of labile to non-labile particulate silicate
         si_lab2non(:,:) = si_xlabile(:,:)/(1.0-si_xlabile(:,:))   
        
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
             if (id_nh4.ne.0)   bionh4 = max(0.0,biotr(i,k,id_nh4))
             if (id_po4.ne.0)   biopo4 = max(0.0,biotr(i,k,id_po4))
             if (id_sil.ne.0)   biosil = max(0.0,biotr(i,k,id_sil))
             if (id_poc.ne.0)   biopoc = max(0.0,biotr(i,k,id_poc))
             if (id_dia.ne.0)   biodia = max(0.0,biotr(i,k,id_dia))
             if (id_diz.ne.0)   biodiz = max(0.0,biotr(i,k,id_diz))
             if (id_mes.ne.0)   biomes = max(0.0,biotr(i,k,id_mes))
             if (id_caco3.ne.0) biocaco3 = max(0.0,biotr(i,k,id_caco3))
             if (id_fe.ne.0)    biofe = max(0.0,biotr(i,k,id_fe))
             if (id_pchl.ne.0)  biopchl = max(0.0,biotr(i,k,id_pchl))
             if (id_dchl.ne.0)  biodchl = max(0.0,biotr(i,k,id_dchl))
             if (id_zchl.ne.0)  biozchl = max(0.0,biotr(i,k,id_zchl))
             if (id_phyfe.ne.0) biophyfe = max(0.0,biotr(i,k,id_phyfe))
             if (id_diafe.ne.0) biodiafe = max(0.0,biotr(i,k,id_diafe))
             if (id_dizfe.ne.0) biodizfe = max(0.0,biotr(i,k,id_dizfe))
             if (id_diasi.ne.0) biodiasi = max(0.0,biotr(i,k,id_diasi))
             if (id_zoofe.ne.0) biozoofe = max(0.0,biotr(i,k,id_zoofe))
             if (id_mesfe.ne.0) biomesfe = max(0.0,biotr(i,k,id_mesfe))
             if (id_detfe.ne.0) biodetfe = max(0.0,biotr(i,k,id_detfe))
             if (id_pocfe.ne.0) biopocfe = max(0.0,biotr(i,k,id_pocfe))
             if (id_pocsi.ne.0) biopocsi = max(0.0,biotr(i,k,id_pocsi))

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.2] - Get important ratios and values                                                !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             ! Get important ratios
             if (biophy*biopchl.gt.0.0) phy_chlc(i,k) = biopchl / (biophy * 12)
             if (biodia*biodchl.gt.0.0) dia_chlc(i,k) = biodchl / (biodia * 12)
             if (biodiz*biozchl.gt.0.0) diz_chlc(i,k) = biozchl / (biodiz * 12)
             if (id_fe.ne.0) then 
               if (biophy*biophyfe.gt.0.0) phy_FeC(i,k) = biophyfe/biophy
               if (biodia*biodiafe.gt.0.0) dia_FeC(i,k) = biodiafe/biodia
               if (biodiz*biodizfe.gt.0.0) diz_FeC(i,k) = biodizfe/biodiz
               if (biozoo*biozoofe.gt.0.0) zoo_FeC(i,k) = biozoofe/biozoo
               if (biomes*biomesfe.gt.0.0) mes_FeC(i,k) = biomesfe/biomes
               if (biodet*biodetfe.gt.0.0) det_FeC(i,k) = biodetfe/biodet
               if (biopoc*biopocfe.gt.0.0) poc_FeC(i,k) = biopocfe/biopoc
             endif
             if (id_sil.ne.0) then
               if (biodia*biodiasi.gt.0.0) dia_SiC(i,k) = biodiasi/biodia
               if (biopoc*biopocsi.gt.0.0) poc_SiC(i,k) = biopocsi/biopoc
             endif

             ! Temperature and salinity
             ztem = t_prog(index_temp)%field(i,j,k,time%tau)
             ztemk = ztem + 273.15
             zsal = t_prog(index_salt)%field(i,j,k,time%tau)


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.3] - Dissolved Fe chemistry                                                         !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
          !!! THINGS NEEDING TO BE DONE:
          !!!   - Variable ligand concentration rather than constant?
          !!!   - Variable pH value (currently set to constant 7.9)

             ! 1. Determine equilibrium fractionation of total dFe into Fe` and L-Fe
             ligand = 0.7  ! nM ligand concentration
             fe_keq = 10**( 17.27 - 1565.7 / ztemk ) * 1e-9
             fe_III = ( -( 1. + fe_keq * ligand - fe_keq * biofe )                                 &
                        + SQRT( ( 1. + fe_keq * ligand - fe_keq * biofe )**2                       &
                                + 4. * biofe * fe_keq) ) / ( 2. * fe_keq )
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
             precip = max(0.0, ( fe_III - fe3sol ) ) * 0.01/86400.0

             ! 3. Scavenging of Fe` (Fe` --> DET & POC)
             partic = (biodet + biopoc + biocaco3 + biopocsi)
             scaven = fe_III * (3e-5 + 0.005 * partic) / 86400.0 * 1e-3
             scadet = scaven * (biodet+epsi) / (partic+epsi)
             scapoc = scaven * (biopoc+epsi) / (partic+epsi)

             ! 4. Increased loss of Fe in near-coastal environments
             zval = min(1.0, 1000.0 / (grid%zt(grid%kmt(i,j))+1) )
             feloss = 1e-4 / 86400.0 * zval * biofe

             ! 5. Coagulation of colloidal Fe (nM) into small and large particles (uM)
             fe_col = fe_lig * 0.5
             if (grid%zw(k) .le. mld(i,j)) then
               zval =        ( 0.369 * 0.3 * 40e-6 + 102.4 * biodet*1e-6 ) + ( 114. * 0.3 * 40e-6 )
             else
               zval = 0.01 * ( 0.369 * 0.3 * 40e-6 + 102.4 * biodet*1e-6 ) + ( 114. * 0.3 * 40e-6 )
             endif
             fe2det = fe_col * zval / 86400.0
             if (grid%zw(k) .le. mld(i,j)) then
               zval =        ( 3.53e3 * biopoc*1e-6 )
             else
               zval = 0.01 * ( 3.53e3 * biopoc*1e-6 )
             endif
             fe2poc = fe_col * zval / 86400

          
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.4] - Nutrient limitation terms for phytoplankton growth                             !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             ! Variable half-saturation Ks (change with population biomass)
             phy_conc1 = max(0.0, biophy - phybiot(i,j))  ! phybiot is threshold of phyto biomass
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
             dia_conc1 = max(0.0, biodia - diabiot(i,j))  ! diabiot is threshold of diatom biomass 
             dia_conc2 = biodia - dia_conc1
             dia_k_nit = max( diakn(i,j),                                                          &
                              ( dia_conc2 * diakn(i,j) + dia_conc1 * 3.0 * diakn(i,j) )            &
                              / ( biodia + epsi ) )
             dia_k_fer = max( diakf(i,j),                                                          &
                              ( dia_conc2 * diakf(i,j) + dia_conc1 * 3.0 * diakf(i,j) )            &
                              / ( biodia + epsi ) )
             dia_k_pho = max( diakp(i,j),                                                          &
                              ( dia_conc2 * diakp(i,j) + dia_conc1 * 3.0 * diakp(i,j) )            &
                              / ( biodia + epsi ) )
             dia_k_sil = max(8.0, 1.0 + 7.0*biosil**2 / (biosil**2 + 16.5**2) )
             diz_conc1 = max(0.0, biodiz - dizbiot(i,j))  ! dizbiot is threshold of diztom biomass 
             diz_conc2 = biodiz - diz_conc1
             diz_k_nit = max( dizkn(i,j),                                                          &
                              ( diz_conc2 * dizkn(i,j) + diz_conc1 * 3.0 * dizkn(i,j) )            &
                              / ( biodiz + epsi ) )
             diz_k_fer = max( dizkf(i,j),                                                          &
                              ( diz_conc2 * dizkf(i,j) + diz_conc1 * 3.0 * dizkf(i,j) )            &
                              / ( biodiz + epsi ) )
             diz_k_pho = max( dizkp(i,j),                                                          &
                              ( diz_conc2 * dizkp(i,j) + diz_conc1 * 3.0 * dizkp(i,j) )            &
                              / ( biodiz + epsi ) )


             ! Nitrogen (NH4 & NO3)
             phy_zlimnh4 = bionh4 / (phy_k_nit + bionh4)
             phy_zlimno3 = biono3 / (phy_k_nit + biono3)
             phy_zlimnit = ( bionh4 + biono3 ) / ( phy_k_nit + bionh4 + biono3 )
             phy_limnh4 = phy_zlimnit * 5.0 * phy_zlimnh4 /                                        &
                          ( phy_zlimno3 + 5.0 * phy_zlimnh4 + epsi )
             phy_limno3 = phy_zlimnit * phy_zlimno3       /                                        &
                          ( phy_zlimno3 + 5.0 * phy_zlimnh4 + epsi )
             phy_limnit(i,k) = phy_limnh4 + phy_limno3

             dia_zlimnh4 = bionh4 / (dia_k_nit + bionh4)
             dia_zlimno3 = biono3 / (dia_k_nit + biono3)
             dia_zlimnit = ( bionh4 + biono3 ) / ( dia_k_nit + bionh4 + biono3 )
             dia_limnh4 = dia_zlimnit * 5.0 * dia_zlimnh4 /                                        & 
                          ( dia_zlimno3 + 5.0 * dia_zlimnh4 + epsi )
             dia_limno3 = dia_zlimnit * dia_zlimno3       /                                        &
                          ( dia_zlimno3 + 5.0 * dia_zlimnh4 + epsi )
             dia_limnit(i,k) = dia_limnh4 + dia_limno3

             diz_zlimnh4 = bionh4 / (diz_k_nit + bionh4)
             diz_zlimno3 = biono3 / (diz_k_nit + biono3)
             diz_zlimnit = ( bionh4 + biono3 ) / ( diz_k_nit + bionh4 + biono3 )
             diz_limnh4 = diz_zlimnit * 5.0 * diz_zlimnh4 /                                        & 
                          ( diz_zlimno3 + 5.0 * diz_zlimnh4 + epsi )
             diz_limno3 = diz_zlimnit * diz_zlimno3       /                                        &
                          ( diz_zlimno3 + 5.0 * diz_zlimnh4 + epsi )
             diz_limnit(i,k) = diz_limnh4 + diz_limno3


             ! Iron
             if (id_fe.ne.0) then 
               if (biopchl*biophy.gt.0.0) then
                 phy_minFeC = 0.0016 / 55.85 * biopchl/biophy +                                    &
                              1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * phy_limnit(i,k) +       &
                              1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * phy_limno3
               else
                 phy_minFeC = 0.0016 / 55.85 * phyminqc(i,j)*12 +                                  &
                              1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * phy_limnit(i,k) +       &
                              1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * phy_limno3
               endif
               if (biodchl*biodia.gt.0.0) then
                 dia_minFeC = 0.0016 / 55.85 * biodchl/biodia +                                    &
                              1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * dia_limnit(i,k) +       &
                              1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * dia_limno3
               else
                 dia_minFeC = 0.0016 / 55.85 * diaminqc(i,j)*12 +                                  &
                              1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * dia_limnit(i,k) +       &
                              1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * dia_limno3
               endif
               if (biozchl*biodiz.gt.0.0) then
                 diz_minFeC = 0.0016 / 55.85 * biozchl/biodiz +                                    &
                              1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * diz_limnit(i,k) +       &
                              1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * diz_limno3
               else
                 diz_minFeC = 0.0016 / 55.85 * dizminqc(i,j)*12 +                                  &
                              1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * diz_limnit(i,k) +       &
                              1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * diz_limno3
               endif
               phy_limfe(i,k) = min(1.0, max(0.0, (phy_FeC(i,k) - phy_minFeC)/phyoptqf(i,j) ))
               dia_limfe(i,k) = min(1.0, max(0.0, (dia_FeC(i,k) - dia_minFeC)/diaoptqf(i,j) ))
               diz_limfe(i,k) = min(1.0, max(0.0, (diz_FeC(i,k) - diz_minFeC)/dizoptqf(i,j) ))
             else
               phy_limfe(i,k) = 1.0
               dia_limfe(i,k) = 1.0
               diz_limfe(i,k) = 1.0
             endif


             ! Phosphorus (PO4)
             if (id_po4.gt.0) then
               phy_limpho(i,k) = biopo4 / ( biopo4 + phy_k_pho )
               dia_limpho(i,k) = biopo4 / ( biopo4 + dia_k_pho )
               diz_limpho(i,k) = biopo4 / ( biopo4 + diz_k_pho )
             endif


             ! Silicic acid (Sil)
             dia_limsil(i,k) = biosil / ( biosil + dia_k_sil )


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
             phy_pisl2= phy_pisl(i,k) / ( (1.0 + phylmor(i,j)*86400.0) * max(daylight(i),1.0)/24.0 )
             phy_lpar(i,k) = (1.0 - EXP(-phy_pisl2 * par_phy(i,k))) * phy_lday(i,k)
             mu_phy(i,k) = mumax_phy(i,k) * phy_lpar(i,k)

             dia_pisl(i,k) = MAX(alphabio(i,j) * dia_chlc(i,k), alphabio(i,j) * diaminqc(i,j))
             dia_pisl2= dia_pisl(i,k) / ( (1.0 + phylmor(i,j)*86400.0) * max(daylight(i),1.0)/24.0 )
             dia_lpar(i,k) = (1.0 - EXP(-dia_pisl2 * par_dia(i,k))) * dia_lday(i,k)
             mu_dia(i,k) = mumax_dia(i,k) * dia_lpar(i,k)

             diz_pisl(i,k) = MAX(alphabio(i,j) * diz_chlc(i,k), alphabio(i,j) * dizminqc(i,j))
             diz_pisl2= diz_pisl(i,k) / ( (1.0 + phylmor(i,j)*86400.0) * max(daylight(i),1.0)/24.0 )
             diz_lpar(i,k) = (1.0 - EXP(-diz_pisl2 * par_diz(i,k))) * diz_lday(i,k)
             mu_diz(i,k) = mumax_diz(i,k) * diz_lpar(i,k)


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
             dia_limnut = min( dia_limnit(i,k), dia_limpho(i,k), dia_limfe(i,k) )
             diz_limnut = min( 1.0, diz_limpho(i,k), diz_limfe(i,k) )

             u_phy = phy_limnut * mu_phy(i,k)
             u_dia = min(dia_limnut, dia_limsil(i,k)) * mu_dia(i,k)
             u_diz = diz_limnut * mu_diz(i,k)


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.8] - Si:C ratio for diatom productivity                                             !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             zval = biosil / (biosil + 2.0)
             sicfac1 = 1.0 + 4.4 * exp(-4.23 * min(dia_lday(i,k), dia_limnut)) *                   &
                      max(0.0, min(1.0, 2.2 * (zval-0.5) ) )
             sicfac2 = 1.0 + 2.0 * biosil**3 / (biosil**3 + 20.0**3) 
             ! sicfac1 --> [1.0,5.4] as SiOH4 increases and limitation strengthens
             ! sicfac2 --> [1.0,3.0] as SiOH4 increases 
             dia_SiCupt(i,k) = 0.14 * zval * sicfac1 * sicfac2


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.9] - Silicate dissolution rate                                                      !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             sio3eq = exp( log(10.0) * (6.44 - 968. / ztemk) )
             sisat1 = max( epsi, (sio3eq - biosil) / (sio3eq + epsi) )
             sisat2 = (1. + ztem / 400.0)**37.0
             sidisr = 0.225 * (1.0 + ztem/15.0) * sisat1 + 0.775*sisat2*sisat1**9.25 
             if (grid%zt(k) .gt. max(zeuphot(i,j), mld(i,j))) then
               si_lab2non(i,k) = si_lab2non(i,k-1) * exp( -0.5 * (0.03 - 0.003) * sidisr *         &
                                 thickness%dzt(i,j,k) / (wdetbio(i,j)*10*86400.0) )
               si_xlabile(i,k) = si_lab2non(i,k) / (1.0 + si_lab2non(i,k))
               si_lab2non(i,k) = si_lab2non(i,k) * exp( -0.5 * (0.03 - 0.003) * sidisr *           &
                                 thickness%dzt(i,j,k) / (wdetbio(i,j)*10*86400.0) )
             endif 
             siremin = 0.03 * si_xlabile(i,k) + 0.003 * (1.-si_xlabile(i,k)) * sidisr


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.10] - Growth in chlorophyll (mg Chl / m3)                                            !
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

             dia_chl_pisl = dia_pisl(i,k) / ( mumax_dia(i,k) * 86400.0                             &
                                              * max(daylight(i),1.0)/24.0                          &
                                              * (1.0 - dia_limnut) + epsi)
             dia_chl_lpar(i,k) = (1.0 - EXP( -dia_chl_pisl * par_diamld(i,k))) * dia_lday(i,k)
             mumin_dchl = diaminqc(i,j) * u_dia * biodia * 12.
             muopt_dchl = diaoptqc(i,j) * u_dia * biodia * 12.
             mu_dchl(i,k) = (muopt_dchl - mumin_dchl) * dia_chl_lpar(i,k) * dia_limnut
             if (dia_pisl(i,k)*par_diamld(i,k).gt.0.0) then
               mu_dchl(i,k) = mumin_dchl + mu_dchl(i,k) / (dia_pisl(i,k) * par_diamld(i,k))
             endif

             diz_chl_pisl = diz_pisl(i,k) / ( mumax_diz(i,k) * 86400.0                             &
                                              * max(daylight(i),1.0)/24.0                          &
                                              * (1.0 - diz_limnut) + epsi)
             diz_chl_lpar(i,k) = (1.0 - EXP( -diz_chl_pisl * par_dizmld(i,k))) * diz_lday(i,k)
             mumin_zchl = dizminqc(i,j) * u_diz * biodiz * 12.
             muopt_zchl = dizoptqc(i,j) * u_diz * biodiz * 12.
             mu_zchl(i,k) = (muopt_zchl - mumin_zchl) * diz_chl_lpar(i,k) * diz_limnut
             if (diz_pisl(i,k)*par_dizmld(i,k).gt.0.0) then
               mu_zchl(i,k) = mumin_zchl + mu_zchl(i,k) / (diz_pisl(i,k) * par_dizmld(i,k))
             endif


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.11] - Phytoplankton uptake of Fe (nmol Fe / m3)                                     !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
         
             ! 1. Maximum iron content of phytoplankton biomass
             ! 2. Ensure that Fe uptake increases or decreases in response to standing quotas
             ! 3. Iron uptake by phytoplankton

             if (id_fe.ne.0) then 
               if (phy_FeC(i,k).gt.0.0) then 
                 phy_maxQFe = biophy * phymaxqf(i,j)
                 phy_Feupt_upreg   = (4.0 - 4.5 * phy_limfe(i,k) / (phy_limfe(i,k) + 0.5)) 
                 phy_Feupt_downreg = max(0.0, (1.0 - biophyfe/phy_maxQFe) /                        &
                                              abs(1.05 - biophyfe/phy_maxQFe) )
                 phy_Feupt(i,k) = biophy * mumax_phy(i,k) * phymaxqf(i,j) *                        &
                                  biofe / (biofe + phy_k_fer) * phy_Feupt_downreg * phy_Feupt_upreg
               endif
               if (dia_FeC(i,k).gt.0.0) then 
                 dia_maxQFe = biodia * diamaxqf(i,j)
                 dia_Feupt_upreg   = (4.0 - 4.5 * dia_limfe(i,k) / (dia_limfe(i,k) + 0.5)) 
                 dia_Feupt_downreg = max(0.0, (1.0 - biodiafe/dia_maxQFe) /                        &
                                              abs(1.05 - biodiafe/dia_maxQFe) )
                 dia_Feupt(i,k) = biodia * mumax_dia(i,k) * diamaxqf(i,j) *                        &
                                  biofe / (biofe + dia_k_fer) * dia_Feupt_downreg * dia_Feupt_upreg
               endif
               if (diz_FeC(i,k).gt.0.0) then 
                 diz_maxQFe = biodiz * dizmaxqf(i,j)
                 diz_Feupt_upreg   = (4.0 - 4.5 * diz_limfe(i,k) / (diz_limfe(i,k) + 0.5)) 
                 diz_Feupt_downreg = max(0.0, (1.0 - biodizfe/diz_maxQFe) /                        &
                                              abs(1.05 - biodizfe/diz_maxQFe) )
                 diz_Feupt(i,k) = biodiz * mumax_diz(i,k) * dizmaxqf(i,j) *                        &
                                  biofe / (biofe + diz_k_fer) * diz_Feupt_downreg * diz_Feupt_upreg
               endif
             endif
          

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.12] - Temperature dependence of heterotrophy / respiration                          !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             fbc = abioh(i,j) * bbioh(i,j)**(cbioh(i,j) * ztem)


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.13] - Grazing functions                                                             !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             biophy_min = MAX(0.0, biophy - 0.001)  ! Minimum Phy beneath which no grazing 
             biodia_min = MAX(0.0, biodia - 0.001)  ! Minimum dia beneath which no grazing 
             biodiz_min = MAX(0.0, biodiz - 0.001)  ! Minimum diz beneath which no grazing 
             biodet_min = MAX(0.0, biodet - 0.001)  ! Minimum Det beneath which no grazing 
             biopoc_min = MAX(0.0, biopoc - 0.001)  ! Minimum poc beneath which no grazing 
             biozoo_min = MAX(0.0, biozoo - 0.001)  ! Minimum mes beneath which no grazing 

             ! preferences and prey availability
             zooprey = zprefphy(i,j)*biophy_min + zprefdia(i,j)*biodia_min                         &
                       + zprefdiz(i,j)*biodiz_min +                                                &
                       + zprefdet(i,j)*biodet_min + zprefpoc(i,j)*biopoc_min
             mesprey = mprefphy(i,j)*biophy_min + mprefdia(i,j)*biodia_min                         &
                       + mprefdet(i,j)*biodet_min + mprefpoc(i,j)*biopoc_min                       &
                       + mprefzoo(i,j)*biozoo_min + mprefdiz(i,j)*biodiz_min

             g_zoo = zoogmax(i,j)*fbc * epszoo(i,j) * zooprey*zooprey /                            &
                     ( zoogmax(i,j)*fbc + epszoo(i,j) * zooprey*zooprey )
             g_mes = mesgmax(i,j)*fbc * epsmes(i,j) * mesprey*mesprey /                            &
                     ( mesgmax(i,j)*fbc + epsmes(i,j) * mesprey*mesprey )


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.14] - Mortality scalers                                                             !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             ! Depress mortality in nutrient limited regime (i.e., investment in cell maintenance)
             phy_slmor = MIN(1.0, MAX(0.0,                                                         &
                         min(phy_limnit(i,k), phy_limpho(i,k), phy_limfe(i,k))/0.3))
             dia_slmor = MIN(1.0, MAX(0.0,                                                         &
                         min(dia_limnit(i,k), dia_limpho(i,k), dia_limfe(i,k))/0.3))
             diz_slmor = MIN(1.0, MAX(0.0,                                                         &
                         min(diz_limnit(i,k), diz_limpho(i,k), diz_limfe(i,k))/0.3))


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.15] - Bacterial processes                                                           !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             ! Denitrification and nitrification (set at maximum of 5% of NH4 per day)
             nitrif(i,k) = 0.05/86400.0 * bionh4
             denitr(i,k) = 0.05/86400.0 * biono3 * min(1.0, exp(-biooxy))


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.16] - Collect major source and sink terms for tendency equations                    !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             ! PRIMARY PRODUCTION
             f11phy = u_phy * biophy                                        * biomask(i,k)         ! DIC --> PHY
             f11dia = u_dia * biodia                                        * biomask(i,k)         ! DIC --> DIA
             f11diz = u_diz * biodiz                                        * biomask(i,k)         ! DIC --> DIZ

             ! GRAZING
             if (zooprey.gt.0.0) then
               z21phy = g_zoo * biozoo * (zprefphy(i,j)*biophy_min/zooprey) * biomask(i,k)         ! PHY --> ZOO & DET & POC
               z21dia = g_zoo * biozoo * (zprefdia(i,j)*biodia_min/zooprey) * biomask(i,k)         ! DIA --> ZOO & DET & POC
               z21diz = g_zoo * biozoo * (zprefdiz(i,j)*biodiz_min/zooprey) * biomask(i,k)         ! DIZ --> ZOO & DET & POC
               z21det = g_zoo * biozoo * (zprefdet(i,j)*biodet_min/zooprey) * biomask(i,k)         ! DET --> ZOO & DET
               z21poc = g_zoo * biozoo * (zprefpoc(i,j)*biopoc_min/zooprey) * biomask(i,k)         ! POC --> ZOO & POC
             else
               z21phy = 0.0 
               z21dia = 0.0 
               z21diz = 0.0 
               z21det = 0.0
               z21poc = 0.0
             endif
             if (mesprey.gt.0.0) then
               m21phy = g_mes * biomes * (mprefphy(i,j)*biophy_min/mesprey) * biomask(i,k)         ! PHY --> MES & DET & POC
               m21dia = g_mes * biomes * (mprefdia(i,j)*biodia_min/mesprey) * biomask(i,k)         ! DIA --> MES & DET & POC
               m21diz = g_mes * biomes * (mprefdiz(i,j)*biodiz_min/mesprey) * biomask(i,k)         ! DIZ --> MES & DET & POC
               m21det = g_mes * biomes * (mprefdet(i,j)*biodet_min/mesprey) * biomask(i,k)         ! DET --> MES & DET
               m21poc = g_mes * biomes * (mprefpoc(i,j)*biopoc_min/mesprey) * biomask(i,k)         ! POC --> MES & POC
               m21zoo = g_mes * biomes * (mprefzoo(i,j)*biozoo_min/mesprey) * biomask(i,k)         ! ZOO --> MES & POC
             else
               m21phy = 0.0 
               m21dia = 0.0 
               m21diz = 0.0 
               m21det = 0.0
               m21poc = 0.0
               m21zoo = 0.0
             endif

             ! MORTALITY
             f22phy = phylmor(i,j) * fbc * phy_slmor * biophy               * biomask(i,k)         ! PHY --> DIC
             f22dia = phylmor(i,j) * fbc * dia_slmor * biodia               * biomask(i,k)         ! DIA --> DIC
             f22diz = phylmor(i,j) * fbc * diz_slmor * biodiz               * biomask(i,k)         ! DIZ --> DIC
             f23phy = phyqmor(i,j) * biophy * biophy                        * biomask(i,k)         ! PHY --> DET & POC
             f23dia = phyqmor(i,j) * biodia * biodia                        * biomask(i,k)         ! DIA --> DET & POC
             f23diz = phyqmor(i,j) * biodiz * biodiz                        * biomask(i,k)         ! DIZ --> DET & POC
             z31 = zooexcr(i,j) * biozoo * fbc                              * biomask(i,k)         ! ZOO --> DIC
             z32 = zooqmor(i,j) * biozoo*biozoo                             * biomask(i,k)         ! ZOO --> DET & POC
             m31 = zooexcr(i,j) * biomes * fbc                              * biomask(i,k)         ! MES --> DIC
             m32 = mesqmor(i,j) * biomes*biomes                             * biomask(i,k)         ! MES --> DET & POC

             ! REMINERALISATION
             f41det = detlrem(i,j) * biodet * fbc                           * biomask(i,k)         ! DET --> DIC
             f41poc = detlrem(i,j) * biopoc * fbc                           * biomask(i,k)         ! POC --> DIC
             f41pocsi = siremin/86400.0 * biopocsi                          * biomask(i,k)         ! POCSi --> Sil
             if (grid%zw(k) .ge. 180) then
               f41det = f41det*.5   ! reduce decay below 180m
               f41poc = f41poc*.5   ! reduce decay below 180m
             endif
             if (id_caco3.ne.0) f51 = caco3rem(i,j)* biocaco3 * biomask(i,k)


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.17] - Apply sources and sinks to tendency equations                                 !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             !!!!!!!!!!!!!!!!!!!
             ! oxygen equation !
             !!!!!!!!!!!!!!!!!!!
             biotr(i,k,id_o2)   = biotr(i,k,id_o2) - biomask(i,k) * 172./122. * dtsb * (           &
                                                     f41det + f41poc + z31 + m31                   &
                                                     + f22phy + f22dia + f22diz                    &
                                                     - f11phy - f11dia - f11diz )

             !!!!!!!!!!!!!!!!!!!!!!
             ! phosphate equation !
             !!!!!!!!!!!!!!!!!!!!!!
             if (id_po4.ne.0)                                                                      &
             biotr(i,k,id_po4)  = biotr(i,k,id_po4) + dtsb * 1./122. * (                           &
                                                      f41det + f41poc                              &
                                                    + z31 + m31 + f22phy + f22dia + f22diz         &
                                                    - f11phy - f11dia - f11diz)

             !!!!!!!!!!!!!!!!!!!!!!!!!
             ! silicic acid equation !
             !!!!!!!!!!!!!!!!!!!!!!!!!
             if (id_sil.ne.0)                                                                      &
             biotr(i,k,id_sil)  = biotr(i,k,id_sil) + dtsb * ( f41pocsi + f22dia * dia_SiC(i,k)    &
                                                               - f11dia * dia_SiCupt(i,k) ) 

             !!!!!!!!!!!!!!!!!!!!!
             ! ammonium equation !
             !!!!!!!!!!!!!!!!!!!!!
             if (id_nh4.ne.0)                                                                      &
             biotr(i,k,id_nh4)  = biotr(i,k,id_nh4) - dtsb * nitrif(i,k) + dtsb * 16./122. * (     &
                                                      f41det + f41poc                              &
                                                      + z31 + m31 + f22phy + f22dia + f22diz       &
                                                    - f11phy *(phy_limnh4)/(phy_limnit(i,k)+epsi)  &
                                                    - f11dia *(dia_limnh4)/(dia_limnit(i,k)+epsi)  &
                                                    - f11diz * diz_limnh4 )


             !!!!!!!!!!!!!!!!!!!!
             ! nitrate equation !
             !!!!!!!!!!!!!!!!!!!!
             biotr(i,k,id_no3)  = biotr(i,k,id_no3) + dtsb * ( nitrif(i,k) - denitr(i,k) )         &
                                                    - dtsb * 16./122. * (                          &
                                                      f11phy *(phy_limno3)/(phy_limnit(i,k)+epsi)  &
                                                    + f11dia *(dia_limno3)/(dia_limnit(i,k)+epsi)  &
                                                    + f11diz * diz_limno3 )  

             !!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! dissolved iron equation !
             !!!!!!!!!!!!!!!!!!!!!!!!!!!
             if (id_fe.ne.0)                                                                       &
             biotr(i,k,id_fe)   = biotr(i,k,id_fe)  + dtsb * (f41det*det_FeC(i,k) +                &
                                                              f41poc*poc_FeC(i,k) +                &
                                                              z31*zoo_FeC(i,k) +                   &
                                                              m31*mes_FeC(i,k) +                   &
                                                              f22phy*phy_FeC(i,k) +                &
                                                              f22dia*dia_FeC(i,k) -                &
                                                              f22diz*diz_FeC(i,k) -                &
                                                              phy_Feupt(i,k) -                     & 
                                                              dia_Feupt(i,k) - diz_Feupt(i,k))*1e3 &
                                             - dtsb * (precip + scaven + feloss + fe2det + fe2poc)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! nanophytoplankton equation !
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             biotr(i,k,id_phy)  = biotr(i,k,id_phy) + dtsb * (f11phy-z21phy-m21phy-f22phy-f23phy)
             if (id_phyfe.ne.0)                                                                    &
             biotr(i,k,id_phyfe)= biotr(i,k,id_phyfe)+dtsb * (phy_Feupt(i,k) - (z21phy + m21phy    &
                                                             + f22phy + f23phy)*phy_FeC(i,k) )

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! microphytoplankton equation !
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             if (id_dia.ne.0)                                                                      &
             biotr(i,k,id_dia)  = biotr(i,k,id_dia) + dtsb * (f11dia-z21dia-m21dia-f22dia-f23dia)
             if (id_diafe.ne.0)                                                                    &
             biotr(i,k,id_diafe)= biotr(i,k,id_diafe)+dtsb * (dia_Feupt(i,k) - (z21dia + m21dia    &
                                                             + f22dia + f23dia)*dia_FeC(i,k) )
             if (id_diasi.ne.0)                                                                    &
             biotr(i,k,id_diasi)= biotr(i,k,id_diasi)+dtsb * ( f11dia * dia_SiCupt(i,k) - (z21dia  &
                                                            + m21dia + f22dia + f23dia)*dia_SiC(i,k))

             !!!!!!!!!!!!!!!!!!!!!!!
             ! diazotroph equation !
             !!!!!!!!!!!!!!!!!!!!!!!
             if (id_diz.ne.0)                                                                      &
             biotr(i,k,id_diz)  = biotr(i,k,id_diz) + dtsb * (f11diz-z21diz-m21diz-f22diz-f23diz)
             if (id_dizfe.ne.0)                                                                    &
             biotr(i,k,id_dizfe)= biotr(i,k,id_dizfe)+dtsb * (diz_Feupt(i,k) - (z21diz + m21diz    &
                                                             + f22diz + f23diz)*diz_FeC(i,k) )

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! microzooplankton equation !
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             biotr(i,k,id_zoo)  = biotr(i,k,id_zoo) + dtsb * (zooassi(i,j) *                       &
                                                             (z21phy + z21dia + z21diz +           &
                                                              z21det + z21poc) -                   &
                                                              z31 - z32 - m21zoo)
             if (id_zoofe.ne.0)                                                                    &
             biotr(i,k,id_zoofe)= biotr(i,k,id_zoofe)+dtsb * (zooassi(i,j)*(z21phy*phy_FeC(i,k) +  &
                                                                            z21dia*dia_FeC(i,k) +  &
                                                                            z21diz*diz_FeC(i,k) +  &
                                                                            z21det*det_FeC(i,k) +  &
                                                                            z21poc*poc_FeC(i,k)) - &
                                                              (z31+z32+m21zoo)*zoo_FeC(i,k) )

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! mesozooplankton equation !
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!
             if (id_mes.ne.0)                                                                      &
             biotr(i,k,id_mes)  = biotr(i,k,id_mes) + dtsb * (zooassi(i,j) *                       &
                                                             (m21phy + m21dia + m21diz +           &
                                                              m21det + m21poc + m21zoo)            &
                                                              - m31 - m32)
             if (id_mesfe.ne.0)                                                                    &
             biotr(i,k,id_mesfe)= biotr(i,k,id_mesfe)+dtsb * (zooassi(i,j)*(m21phy*phy_FeC(i,k) +  &
                                                                            m21dia*dia_FeC(i,k) +  &
                                                                            m21diz*diz_FeC(i,k) +  &
                                                                            m21det*det_FeC(i,k) +  &
                                                                            m21zoo*zoo_FeC(i,k) +  &
                                                                            m21poc*poc_FeC(i,k)) - &
                                                              (m31+m32)*mes_FeC(i,k) )
        

             ! Determine detrital routing to slow or fast sinking detritus
             if (id_poc.ne.0) then
               phyxsize(i,j,k) = ( min(phybiot(i,j),biophy) + 3*max(0.0,biophy-phybiot(i,j)) )     &
                                  / (biophy+epsi)
               diaxsize(i,j,k) = ( min(diabiot(i,j),biodia) + 9*max(0.0,biodia-diabiot(i,j)) )     &
                                  / (biodia+epsi)
               phy_det2poc = min(0.80, 1.0/(phyxsize(i,j,k) + epsi))
               dia_det2poc = min(0.66, 1.0/(diaxsize(i,j,k) + epsi))
               zoo_det2poc = max(0.0, min( 0.7, (phy_det2poc*biophy + dia_det2poc*biodia)          &
                                                / (biophy + biodia + epsi) ))
               mes_det2poc = max(0.0, min( 0.5,                                                    &
                             (phy_det2poc*biophy + dia_det2poc*biodia + zoo_det2poc*biozoo)        &
                                                / (biophy + biodia + biozoo + epsi) ))
             else
               phy_det2poc = 1.0
               dia_det2poc = 1.0
               zoo_det2poc = 1.0
               mes_det2poc = 1.0
             endif
             diz_det2poc = 0.8


             !!!!!!!!!!!!!!!!!!!!!!!!!
             ! Slow-sinking detritus !
             !!!!!!!!!!!!!!!!!!!!!!!!!
             biotr(i,k,id_det)  = biotr(i,k,id_det) + dtsb * ( (1.-zooassi(i,j)) * z21det          &
                                                      +(1.-zooassi(i,j)) * z21phy*phy_det2poc      &
                                                      +(1.-zooassi(i,j)) * z21dia*dia_det2poc      &
                                                      +(1.-zooassi(i,j)) * z21diz*diz_det2poc      &
                                                      +(1.-zooassi(i,j)) * m21det                  &
                                                      +(1.-zooassi(i,j)) * m21phy*phy_det2poc      &
                                                      +(1.-zooassi(i,j)) * m21dia*dia_det2poc      &
                                                      +(1.-zooassi(i,j)) * m21diz*diz_det2poc      &
                                                      +(1.-zooassi(i,j)) * m21zoo*zoo_det2poc      &
                                                      + f23phy * phy_det2poc                       &
                                                      + f23dia * dia_det2poc                       &
                                                      + f23diz * diz_det2poc                       &
                                                      + z32 * zoo_det2poc                          &
                                                      + m32 * mes_det2poc                          &
                                                      - (z21det + m21det + f41det) )
             if (id_detfe.ne.0)                                                                    &
             biotr(i,k,id_detfe)= biotr(i,k,id_detfe) + dtsb * (                                   &
                                               (1.-zooassi(i,j)) * z21det*det_FeC(i,k)             &
                                              +(1.-zooassi(i,j)) * z21phy*phy_det2poc*phy_FeC(i,k) &
                                              +(1.-zooassi(i,j)) * z21dia*dia_det2poc*dia_FeC(i,k) &
                                              +(1.-zooassi(i,j)) * z21diz*diz_det2poc*diz_FeC(i,k) &
                                              +            (1.-zooassi(i,j)) * m21det*det_FeC(i,k) &
                                              +(1.-zooassi(i,j)) * m21phy*phy_det2poc*phy_FeC(i,k) &
                                              +(1.-zooassi(i,j)) * m21dia*dia_det2poc*dia_FeC(i,k) &
                                              +(1.-zooassi(i,j)) * m21diz*diz_det2poc*diz_FeC(i,k) &
                                              +(1.-zooassi(i,j)) * m21zoo*zoo_det2poc*zoo_FeC(i,k) &
                                              + f23phy * phy_FeC(i,k) * phy_det2poc                & 
                                              + f23dia * dia_FeC(i,k) * dia_det2poc                &
                                              + f23diz * diz_FeC(i,k) * diz_det2poc                &
                                              + z32 * zoo_FeC(i,k) * zoo_det2poc                   &
                                              + m32 * mes_FeC(i,k) * mes_det2poc                   &
                                              - (z21det + m21det + f41det)*det_FeC(i,k)            &
                                              + (scadet + fe2det)*1e-3 )

             !!!!!!!!!!!!!!!!!!!!!!!!!
             ! Fast-sinking detritus !
             !!!!!!!!!!!!!!!!!!!!!!!!!
             if (id_poc.ne.0)                                                                      &
             biotr(i,k,id_poc)  = biotr(i,k,id_poc) + dtsb * (                                     &
                                                      (1.-zooassi(i,j)) * z21poc                   &
                                                     +(1.-zooassi(i,j)) * z21phy*(1.-phy_det2poc)  &
                                                     +(1.-zooassi(i,j)) * z21dia*(1.-dia_det2poc)  &
                                                     +(1.-zooassi(i,j)) * z21diz*(1.-diz_det2poc)  &
                                                     +      (1.-zooassi(i,j)) * m21poc             &
                                                     +(1.-zooassi(i,j)) * m21phy*(1.-phy_det2poc)  &
                                                     +(1.-zooassi(i,j)) * m21dia*(1.-dia_det2poc)  &
                                                     +(1.-zooassi(i,j)) * m21diz*(1.-diz_det2poc)  &
                                                     +(1.-zooassi(i,j)) * m21zoo*(1.-zoo_det2poc)  &
                                                     + f23phy * (1.-phy_det2poc)                   &
                                                     + f23dia * (1.-dia_det2poc)                   &
                                                     + f23diz * (1.-diz_det2poc)                   &
                                                     + z32 * (1.-zoo_det2poc)                      &
                                                     + m32 * (1.-mes_det2poc)                      &
                                                     - (z21poc + m21poc + f41poc) )
             if (id_pocfe.ne.0)                                                                    &
             biotr(i,k,id_pocfe)= biotr(i,k,id_pocfe) + dtsb * (                                   &
                                          (1.-zooassi(i,j))*z21poc*poc_FeC(i,k)                    &
                                         +(1.-zooassi(i,j)) * z21phy*(1.-phy_det2poc)*phy_FeC(i,k) &
                                         +(1.-zooassi(i,j)) * z21dia*(1.-dia_det2poc)*dia_FeC(i,k) &
                                         +(1.-zooassi(i,j)) * z21diz*(1.-diz_det2poc)*diz_FeC(i,k) &
                                         +                 (1.-zooassi(i,j)) * m21poc*poc_FeC(i,k) &
                                         +(1.-zooassi(i,j)) * m21phy*(1.-phy_det2poc)*phy_FeC(i,k) &
                                         +(1.-zooassi(i,j)) * m21dia*(1.-dia_det2poc)*dia_FeC(i,k) &
                                         +(1.-zooassi(i,j)) * m21diz*(1.-diz_det2poc)*diz_FeC(i,k) &
                                         +(1.-zooassi(i,j)) * m21zoo*(1.-zoo_det2poc)*zoo_FeC(i,k) &
                                         + f23phy * phy_FeC(i,k) * (1.-phy_det2poc)                &
                                         + f23dia * dia_FeC(i,k) * (1.-dia_det2poc)                &
                                         + f23diz * diz_FeC(i,k) * (1.-diz_det2poc)                &
                                         + z32 * zoo_FeC(i,k) * (1.-zoo_det2poc)                   &
                                         + m32 * mes_FeC(i,k) * (1.-mes_det2poc)                   &
                                         - (z21poc + m21poc + f41poc)*poc_FeC(i,k)                 &
                                         + (scapoc + fe2poc)*1e-3 )
             if (id_pocsi.ne.0)                                                                    &
             biotr(i,k,id_pocsi)= biotr(i,k,id_pocsi)+dtsb * ( (z21dia+m21dia+f23dia)              &
                                                               * dia_SiC(i,k) - f41pocsi )
             
             !!!!!!!!!!!!!!!!!!
             ! CaCO3 equation !
             !!!!!!!!!!!!!!!!!!
             if (id_caco3.ne.0)                                                                    &
             biotr(i,k,id_caco3) = biotr(i,k,id_caco3) + dtsb * ( f_inorg(i,j) * ( f23phy          &
                                   + (1.-zooassi(i,j))*(z21phy + m21phy)                           &
                                   + z32*(biophy*zprefphy(i,j)/(zooprey+epsi))                     &
                                   + m32*(biophy*mprefphy(i,j)/(mesprey+epsi)) ) - f51)   


             !!!!!!!!!!!!!!!!
             ! DIC equation !
             !!!!!!!!!!!!!!!!
             biotr(i,k,id_dic)  = biotr(i,k,id_dic) + dtsb * ( f41det + f41poc + f51               &
                                                             + z31 + m31 + f22phy + f22dia + f22diz&
                                                             - f11phy - f11dia - f11diz            &
                                  - f_inorg(i,j) * ( f23phy + (1.-zooassi(i,j))*(z21phy + m21phy)  &
                                                     + z32*(biophy*zprefphy(i,j)/(zooprey+epsi))   &
                                                     + m32*(biophy*mprefphy(i,j)/(mesprey+epsi)) ) )

             !!!!!!!!!!!!!!!!!!!!!!!!!
             ! Chlorophyll equations !
             !!!!!!!!!!!!!!!!!!!!!!!!!
             if (id_pchl.ne.0)                                                                     &
             biotr(i,k,id_pchl) = biotr(i,k,id_pchl) + dtsb * ( mu_pchl(i,k) - phy_chlc(i,k) *     &
                                                       ( z21phy + m21phy + f22phy + f23phy) * 12.0 )
             if (id_dchl.ne.0)                                                                     &
             biotr(i,k,id_dchl) = biotr(i,k,id_dchl) + dtsb * ( mu_dchl(i,k) - dia_chlc(i,k) *     &
                                                       ( z21dia + m21dia + f22dia + f23dia) * 12.0 )
             if (id_zchl.ne.0)                                                                     &
             biotr(i,k,id_zchl) = biotr(i,k,id_zchl) + dtsb * ( mu_zchl(i,k) - diz_chlc(i,k) *     &
                                                       ( z21diz + m21diz + f22diz + f23diz) * 12.0 )


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.18] - Collect key diagnostics                                                       !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

             ! Gross NPP
             pprod_gross(i,j,k) = pprod_gross(i,j,k) + dtsb*(f11phy+f11dia+f11diz)

             ! Net NPP (gross minus linear mortality (respiration))
	         npp3d(i,j,k) = npp3d(i,j,k) + dtsb*(f11phy+f11dia+f11diz - f22phy-f22dia-f22diz)

             ! Secondary productivity
             zprod_gross(i,j,k) = zprod_gross(i,j,k) + dtsb * zooassi(i,j) *                       &
                                  (z21phy+z21dia+z21dia+z21det+z21poc)
             mprod_gross(i,j,k) = mprod_gross(i,j,k) + dtsb * zooassi(i,j) *                       &
                                  (m21phy+m21dia+m21diz+m21det+m21poc+m21zoo)

             ! Net Secondary production (gross minus respiration losses)
             nsp3d(i,j,k) = nsp3d(i,j,k) + dtsb*((z21phy+z21dia+z21diz+z21det+z21poc +             &
                                                  m21phy+m21dia+m21diz+m21det+m21poc+m21zoo)       &
                                                  * zooassi(i,j) - z31 - m31)

             ! Light limitation of phytoplankton
             phy_parlimit(i,j,k) = phy_parlimit(i,j,k) + phy_lpar(i,k)*dtsb/dtts
             dia_parlimit(i,j,k) = dia_parlimit(i,j,k) + dia_lpar(i,k)*dtsb/dtts
             diz_parlimit(i,j,k) = diz_parlimit(i,j,k) + diz_lpar(i,k)*dtsb/dtts
             ! Iron limitation of phytoplankton
             phy_Felimit(i,j,k) = phy_Felimit(i,j,k) + phy_limfe(i,k)*dtsb/dtts
             dia_Felimit(i,j,k) = dia_Felimit(i,j,k) + dia_limfe(i,k)*dtsb/dtts
             diz_Felimit(i,j,k) = diz_Felimit(i,j,k) + diz_limfe(i,k)*dtsb/dtts
             ! Nitrogen limitation of phytoplankton
             phy_Nlimit(i,j,k) = phy_Nlimit(i,j,k) + phy_limnit(i,k)*dtsb/dtts
             dia_Nlimit(i,j,k) = dia_Nlimit(i,j,k) + dia_limnit(i,k)*dtsb/dtts
             diz_Nlimit(i,j,k) = diz_Nlimit(i,j,k) + diz_limnit(i,k)*dtsb/dtts
             ! Phosphorus limitation of phytoplankton
             phy_Plimit(i,j,k) = phy_Plimit(i,j,k) + phy_limpho(i,k)*dtsb/dtts
             dia_Plimit(i,j,k) = dia_Plimit(i,j,k) + dia_limpho(i,k)*dtsb/dtts
             diz_Plimit(i,j,k) = diz_Plimit(i,j,k) + diz_limpho(i,k)*dtsb/dtts
             ! Silicic acid limitation of microphytoplankton
             dia_Silimit(i,j,k) = dia_Silimit(i,j,k) + dia_limsil(i,k)*dtsb/dtts
             dia_SiCupta(i,j,k) = dia_SiCupta(i,j,k) + dia_SiCupt(i,k)*dtsb/dtts
             
             ! Specific grazing pressure of zooplankton (µM Z per µM P per second)
             if (zooprey.gt.0.0) then
               zoo_grazpres(i,j,k) = zoo_grazpres(i,j,k) + dtsb *                                  &
                                    (z21phy+z21dia+z21diz+z21det+z21poc) /                         & 
                                    (biophy+biodia+biodiz+biodet+biopoc) 
             endif
             if (mesprey.gt.0.0) then
               mes_grazpres(i,j,k) = mes_grazpres(i,j,k) + dtsb *                                  &
                                    (m21phy+m21dia+m21diz+m21det+m21poc+m21zoo) /                  & 
                                    (biophy+biodia+biodiz+biodet+biopoc+biozoo) 
             endif
             
             ! Denitrification and Nitrification
             nitrific(i,j,k) = nitrific(i,j,k) + dtsb * nitrif(i,k)
             denitrif(i,j,k) = denitrif(i,j,k) + dtsb * denitr(i,k)

             ! Nitrogen fixation by diazotrophs 
             diazofix(i,j,k) = diazofix(i,j,k) + dtsb * f11diz * 16./122.0 * (1.0 - diz_limnit(i,k))

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!     Step [11.19] - Check for conservation of mass by ecosystem component                         !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!             

             n_pools(i,k,1) = biotr(i,k,id_nh4) + biotr(i,k,id_no3)                                & 
                              + (biotr(i,k,id_phy) + biotr(i,k,id_dia) + biotr(i,k,id_diz)         &
                              +  biotr(i,k,id_zoo) + biotr(i,k,id_mes)                             &
                              +  biotr(i,k,id_det) + biotr(i,k,id_poc) ) * 16/122.0
             c_pools(i,k,1) = biotr(i,k,id_dic) + biotr(i,k,id_caco3)                              & 
                              + (biotr(i,k,id_phy) + biotr(i,k,id_dia) + biotr(i,k,id_diz)         & 
                              +  biotr(i,k,id_zoo) + biotr(i,k,id_mes)                             &
                              +  biotr(i,k,id_det) + biotr(i,k,id_poc) )

             if (tn.gt.1) then
               !if (abs(n_pools(i,k,2) - n_pools(i,k,1)) .gt. 1e-6 ) then
               !  print*, "Error: Ecosystem model is not conserving nitrogen"
               !  print*, "       Longitude = ", grid%xt(i,j)
               !  print*, "       Latitude  = ", grid%yt(i,j)
               !  print*, "       Depth     = ", grid%zt(k)
               !  print*, " "
               !  print*, "       Biological N budget at two timesteps ", n_pools(i,k,:)
               !  stop
               !endif
               if (abs(c_pools(i,k,2) - c_pools(i,k,1)) .gt. 1e-6 ) then
                 print*, "Error: Ecosystem model is not conserving carbon"
                 print*, "       Longitude = ", grid%xt(i,j)
                 print*, "       Latitude  = ", grid%yt(i,j)
                 print*, "       Depth     = ", grid%zt(k)
                 print*, " "
                 print*, "       Biological C budget at two timesteps ", c_pools(i,k,:)
                 stop
               endif
             endif

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
         
           nh4_bgc_change = grid%tmask(i,j,k) *                                                    &
                            (biotr(i,k,id_nh4) - max(0.0,t_prog(ind_nh4)%field(i,j,k,time%taum1)))
         
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
           if ( id_sil.ne.0) t_prog(ind_sil)%field(i,j,k,time%taum1) = biotr(i,k,id_sil)
           if ( id_nh4.ne.0) t_prog(ind_nh4)%field(i,j,k,time%taum1) = biotr(i,k,id_nh4)
           if ( id_dia.ne.0) t_prog(ind_dia)%field(i,j,k,time%taum1) = biotr(i,k,id_dia)
           if ( id_diz.ne.0) t_prog(ind_diz)%field(i,j,k,time%taum1) = biotr(i,k,id_diz)
           if ( id_mes.ne.0) t_prog(ind_mes)%field(i,j,k,time%taum1) = biotr(i,k,id_mes)
           if ( id_poc.ne.0)   t_prog(ind_poc)%field(i,j,k,time%taum1) = biotr(i,k,id_poc)
           if ( id_caco3.ne.0) t_prog(ind_caco3)%field(i,j,k,time%taum1) = biotr(i,k,id_caco3)
           if ( id_fe.ne.0) t_prog(ind_fe)%field(i,j,k,time%taum1) = biotr(i,k,id_fe)
           if (id_adic.ne.0) &
           t_prog(ind_adic)%field(i,j,k,time%taum1) = t_prog(ind_adic)%field(i,j,k,time%taum1) +   &
                                                      122./16. * (nh4_bgc_change + no3_bgc_change) &
                                                      - caco3_bgc_change 
           if (id_alk.ne.0) &
           t_prog(ind_alk)%field(i,j,k,time%taum1) = t_prog(ind_alk)%field(i,j,k,time%taum1) +     &  
                                                     ( -2.0 * caco3_bgc_change                     &
                                                     + nh4_bgc_change - no3_bgc_change)
           if (id_pchl.ne.0)  t_prog(ind_pchl)%field(i,j,k,time%taum1) = biotr(i,k,id_pchl) 
           if (id_dchl.ne.0)  t_prog(ind_dchl)%field(i,j,k,time%taum1) = biotr(i,k,id_dchl) 
           if (id_zchl.ne.0)  t_prog(ind_zchl)%field(i,j,k,time%taum1) = biotr(i,k,id_zchl) 
           if (id_phyfe.ne.0) t_prog(ind_phyfe)%field(i,j,k,time%taum1) = biotr(i,k,id_phyfe) 
           if (id_diafe.ne.0) t_prog(ind_diafe)%field(i,j,k,time%taum1) = biotr(i,k,id_diafe) 
           if (id_dizfe.ne.0) t_prog(ind_dizfe)%field(i,j,k,time%taum1) = biotr(i,k,id_dizfe) 
           if (id_diasi.ne.0) t_prog(ind_diasi)%field(i,j,k,time%taum1) = biotr(i,k,id_diasi) 
           if (id_zoofe.ne.0) t_prog(ind_zoofe)%field(i,j,k,time%taum1) = biotr(i,k,id_zoofe) 
           if (id_mesfe.ne.0) t_prog(ind_mesfe)%field(i,j,k,time%taum1) = biotr(i,k,id_mesfe) 
           if (id_detfe.ne.0) t_prog(ind_detfe)%field(i,j,k,time%taum1) = biotr(i,k,id_detfe) 
           if (id_pocfe.ne.0) t_prog(ind_pocfe)%field(i,j,k,time%taum1) = biotr(i,k,id_pocfe) 
           if (id_pocsi.ne.0) t_prog(ind_pocsi)%field(i,j,k,time%taum1) = biotr(i,k,id_pocsi)
         

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
           mprod_gross(i,j,k) = rdtts * mprod_gross(i,j,k) * grid%tmask(i,j,k)
           npp3d(i,j,k)       = rdtts * npp3d(i,j,k)       * grid%tmask(i,j,k)
           nsp3d(i,j,k)       = rdtts * nsp3d(i,j,k)       * grid%tmask(i,j,k)
           zoo_grazpres(i,j,k)= rdtts * zoo_grazpres(i,j,k)* grid%tmask(i,j,k)
           mes_grazpres(i,j,k)= rdtts * mes_grazpres(i,j,k)* grid%tmask(i,j,k)

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

       ! initialise advective sinking arrays
       adv_fb1(:,:) = 0.0
       adv_fb2(:,:) = 0.0
       adv_fb3(:,:) = 0.0
       adv_fb4(:,:) = 0.0
       adv_fb5(:,:) = 0.0

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! RASF upstream sinking of detritus and particulates !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do k=2,grid%nk+1
         do i=isc,iec
           adv_fb1(i,k) = wdetbio(i,j)*phyxsize(i,j,k) * biotr(i,k-1,id_det)
           adv_fb2(i,k) = wdetbio(i,j)*phyxsize(i,j,k) * biotr(i,k-1,id_detfe)
           adv_fb3(i,k) = wdetbio(i,j)*10*diaxsize(i,j,k) * biotr(i,k-1,id_poc)
           adv_fb4(i,k) = wdetbio(i,j)*10*diaxsize(i,j,k) * biotr(i,k-1,id_pocfe)
           adv_fb5(i,k) = wdetbio(i,j)*10*diaxsize(i,j,k) * biotr(i,k-1,id_pocsi)
         enddo  ! i
       enddo  ! k

       ! Deposit tracer to sediment as tracer sinks through base of column
       do i = isc,iec 
         k = grid%kmt(i,j)
         if (k.gt.0) then 
           biotic(n)%det_sed_depst(i,j) = adv_fb1(i,k+1) + adv_fb3(i,k+1)
           biotic(n)%detfe_sed_depst(i,j) = adv_fb2(i,k+1) + adv_fb4(i,k+1)
           biotic(n)%detsi_sed_depst(i,j) = adv_fb5(i,k+1)
         endif  ! k.gt.0
       enddo  ! i

       do k =1,grid%nk  !{
         do i =isc,iec  !{
           t_prog(ind_det)%field(i,j,k,time%taum1) = t_prog(ind_det)%field(i,j,k,time%taum1) +     &
               grid%tmask(i,j,k) * dtts * (-adv_fb1(i,k+1) + adv_fb1(i,k))/Thickness%dzt(i,j,k)
           t_prog(ind_detfe)%field(i,j,k,time%taum1) = t_prog(ind_detfe)%field(i,j,k,time%taum1) + &
               grid%tmask(i,j,k) * dtts * (-adv_fb2(i,k+1) + adv_fb2(i,k))/Thickness%dzt(i,j,k)
           t_prog(ind_poc)%field(i,j,k,time%taum1) = t_prog(ind_poc)%field(i,j,k,time%taum1) +     &
               grid%tmask(i,j,k) * dtts * (-adv_fb3(i,k+1) + adv_fb3(i,k))/Thickness%dzt(i,j,k)
           t_prog(ind_pocfe)%field(i,j,k,time%taum1) = t_prog(ind_pocfe)%field(i,j,k,time%taum1) + &
               grid%tmask(i,j,k) * dtts * (-adv_fb4(i,k+1) + adv_fb4(i,k))/Thickness%dzt(i,j,k)
           t_prog(ind_pocsi)%field(i,j,k,time%taum1) = t_prog(ind_pocsi)%field(i,j,k,time%taum1) + &
               grid%tmask(i,j,k) * dtts * (-adv_fb5(i,k+1) + adv_fb5(i,k))/Thickness%dzt(i,j,k)
         enddo  !} i
       enddo  !} k


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

         ! Deposit tracer to sediment as tracer sinks through base of column
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
end subroutine bio_v4
