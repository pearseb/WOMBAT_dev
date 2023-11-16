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
!
!   Purpose
!   -------
!       This code calculates...
!
!
!   Sequence
!   --------
!   - Set up
!       1.  Retrieve key biotic parameters from input file
!       2.  Get time information (timestep length and nested timestep for source and sink equations)
!       3.  Calculate mixed layer depth
!       4.  Initialise diagnostic arrays
!       5.  Retrieve indices of the active tracers
!       6.  Retrieve active tracers
!       7.  Create tracer mask
!
!               [ BEGIN do loop over k and i indices ]
!
!   - Light
!       8.  Calculate short-wave irradiance intensity as a function of depth
!       9.  Include effect of sea-ice concentration
!   - Phytoplankton growth
!       10. Calculate phytoplankton growth rate as a function of light and temperature
!       11. Collect the average light limitition over the mixed layer depth ["light_limit"]
!
!               [ END do loop over k and i indices ]
!
!
!               [ BEGIN do loop over t, k and i indices ]
!
!       12. Apply nutrient limitation terms to maximum phytoplankton growth, taking the
!           minimum of light, no3 and dfe limitations
!   - Zooplankton grazing
!       13. Apply typeIII disk formulation for grazing of phytoplankton biomass
!       14. Calculate a scaling term according to the temperature
!   - Source and sinks
!       15. Collect all the source and sink terms for:
!           i.   Phytoplankton growth and mortality
!           ii.  Zooplankton growth and mortality
!           iii. And how these route to detritus, CaCO3 and inorganic nutrients
!       16. Reduce decay rates of detrius below 180 metres
!       17. Apply source and sink terms to the tracers
!
!               [ END do loop over t, k and i indices ]
!
!
!               [ BEGIN do loop over k and i indices ]
!
!   - Tracer tendencies and diagnostics
!       18. Apply the tracer tendencies to the prognostic tracers in the model
!       19. Collect the diagnostics
!
!               [ END do loop over k and i indices ]
!
!
!               [ BEGIN do loop over k and i indices ]
!
!   - Sinking
!       20. Sink detritus and CaCO3
!       21. Deposit some tracer into the sediment
!       22. Apply change to detritus and CaCO3
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
real, intent(in), dimension(isd:,jsd:,:)                    :: sw_frac_zt        ! short wave fraction on T grid (none)
real, intent(in)                                                :: days_in_this_year

!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------

integer :: i
integer :: j
integer :: k
integer :: n

logical :: used
integer :: index_temp, index_salt

! BGC parameters now read from bgc_param.nc at run time. mac, aug11.  
! include "rjm_param_201005.inc"

      integer :: ts_npzd     ! number of time steps within NPZD model
      integer :: tn, trn
      integer :: ichl

      real :: pi    = 3.14159265358979  !yes, this is pi
      real :: biotr(isc:iec,grid%nk,ntr_bmax), bioma(isc:iec,grid%nk,ntr_bmax), pprod(isc:iec,jsc:jec,grid%nk)
      real :: biono3, biophy, biozoo, biodet, biooxy, biocaco3, biofe 
      real :: biopchl, biophyfe, biozoofe, biodetfe 
      real :: biophy_min, biodet_min
      real :: limno3, limfe, mort_s, lmor_s
      real :: u_npz, gmax, g_npz, prey
      real :: fx1,fx2,fx3,fx4,fu1,fu2
      real :: mumax_phy(isc:iec,grid%nk)
      real :: mu_phy(isc:iec,grid%nk)
!chd auxiliary variables to prevent unnecessary computation
      real :: fbc
      real :: f11,f21p,f21d,f22,f23,f31,f32,f41,f51 
      real :: epsi      = 1e-15
      real :: rdtts      !1/dtts
      real :: dtsb
      real :: adv_fb(isc:iec,1:grid%nk+1)
      real, dimension(isd:ied,jsd:jed) :: mld
      real :: caco3_bgc_change, no3_bgc_change
! PJB additions
      real :: day_angle, declination, cos_hour_angle, daylight_hours, zchl
      real, dimension(4,61) :: zrgb
      real, dimension(isc:iec)           :: keuphot
      real, dimension(isc:iec,grid%nk,3) :: ek_rgb, par_rgb
      real, dimension(isc:iec,grid%nk)   :: par_tot, par_phy, par_phymld, par_eup
      real, dimension(isc:iec,grid%nk)   :: phy_pisl, phy_lday, chl_lday, phy_lpar, chl_lpar
      real, dimension(isc:iec,grid%nk)   :: phy_chlc, phy_chl_lpar, phy_chl_prod, mu_chl
      real, dimension(isc:iec,grid%nk)   :: phy_FeC, zoo_FeC, det_FeC
      real, dimension(isc:iec) :: par_phy_mldsum, par_tot_mldsum, par_z_mldsum
      real :: z1_dep, zval
      real :: phy_pisl2
      real :: phy_chl_pisl, max_chl, mumin_chl


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

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! read biotic parameters !
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  call time_interp_external(alphabio_id, time%model_time, alphabio)
  call time_interp_external(parbio_id, time%model_time, parbio)
  call time_interp_external(kwbio_id, time%model_time, kwbio)
  call time_interp_external(kcbio_id, time%model_time, kcbio)
  call time_interp_external(abio_id, time%model_time, abio)
  call time_interp_external(bbio_id, time%model_time, bbio)
  call time_interp_external(cbio_id, time%model_time, cbio)
  call time_interp_external(k1bio_id, time%model_time, k1bio)
  call time_interp_external(muepbio_id, time%model_time, muepbio)
  call time_interp_external(muepsbio_id, time%model_time, muepsbio)
  call time_interp_external(gam1bio_id, time%model_time, gam1bio)
  call time_interp_external(gbio_id, time%model_time, gbio)
  call time_interp_external(epsbio_id, time%model_time, epsbio)
  call time_interp_external(muezbio_id, time%model_time, muezbio)
  call time_interp_external(gam2bio_id, time%model_time, gam2bio)
  call time_interp_external(muedbio_id, time%model_time, muedbio)
  call time_interp_external(wdetbio_id, time%model_time, wdetbio)
  call time_interp_external(muecaco3_id, time%model_time, muecaco3)
  call time_interp_external(wcaco3_id, time%model_time, wcaco3)
  call time_interp_external(tscav_fe_id, time%model_time, tscav_fe)
  call time_interp_external(fe_bkgnd_id, time%model_time, fe_bkgnd)
  call time_interp_external(f_inorg_id, time%model_time, f_inorg)



  !!!!!!!!!!!!!!!!!!!!!!!!
  ! get basic parameters !
  !!!!!!!!!!!!!!!!!!!!!!!!

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
  call calc_mixed_layer_depth(Thickness,                                                      &
                              T_prog(index_salt)%field(isd:ied,jsd:jed,:,Time%tau),           &
                              T_prog(index_temp)%field(isd:ied,jsd:jed,:,Time%tau),           &
                              Dens%rho(isd:ied,jsd:jed,:,Time%tau),                           &
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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialise diagnostic arrays for saving later !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pprod_gross(:,:,:) = 0.0
  zprod_gross(:,:,:) = 0.0
  light_limit(:,:) = 0.0
  radbio3d(:,:,:) = 0.0
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
  zeuphot(:,:) = 20.0       ! minimum euphotic zone depth is 20 metres
  phy_parlimit(:,:,:) = 0.0
  phy_chl2c(:,:,:) = 0.0
  zoo_grazpres(:,:,:) = 0.0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! possible tracers in model !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! transfer tracers at tau-1 to temporary arrays !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do k=1,grid%nk
          do i=isc,iec

            biotr(i,k,id_no3) = max(0.0,t_prog(ind_no3)%field(i,j,k,Time%taum1)) 
            biotr(i,k,id_phy) = max(0.0,t_prog(ind_phy)%field(i,j,k,Time%taum1))
            biotr(i,k,id_zoo) = max(0.0,t_prog(ind_zoo)%field(i,j,k,Time%taum1))
            biotr(i,k,id_det) = max(0.0,t_prog(ind_det)%field(i,j,k,Time%taum1))
            biotr(i,k,id_o2)  = max(0.0,t_prog(ind_o2)%field(i,j,k,Time%taum1)) 
            if (id_caco3.ne.0) biotr(i,k,id_caco3) = max(0.0,t_prog(ind_caco3)%field(i,j,k,Time%taum1)) 
            if (id_fe.ne.0)    biotr(i,k,id_fe) = max(0.0,t_prog(ind_fe)%field(i,j,k,Time%taum1)) 
            if (id_pchl.ne.0)  biotr(i,k,id_pchl) = max(0.0,t_prog(ind_pchl)%field(i,j,k,Time%taum1))
            if (id_phyfe.ne.0) biotr(i,k,id_phyfe) = max(0.0,t_prog(ind_phyfe)%field(i,j,k,Time%taum1))
            if (id_zoofe.ne.0) biotr(i,k,id_zoofe) = max(0.0,t_prog(ind_zoofe)%field(i,j,k,Time%taum1))
            if (id_detfe.ne.0) biotr(i,k,id_detfe) = max(0.0,t_prog(ind_detfe)%field(i,j,k,Time%taum1))

          enddo  ! i
        enddo  ! k


! pjb - I don't know why you would do this for every tracer... Just use one? One we trust, like DIC?
!       Because otherwise we could be masking places where tracers go negative, like O2 in the anoxic zones

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! create mask bioma, which is 0 for tracer concentration < epsi and 1 when > epsi !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do trn=1,ntr_bmax   
        do k=1,grid%nk
          do i=isc,iec

           if (biotr(i,k,trn) .gt. epsi) then
            bioma(i,k,trn) = 1.
           else
            bioma(i,k,trn) = 0.
           endif 

          enddo
        enddo
       enddo


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Light attenuation (Red, Green and Blue) and Euphotic zone depth !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
       keuphot(:) = 6  ! set k level of euphotic zone to 20 metres

       do k=1,grid%nk
         do i=isc,iec

! pjb - Calculate the daylength (hours per day) given the latitude and day of year
           day_angle = 2.0*pi*(days_in_this_year-1)/365.0
           declination = 0.409 * sin(2*pi*days_in_this_year/365 - 1.39)
           cos_hour_angle = max(min(-tan(declination) * tan(grid%yt(i,j) * pi/180.0), 1.0), -1.0)
           daylight_hours = 2 * (acos(cos_hour_angle) * 180 / pi) / 15.0

! pjb - Get the attenuation coefficients for Red, Green and Blue light given chlorophyll concentration
           zchl = MAX(0.05, MIN(10.0, biotr(i,k,id_pchl) ))
           ichl = NINT( 41 + 20.0*LOG10(zchl) + epsi )
           ek_rgb(i,k,1) = zrgb(2,ichl) * thickness%dzt(i,j,k)
           ek_rgb(i,k,2) = zrgb(3,ichl) * thickness%dzt(i,j,k)
           ek_rgb(i,k,3) = zrgb(4,ichl) * thickness%dzt(i,j,k)

! pjb - Estimate the amount of RGB light available in the water column given RGB attenuation
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

! pjb - Calculate the amount of light available for phytoplankton
           if (swflx(i,j).gt.0.0) then
             par_tot(i,k) =        par_rgb(i,k,1) +        par_rgb(i,k,2) +        par_rgb(i,k,3)
             par_phy(i,k) = 1.85 * par_rgb(i,k,1) + 0.68 * par_rgb(i,k,2) + 0.46 * par_rgb(i,k,3)
           else
             par_tot(i,k) = 0.0
             par_phy(i,k) = 0.0
           endif
           par_eup(i,k) = par_rgb(i,k,1) + par_rgb(i,k,2) + par_rgb(i,k,3)

! pjb - Calculate the euphotic zone depth
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

! pjb - Calculate the integrated light level in the mixed layer
           if (grid%zw(k) .le. mld(i,j)) then
             par_tot_mldsum(i) = par_tot_mldsum(i) + par_tot(i,k) * thickness%dzt(i,j,k)
             par_phy_mldsum(i) = par_phy_mldsum(i) + par_phy(i,k) * thickness%dzt(i,j,k)
             par_z_mldsum(i) = par_z_mldsum(i) + thickness%dzt(i,j,k)
           endif

         enddo  ! i
       enddo  ! k

       mumax_phy(:,:) = 0.0     ! phytoplankton maximum growth rate (temperature-dependent)

       do k=1,grid%nk
         do i=isc,iec

! pjb - Calculate average light level in the mixed layer
           if (grid%zw(k) .le. mld(i,j)) then
             z1_dep = 1.0/par_z_mldsum(i)
             par_phymld(i,k) = par_phy_mldsum(i) * z1_dep
           else
             par_phymld(i,k) = par_phy(i,k)
           endif

! pjb - Calculate impact of daylength on phytoplankton and chlorophyll production
           zval = MAX(1.0, daylight_hours)
           if (grid%zw(k) .le. mld(i,j)) then
             zval = zval * MIN(1.0, (zeuphot(i,j)+epsi)/(mld(i,j)+epsi))
           endif
           chl_lday(i,k) = zval / 24.0
           phy_lday(i,k) = 1.5 * zval / (12.0 + zval)

           ! Save total PAR to radbio array for diagnostic output
           radbio3d(i,j,k) = par_tot(i,k)

           ! Temperature-dependent maximum growth rate (Eppley curve)
           zval = t_prog(index_temp)%field(i,j,k,time%tau)
           mumax_phy(i,k) = abio(i,j) * bbio(i,j)**(cbio(i,j) * zval)

         enddo  ! i
       enddo  ! k

      
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Begin nested time-stepping for ecosystem model !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! Time stepping over dtts with ecosystem model
       !  
       !    (P: phytoplankton, Z: Zooplankton, N: Nitrate and D: Detritus)
       !
       !    dP/dt = u(Temp,Light(Chl),N,Fe) P - p_P P - g(P) P Z
       !
       !    dZ/dt = a g(P) P Z - d Z - p_Z Z^2
       !
       !    dN/dt = r D + d Z - u(Temp,Light(Chl),N,Fe) P  [ + r_d DOC ]
       !
       !    dD/dt = (1-s)[ (1-a) g(P) P Z + p_P P + p_Z Z^2] -r D + w_D dD/dz


       do tn = 1,ts_npzd

        ! Initialise some arrays to zero
        phy_chlc(:,:) = 0.0      ! chlorophyll to carbon ratio of phytoplankton
        phy_FeC(:,:) = 0.0       ! iron to carbon ratio of phytoplankton
        zoo_FeC(:,:) = 0.0       ! iron to carbon ratio of phytoplankton
        det_FeC(:,:) = 0.0       ! iron to carbon ratio of phytoplankton
        phy_lpar(:,:) = 0.0      ! light limitation of phytoplankton growth
        mu_phy(:,:) = 0.0        ! phytoplankton realised growth rate (light limited)
        phy_chl_lpar(:,:) = 0.0  ! light limitation of chlorophyll production by phytoplankton
        phy_chl_prod(:,:) = 0.0  ! phytoplankton chlorophyll realised growth rate
        mu_chl(:,:) = 0.0

        do k=1,grid%nk
         do i=isc,iec

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Get the tracers and save to a scalar variable !
          !   - Remember, tracer values are updated at    !
          !     each of the nested timesteps of the       !
          !     ecosystem model.                          !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          biono3 = max(0.0,biotr(i,k,id_no3))
          biophy = max(0.0,biotr(i,k,id_phy))
          biozoo = max(0.0,biotr(i,k,id_zoo))
          biodet = max(0.0,biotr(i,k,id_det))
          biooxy = max(0.0,biotr(i,k,id_o2))
          if (id_caco3.ne.0) biocaco3 = max(0.0,biotr(i,k,id_caco3))
          if (id_fe.ne.0)    biofe = max(0.0,biotr(i,k,id_fe))
          if (id_pchl.ne.0)  biopchl = max(0.0,biotr(i,k,id_pchl))
          if (id_phyfe.ne.0) biophyfe = max(0.0,biotr(i,k,id_phyfe))
          if (id_zoofe.ne.0) biozoofe = max(0.0,biotr(i,k,id_zoofe))
          if (id_detfe.ne.0) biodetfe = max(0.0,biotr(i,k,id_detfe))

          ! Get in situ temperature
          zval = t_prog(index_temp)%field(i,j,k,time%tau)


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Limitation terms for resource availability !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          limno3 = biono3 / (k1bio(i,j) + biono3)
          if (id_fe.ne.0) then 
            limfe  = biofe / (0.1 + biofe)
          else
            limfe = 1.0
          endif


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Limitation term for light availability !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          !   1. Get chlorophyll:C ratio
          if (biophy*biopchl.gt.0.0) phy_chlc(i,k) = biopchl / (biophy * 12)
          !   2. Get initial slope of Photosynthesis-Irradiance curve
          phy_pisl(i,k) = MAX(alphabio(i,j) * phy_chlc(i,k), alphabio(i,j) * 0.004)
          !   3. Alter the slope to account for respiration and daylength limitation
          phy_pisl2= phy_pisl(i,k) / ( (1.0 + muepbio(i,j)*86400.0) * phy_lday(i,k) )
          !   4. Calculate light limitation
          phy_lpar(i,k) = 1.0 - EXP(-phy_pisl2 * par_phy(i,k))
          !   5. Apply to maximum growth rate
          mu_phy(i,k) = mumax_phy(i,k) * phy_lpar(i,k)

          !if (i.eq.152) then
          !  if (j.eq.112) then
          !    if (k.eq.1) print*, "Grid =", grid%xt(i,j), grid%yt(i,j), grid%zt(k)
          !    print*, phy_chlc(i,k), mumax_phy(i,k), phy_lpar(i,k), mu_phy(i,k)
          !  endif
          !endif


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Realised growth rate given resource availability !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          u_npz = min(limno3, limfe) * mu_phy(i,k)


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Calculate growth in chlorophyll (mg Chl / m3) !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          !     1. Light-limitation of chlorophyll production
          phy_chl_pisl = phy_pisl(i,k) / (mumax_phy(i,k) * 86400.0 * chl_lday(i,k)) 
          phy_chl_lpar(i,k) = 1.0 - EXP( -phy_chl_pisl * par_phymld(i,k))
          phy_chl_prod(i,k) = mumax_phy(i,k) * phy_chl_lpar(i,k)
          !     2. Minimum rate of chlorophyll growth
          mumin_chl = 0.004 * u_npz * biophy * 12. 
          !     3. Maximum chlorophyll:C ratio based on temperature
          max_chl = MIN(0.033, (0.033 / (1.0 - 1.14 / 43.4 * zval)) * (1.0 - 1.14 / 43.4 * 20.0))
          !     4. Calculate mg Chl produced per mg C Phy 
          mu_chl(i,k) = u_npz * biophy * 12. * 86400.0 * (phy_chl_prod(i,k) * min(limno3, limfe)) 
          if (mu_chl(i,k).gt.0.0) then
            mu_chl(i,k) = mumin_chl + ( (max_chl - 0.004) * mu_chl(i,k) )                          &
                                      / ( phy_pisl(i,k) * par_phymld(i,k)/chl_lday(i,k) )
          endif


          !!!!!!!!!!!!!!!!!!!!
          ! Grazing function !
          !!!!!!!!!!!!!!!!!!!!

          biophy_min = MAX(0.0, biophy - 0.001)  ! Minimum Phy beneath which grazing doesn't happen
          biodet_min = MAX(0.0, biodet - 0.001)  ! Minimum Det beneath which grazing doesn't happen
          prey = biophy_min + biodet_min

          gmax = gbio(i,j) * bbio(i,j)**(cbio(i,j) * zval)
!          g_npz = gmax * epsbio(i,j) * biophy*biophy_min /                                         &
!                  ( gmax + epsbio(i,j)*  biophy*biophy_min )
          g_npz = gmax * 0.1/86400.0 * prey*prey / ( gmax + 0.1/86400.0 * prey*prey )


          !!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Temperature dependency !
          !!!!!!!!!!!!!!!!!!!!!!!!!!

          fbc = bbio(i,j)**(cbio(i,j) * zval)


          !!!!!!!!!!!!!!!!!!!!!
          ! Mortality scalers !
          !!!!!!!!!!!!!!!!!!!!!

          ! Depress mortality in nutrient limited regime. Assumes investment in cell maintenance over growth.
          mort_s = MIN(1.0, MAX(0.0, limno3/0.3))
          ! Depress mortality as phytoplankton population dilutes. Assumes density-dependent transfer
          lmor_s = biophy_min / (muepsbio(i,j)*86400 + biophy)

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Collect biomass gain and loss terms !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          f11 = u_npz * biophy                                 * bioma(i,k,id_no3)    ! DIC --> PHY
          if (prey.gt.0.0) then
            f21p = g_npz * biozoo * (biophy_min/prey)          * bioma(i,k,id_phy)            ! PHY --> ZOO & DET
            f21d = g_npz * biozoo * (biodet_min/prey)          * bioma(i,k,id_det)            ! DET --> ZOO & DET
          else
            f21p = 0.0 
            f21d = 0.0
          endif
          !f22 = muepbio(i,j) * lmor_s * mort_s * biophy * fbc * bioma(i,k,id_phy)    ! PHY --> DIC
          f22 = 0.04/86400.0 * mumax_phy(i,k) * mort_s * biophy* bioma(i,k,id_phy)    ! PHY --> DIC
          f23 = 0.5/86400.0 * lmor_s * biophy                  * bioma(i,k,id_phy)    ! PHY --> DET
          f31 = gam2bio(i,j) * biozoo * fbc                    * bioma(i,k,id_zoo)    ! ZOO --> DIC
          f32 = muezbio(i,j) * biozoo*biozoo                   * bioma(i,k,id_zoo)    ! ZOO --> DET
          f41 = muedbio(i,j) * biodet * fbc                    * bioma(i,k,id_det)    ! DET --> DIC

          ! Decrease the remineralisation rate beneath a certain depth
          !     ( This will be irrelevant once explicit bacteria is in )
          if (grid%zw(k) .ge. 180) f41 = f41*.5   ! reduce decay below 180m

          if (id_caco3.ne.0) f51 = muecaco3(i,j)*biocaco3 *bioma(i,k,id_caco3)


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Apply biomass gain and loss terms !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! oxygen equation
          biotr(i,k,id_o2)   = biotr(i,k,id_o2)  - bioma(i,k,id_o2) * 172./122.                    &
                                                 * dtsb * (f41 + f31 + f22 - f11)

          ! nitrate equation
          biotr(i,k,id_no3)  = biotr(i,k,id_no3) + dtsb * (f41 + f31 + f22 - f11) * 16./122. 

          if (id_fe.ne.0)                                                                          &
          biotr(i,k,id_fe)   = biotr(i,k,id_fe)  + dtsb * 2.623e-3 * (f41 + f31 + f22 - f11)

          ! phytoplankton equation
          biotr(i,k,id_phy)  = biotr(i,k,id_phy) + dtsb * (f11 - f21p - f22 - f23)

          ! zooplankton equation
          biotr(i,k,id_zoo)  = biotr(i,k,id_zoo) + dtsb * (gam1bio(i,j)*(f21p+f21d) - f31 - f32)

          ! detritus equation
          biotr(i,k,id_det)  = biotr(i,k,id_det) + dtsb * ((1.-gam1bio(i,j))*(f21p+f21d)           &
                                                          - f21d + f23 + f32 - f41)

          ! caco3 equation
          if (id_caco3.ne.0)                                                                       &
          biotr(i,k,id_caco3) = biotr(i,k,id_caco3) + dtsb * ( (   &
                                (1.-gam1bio(i,j))*f21p + f23 + f32 ) * f_inorg(i,j) - f51)   

          ! phytoplankton chlorophyll equation
          if (id_pchl.ne.0)                                                                        &
          biotr(i,k,id_pchl) = biotr(i,k,id_pchl) + dtsb * ( mu_chl(i,k) - phy_chlc(i,k) *         &
                                                             (f21p + f22 + f23) * 12.0 )


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Estimate primary and secondary production !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! Gross NPP
          pprod_gross(i,j,k) = pprod_gross(i,j,k) + dtsb*f11

          ! Net NPP (gross minus linear mortality (respiration))
	      npp3d(i,j,k) = npp3d(i,j,k) + dtsb*(f11 - f22)

          ! Secondary productivity
          zprod_gross(i,j,k) = zprod_gross(i,j,k) + dtsb*gam1bio(i,j)*(f21p+f21d)

          ! Net Secondary production (gross minus respiration losses)
          nsp3d(i,j,k) = nsp3d(i,j,k) + dtsb*((f21p+f21d)*gam1bio(i,j) - f31)

          ! Light limitation of phytoplankton
          phy_parlimit(i,j,k) = phy_parlimit(i,j,k) + phy_lpar(i,k)*dtsb/dtts
          
          ! Chlorophyll : C ratio of phytoplankton
          phy_chl2c(i,j,k) = phy_chl2c(i,j,k) + phy_chlc(i,k)*dtsb/dtts
          ! Fe : C ratio of phytoplankton
          phy_Fe2C(i,j,k) = phy_Fe2C(i,j,k) + phy_FeC(i,k)*dtsb/dtts
          ! Fe : C ratio of zooplankton
          zoo_Fe2C(i,j,k) = zoo_Fe2C(i,j,k) + zoo_FeC(i,k)*dtsb/dtts
          ! Fe : C ratio of detritus 
          det_Fe2C(i,j,k) = det_Fe2C(i,j,k) + det_FeC(i,k)*dtsb/dtts
          
          ! Specific grazing pressure of zooplankton (µM Z per µM P per second)
          if (biophy*biodet.gt.0.0) then
            zoo_grazpres(i,j,k) = zoo_grazpres(i,j,k) + (f21p+f21d)/(biophy+biodet) * dtsb/dtts
          endif

         enddo  ! i
        enddo  ! k
       enddo  ! tn 


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! add biotically induced tendency to biotracers !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do k = 1,grid%nk  
        do i = isc,iec 

         ! Collect total change in NO3, but add back NO3 if negative at tau-1
         no3_bgc_change = grid%tmask(i,j,k) *                                                      &
             (biotr(i,k,id_no3) - max(0.0,t_prog(ind_no3)%field(i,j,k,time%taum1)))
         
         ! Collect total change in CaCO3, but add back CaCO3 if negative at tau-1
         caco3_bgc_change = grid%tmask(i,j,k) *                                                    &
            (biotr(i,k,id_caco3) - max(0.0,t_prog(ind_caco3)%field(i,j,k,time%taum1)))

          
         ! Add sources and sinks of BGC tracers to tracer arrays
         t_prog(ind_o2)%field(i,j,k,time%taum1)  = biotr(i,k,id_o2)
         t_prog(ind_no3)%field(i,j,k,time%taum1) = biotr(i,k,id_no3)
         t_prog(ind_phy)%field(i,j,k,time%taum1) = biotr(i,k,id_phy)
         t_prog(ind_zoo)%field(i,j,k,time%taum1) = biotr(i,k,id_zoo)
         t_prog(ind_det)%field(i,j,k,time%taum1) = biotr(i,k,id_det)
         if ( id_caco3.ne.0)                                                                       &
         t_prog(ind_caco3)%field(i,j,k,time%taum1) = biotr(i,k,id_caco3)
         if ( id_fe.ne.0)                                                                          &
         t_prog(ind_fe)%field(i,j,k,time%taum1) = biotr(i,k,id_fe) - dtts * tscav_fe(i,j) *        &
                                                  max(0.0, (biotr(i,k,id_fe) - fe_bkgnd(i,j)) ) 
         if (id_dic.ne.0) &
         t_prog(ind_dic)%field(i,j,k,time%taum1) = t_prog(ind_dic)%field(i,j,k,time%taum1) +       &
                                                   122./16. * no3_bgc_change - caco3_bgc_change
         if (id_adic.ne.0) &
         t_prog(ind_adic)%field(i,j,k,time%taum1) = t_prog(ind_adic)%field(i,j,k,time%taum1) +     &
                                                    122./16. * no3_bgc_change - caco3_bgc_change 
         if (id_alk.ne.0) &
         t_prog(ind_alk)%field(i,j,k,time%taum1) = t_prog(ind_alk)%field(i,j,k,time%taum1) +       &  
                                                   ( -2.0 * caco3_bgc_change - no3_bgc_change)
         if (id_pchl.ne.0)                                                                         &
         t_prog(ind_pchl)%field(i,j,k,time%taum1) = biotr(i,k,id_pchl) 

         if (id_phyfe.ne.0)                                                                         &
         t_prog(ind_phyfe)%field(i,j,k,time%taum1) = biotr(i,k,id_phyfe) 
         if (id_zoofe.ne.0)                                                                         &
         t_prog(ind_zoofe)%field(i,j,k,time%taum1) = biotr(i,k,id_zoofe) 
         if (id_detfe.ne.0)                                                                         &
         t_prog(ind_detfe)%field(i,j,k,time%taum1) = biotr(i,k,id_detfe) 
         
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


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! RASF upstream sinking of detritus !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do k=2,grid%nk+1
        do i=isc,iec
         adv_fb(i,k) = wdetbio(i,j) * biotr(i,k-1,id_det)
        enddo  ! i
       enddo  ! k

       ! Set the flux coming through the top face of surface box to zero
       do i=isc,iec
        adv_fb(i,1)  = 0.0
       enddo  ! i 

       ! Deposit tracer to sediment as tracer sinks through base of column.  mac, nov12
       do i = isc,iec 
        k = grid%kmt(i,j)
        if (k.gt.0) then 
         biotic(n)%det_sed_depst(i,j) = adv_fb(i,k+1)
        endif  ! k.gt.0
       enddo  ! i


!       !!mac remineralise tracer sinking through the base of the column.
!       do i = isc, iec 
!         k = grid%kmt(i,j)
!         if (k .gt. 0) then 
!
!          t_prog(ind_no3)%th_tendency(i,j,k) = t_prog(ind_no3)%th_tendency(i,j,k) + &
!           grid%tmask(i,j,k)  * Thickness%rho_dzt(i,j,k,time%tau)  &
!           * 16.0/122.0 * adv_fb(i,k+1)/Thickness%dzt(i,j,k)*0.9
!         
!         t_prog(ind_o2)%th_tendency(i,j,k) = t_prog(ind_o2)%th_tendency(i,j,k) - &
!           172.0 / 122.0 * grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,time%tau)  &
!           * adv_fb(i,k+1)/Thickness%dzt(i,j,k)*0.9
!        
!         if (id_fe .ne. 0) then
!          t_prog(ind_fe)%th_tendency(i,j,k) = t_prog(ind_fe)%th_tendency(i,j,k) + &
!           2.623e-3 * grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,time%tau)  &
!           * adv_fb(i,k+1)/Thickness%dzt(i,j,k)*0.9
!         endif
!        
!         if (id_dic .ne. 0) then
!          t_prog(ind_dic)%th_tendency(i,j,k) = t_prog(ind_dic)%th_tendency(i,j,k) + &
!           grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,time%tau)  &
!           * adv_fb(i,k+1)/Thickness%dzt(i,j,k)*0.9
!         endif
!        
!         if (id_adic .ne. 0) then
!          t_prog(ind_adic)%th_tendency(i,j,k) =t_prog(ind_adic)%th_tendency(i,j,k) + &
!           grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,time%tau)  &
!           * adv_fb(i,k+1)/Thickness%dzt(i,j,k)*0.9
!         endif
!        
!         if (id_alk .ne. 0) then
!          t_prog(ind_alk)%th_tendency(i,j,k) = t_prog(ind_alk)%th_tendency(i,j,k) - &
!           16.0/122.0 * grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,time%tau)  &
!           *adv_fb(i,k+1)/Thickness%dzt(i,j,k)*0.9
!         endif
!
!         endif ! k .gt. 0
!       enddo  ! i


       do k =1,grid%nk  !{
        do i =isc,iec  !{
          t_prog(ind_det)%field(i,j,k,time%taum1) = t_prog(ind_det)%field(i,j,k,time%taum1) + &
              grid%tmask(i,j,k) * dtts * &
              (-adv_fb(i,k+1) + adv_fb(i,k))/Thickness%dzt(i,j,k)
        enddo  !} i
       enddo  !} k


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! RASF upstream sinking of caco3 !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if (id_caco3.ne.0) then
        do k=2,grid%nk+1
         do i=isc,iec
          adv_fb(i,k) = wcaco3(i,j) * biotr(i,k-1,id_caco3)
         enddo
        enddo

        ! Set the flux coming through the top face of surface box to zero
        do i=isc,iec
         adv_fb(i,1)  = 0.0
        enddo 

        ! Deposit tracer to sediment as tracer sinks through base of column.  mac, nov12
        do i = isc, iec 
         k = grid%kmt(i,j)
         if (k .gt. 0) then 
          biotic(n)%caco3_sed_depst(i,j) = adv_fb(i,k+1)
         endif ! k .gt. 0
        enddo ! i


!        !!mac remineralise tracer sinking through the base of the column.
!        do i = isc, iec 
!          k = grid%kmt(i,j)
!          if (k .gt. 0) then 
!          
!           if (id_dic .ne. 0) then
!            t_prog(ind_dic)%th_tendency(i,j,k) = t_prog(ind_dic)%th_tendency(i,j,k) + &
!             grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,time%tau)  &
!             * adv_fb(i,k+1)/Thickness%dzt(i,j,k)*0.9
!           endif
!         
!           if (id_adic .ne. 0) then
!            t_prog(ind_adic)%th_tendency(i,j,k) =t_prog(ind_adic)%th_tendency(i,j,k) + &
!             grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,time%tau)  &
!             * adv_fb(i,k+1)/Thickness%dzt(i,j,k)*0.9
!           endif
!         
!           if (id_alk .ne. 0) then
!            t_prog(ind_alk)%th_tendency(i,j,k) = t_prog(ind_alk)%th_tendency(i,j,k) + &
!             2.0 * grid%tmask(i,j,k) * Thickness%rho_dzt(i,j,k,time%tau)  &
!             *adv_fb(i,k+1)/Thickness%dzt(i,j,k)*0.9
!           endif
!
!          endif ! k .gt. 0
!        enddo ! i


        do k =1,grid%nk  !{
         do i =isc,iec  !{
           t_prog(ind_caco3)%field(i,j,k,time%taum1)=t_prog(ind_caco3)%field(i,j,k,time%taum1) +&
               grid%tmask(i,j,k) * dtts * &
               (-adv_fb(i,k+1) + adv_fb(i,k))/Thickness%dzt(i,j,k)
         enddo  ! i
        enddo  ! k

       endif  ! end loop for caco3

  enddo  !} j

enddo  !} n

return
end subroutine bio_v4
