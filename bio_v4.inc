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
      real :: biotr(isc:iec,grid%nk,ntr_bmax), pprod(isc:iec,jsc:jec,grid%nk)
      real :: bionh4, biono3, biophy, biodia, biozoo, biodet, biopoc, biooxy, biocaco3, biomes
      real :: biofe, biopchl, biodchl, biophyfe, biodiafe, biozoofe, biomesfe, biodetfe, biopocfe 
      real :: biophy_min, biodia_min, biodet_min, biopoc_min, biozoo_min,                          &
              zprefphy, zprefdia, zprefdet, zprefpoc,                                              &
              mprefphy, mprefdia, mprefdet, mprefpoc, mprefzoo
      real :: phy_limnh4, phy_limno3, dia_limnh4, dia_limno3,                                      &
              phy_zlimnit, phy_zlimnh4, phy_zlimno3, dia_zlimnit, dia_zlimnh4, dia_zlimno3,        &
              phy_conc1, phy_conc2, dia_conc1, dia_conc2,                                          &
              phy_k_nit, phy_k_fer, phy_k_pho, dia_k_nit, dia_k_fer, dia_k_pho
      real :: phy_slmor, dia_slmor, phy_sqmor, dia_sqmor
      real :: u_npz, u_ndz, gmax, g_npz, g_npm, zooprey, mesprey
      real :: fx1,fx2,fx3,fx4,fu1,fu2
      real :: mumax_phy(isc:iec,grid%nk), mumax_dia(isc:iec,grid%nk),                              &
              mu_phy(isc:iec,grid%nk), mu_dia(isc:iec,grid%nk)
!chd auxiliary variables to prevent unnecessary computation
      real :: fbc
      real :: f11phy,f11dia,                                                                       &
              z21phy,z21dia,z21det,z21poc,                                                         &
              m21phy,m21dia,m21det,m21poc,m21zoo,                                                  & 
              f22phy,f23phy,f22dia,f23dia,                                                         & 
              z31,z32,m31,m32,f41det,f41poc,f51 
      real :: epsi      = 1e-16
      real :: rdtts      !1/dtts
      real :: dtsb
      real :: adv_fb1(isc:iec,1:grid%nk+1), adv_fb2(isc:iec,1:grid%nk+1),                          &
              adv_fb3(isc:iec,1:grid%nk+1), adv_fb4(isc:iec,1:grid%nk+1)
      real, dimension(isd:ied,jsd:jed) :: mld
      real :: caco3_bgc_change, no3_bgc_change, nh4_bgc_change
! PJB additions
      real :: day_angle, declination, cos_hour_angle, daylight_hours, zchl
      real, dimension(4,61) :: zrgb
      real, dimension(isc:iec)           :: keuphot
      real, dimension(isc:iec,grid%nk,3) :: ek_rgb, par_rgb
      real, dimension(isc:iec,grid%nk)   :: biomask
      real, dimension(isc:iec,grid%nk)   :: par_tot, par_phy, par_dia, par_phymld, par_diamld,     &
                                            par_eup, phy_pisl, dia_pisl, phy_lday, dia_lday,       &
                                            chl_lday, phy_lpar, dia_lpar, pchl_lpar, dchl_lpar,    &
                                            phy_chlc, dia_chlc, phy_chl_lpar, dia_chl_lpar,        &
                                            phy_chl_prod, dia_chl_prod, mu_pchl, mu_dchl,          &
                                            phy_FeC, dia_FeC, zoo_FeC, det_FeC, poc_FeC, mes_FeC,  &
                                            phy_Feupt, dia_Feupt, phy_limfe, dia_limfe,            &
                                            phy_limnit, dia_limnit, nitrif
      real, dimension(isc:iec) :: par_phy_mldsum, par_dia_mldsum, par_tot_mldsum, par_z_mldsum
      real :: z1_dep, zval, ztem, zsal, ztemk
      real :: phy_pisl2, dia_pisl2
      real :: phy_chl_pisl, max_pchl, mumin_pchl, phy_minFeC, phy_maxQFe, phy_KFemin, phy_KFe,     &
              phy_Feupt_upreg, phy_Feupt_downreg, phy_xsize, phy_det2poc,                          &
              dia_chl_pisl, max_dchl, mumin_dchl, dia_minFeC, dia_maxQFe, dia_KFemin, dia_KFe,     &
              dia_Feupt_upreg, dia_Feupt_downreg, dia_xsize, dia_det2poc, zoo_det2poc, mes_det2poc
      real :: ligand, fe_keq, fe_III, fe_lig, fe_col,                                              &
              fesol1, fesol2, fesol3, fesol4, fesol5, hp, fe3sol, precip,                          &
              partic, scaven, scadet, scapoc, feloss, fe2det, fe2poc


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
  mprod_gross(:,:,:) = 0.0
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
  zeuphot(:,:) = 10.0       ! minimum euphotic zone depth is 20 metres
  phy_parlimit(:,:,:) = 0.0
  dia_parlimit(:,:,:) = 0.0
  phy_Felimit(:,:,:) = 0.0
  dia_Felimit(:,:,:) = 0.0
  phy_Nlimit(:,:,:) = 0.0
  dia_Nlimit(:,:,:) = 0.0
  phy_chl2c(:,:,:) = 0.0
  dia_chl2c(:,:,:) = 0.0
  zoo_grazpres(:,:,:) = 0.0
  mes_grazpres(:,:,:) = 0.0
  nitrif1(:,:,:) = 0.0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! possible tracers in model !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ind_po4  = biotic(n)%ind_bgc(id_po4)      ! mmol P / m3
  ind_nh4  = biotic(n)%ind_bgc(id_nh4)      ! mmol N / m3
  ind_no3  = biotic(n)%ind_bgc(id_no3)      ! mmol N / m3
  ind_fe   = biotic(n)%ind_bgc(id_fe)       ! µmol Fe / m3
  ind_o2   = biotic(n)%ind_bgc(id_o2)       ! mmol O2 / m3
  ind_dic  = biotic(n)%ind_bgc(id_dic)      ! mmol DIC / m3
  ind_alk  = biotic(n)%ind_bgc(id_alk)      ! mmol Eq / m3
  ind_caco3= biotic(n)%ind_bgc(id_caco3)    ! mmol C / m3
  ind_adic = biotic(n)%ind_bgc(id_adic)     ! mmol DIC / m3
  ind_phy  = biotic(n)%ind_bgc(id_phy)      ! mmol C / m3
  ind_dia  = biotic(n)%ind_bgc(id_dia)      ! mmol C / m3
  ind_zoo  = biotic(n)%ind_bgc(id_zoo)      ! mmol C / m3
  ind_mes  = biotic(n)%ind_bgc(id_mes)      ! mmol C / m3
  ind_det  = biotic(n)%ind_bgc(id_det)      ! mmol C / m3
  ind_poc  = biotic(n)%ind_bgc(id_poc)      ! mmol C / m3
  ind_pchl = biotic(n)%ind_bgc(id_pchl)     ! mg / m3
  ind_dchl = biotic(n)%ind_bgc(id_dchl)     ! mg / m3
  ind_phyfe = biotic(n)%ind_bgc(id_phyfe)   ! mmol Fe / m3
  ind_diafe = biotic(n)%ind_bgc(id_diafe)   ! mmol Fe / m3
  ind_zoofe = biotic(n)%ind_bgc(id_zoofe)   ! mmol Fe / m3
  ind_mesfe = biotic(n)%ind_bgc(id_mesfe)   ! mmol Fe / m3
  ind_detfe = biotic(n)%ind_bgc(id_detfe)   ! mmol Fe / m3
  ind_pocfe = biotic(n)%ind_bgc(id_pocfe)   ! mmol Fe / m3


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! biotic source-sink terms using Euler forward timestepping !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  do j = jsc, jec  !{ 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! transfer tracers at tau-1 to temporary arrays !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do k=1,grid%nk
          do i=isc,iec

            biotr(i,k,id_dic) = max(0.0,t_prog(ind_dic)%field(i,j,k,Time%taum1)) 
            biotr(i,k,id_no3) = max(0.0,t_prog(ind_no3)%field(i,j,k,Time%taum1)) 
            biotr(i,k,id_phy) = max(0.0,t_prog(ind_phy)%field(i,j,k,Time%taum1))
            biotr(i,k,id_zoo) = max(0.0,t_prog(ind_zoo)%field(i,j,k,Time%taum1))
            biotr(i,k,id_det) = max(0.0,t_prog(ind_det)%field(i,j,k,Time%taum1))
            biotr(i,k,id_o2)  = max(0.0,t_prog(ind_o2)%field(i,j,k,Time%taum1)) 
            if (id_nh4.ne.0)   biotr(i,k,id_nh4) = max(0.0,t_prog(ind_nh4)%field(i,j,k,Time%taum1))
            if (id_dia.ne.0)   biotr(i,k,id_dia) = max(0.0,t_prog(ind_dia)%field(i,j,k,Time%taum1))
            if (id_mes.ne.0)   biotr(i,k,id_mes) = max(0.0,t_prog(ind_mes)%field(i,j,k,Time%taum1))
            if (id_poc.ne.0)   biotr(i,k,id_poc) = max(0.0,t_prog(ind_poc)%field(i,j,k,Time%taum1))
            if (id_caco3.ne.0) biotr(i,k,id_caco3) = max(0.0,t_prog(ind_caco3)%field(i,j,k,Time%taum1)) 
            if (id_fe.ne.0)    biotr(i,k,id_fe) = max(0.0,t_prog(ind_fe)%field(i,j,k,Time%taum1)) 
            if (id_pchl.ne.0)  biotr(i,k,id_pchl) = max(0.0,t_prog(ind_pchl)%field(i,j,k,Time%taum1))
            if (id_dchl.ne.0)  biotr(i,k,id_dchl) = max(0.0,t_prog(ind_dchl)%field(i,j,k,Time%taum1))
            if (id_phyfe.ne.0) biotr(i,k,id_phyfe) = max(0.0,t_prog(ind_phyfe)%field(i,j,k,Time%taum1))
            if (id_diafe.ne.0) biotr(i,k,id_diafe) = max(0.0,t_prog(ind_diafe)%field(i,j,k,Time%taum1))
            if (id_zoofe.ne.0) biotr(i,k,id_zoofe) = max(0.0,t_prog(ind_zoofe)%field(i,j,k,Time%taum1))
            if (id_mesfe.ne.0) biotr(i,k,id_mesfe) = max(0.0,t_prog(ind_mesfe)%field(i,j,k,Time%taum1))
            if (id_detfe.ne.0) biotr(i,k,id_detfe) = max(0.0,t_prog(ind_detfe)%field(i,j,k,Time%taum1))
            if (id_pocfe.ne.0) biotr(i,k,id_pocfe) = max(0.0,t_prog(ind_pocfe)%field(i,j,k,Time%taum1))

          enddo  ! i
        enddo  ! k


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! create mask biomask, which is 0 for tracer concentration < epsi and 1 when > epsi !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       biomask(:,:) = 0.0
       do k=1,grid%nk
         do i=isc,iec

          if (biotr(i,k,id_dic).gt.epsi) biomask(i,k) = 1.0

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
       par_dia(:,:) = 0.0
       par_phymld(:,:) = 0.0
       par_diamld(:,:) = 0.0
       par_tot_mldsum = 0.0
       par_phy_mldsum = 0.0
       par_dia_mldsum = 0.0
       par_z_mldsum = 0.0
       phy_lday(:,:) = 1e-2
       dia_lday(:,:) = 1e-2
       chl_lday(:,:) = 1e-2
       keuphot(:) = 4  ! set k level of euphotic zone to 10 metres

       do k=1,grid%nk
         do i=isc,iec

! pjb - Calculate the daylength (hours per day) given the latitude and day of year
           day_angle = 2.0*pi*(days_in_this_year-1)/365.0
           declination = 0.409 * sin(2*pi*days_in_this_year/365 - 1.39)
           cos_hour_angle = max(min(-tan(declination) * tan(grid%yt(i,j) * pi/180.0), 1.0), -1.0)
           daylight_hours = 2 * (acos(cos_hour_angle) * 180 / pi) / 15.0

! pjb - Get the attenuation coefficients for Red, Green and Blue light given chlorophyll concentration
           zchl = MAX(0.05, MIN(10.0, biotr(i,k,id_pchl)+biotr(i,k,id_dchl) ))
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
             par_dia(i,k) = 1.62 * par_rgb(i,k,1) + 0.74 * par_rgb(i,k,2) + 0.63 * par_rgb(i,k,3)
           else
             par_tot(i,k) = 0.0
             par_phy(i,k) = 0.0
             par_dia(i,k) = 0.0
           endif
           par_eup(i,k) = par_rgb(i,k,1) + par_rgb(i,k,2) + par_rgb(i,k,3)

! pjb - Calculate the integrated light level in the mixed layer
           if (grid%zw(k) .le. mld(i,j)) then
             par_tot_mldsum(i) = par_tot_mldsum(i) + par_tot(i,k) * thickness%dzt(i,j,k)
             par_phy_mldsum(i) = par_phy_mldsum(i) + par_phy(i,k) * thickness%dzt(i,j,k)
             par_dia_mldsum(i) = par_dia_mldsum(i) + par_dia(i,k) * thickness%dzt(i,j,k)
             par_z_mldsum(i) = par_z_mldsum(i) + thickness%dzt(i,j,k)
           endif

         enddo  ! i
       enddo  ! k


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

       mumax_phy(:,:) = 0.0     ! phytoplankton maximum growth rate (temperature-dependent)
       mumax_dia(:,:) = 0.0     ! diatom maximum growth rate (temperature-dependent)

       do k=1,grid%nk
         do i=isc,iec

! pjb - Calculate average light level in the mixed layer
           if (grid%zw(k) .le. mld(i,j)) then
             z1_dep = 1.0/par_z_mldsum(i)
             par_phymld(i,k) = par_phy_mldsum(i) * z1_dep
             par_diamld(i,k) = par_dia_mldsum(i) * z1_dep
           else
             par_phymld(i,k) = par_phy(i,k)
             par_diamld(i,k) = par_dia(i,k)
           endif

! pjb - Calculate impact of daylength on phytoplankton and chlorophyll production
           zval = MAX(1.0, daylight_hours)
           if (grid%zw(k) .le. mld(i,j)) then
             zval = zval * MIN(1.0, (zeuphot(i,j)+epsi)/(mld(i,j)+epsi))
           endif
           chl_lday(i,k) = zval / 24.0
           phy_lday(i,k) = 1.5 * zval / (12.0 + zval)
           dia_lday(i,k) = 1.5 * zval / (12.0 + zval)

           ! Save total PAR to radbio array for diagnostic output
           radbio3d(i,j,k) = par_tot(i,k)

           ! Temperature-dependent maximum growth rate (Eppley curve)
           ztem = t_prog(index_temp)%field(i,j,k,time%tau)
           !mumax_phy(i,k) = abio(i,j) * bbio(i,j)**(cbio(i,j) * ztem)
           !mumax_dia(i,k) = abio(i,j) * bbio(i,j)**(cbio(i,j) * ztem)
           mumax_phy(i,k) = 1.0/86400.0 * bbio(i,j)**(cbio(i,j) * ztem)
           mumax_dia(i,k) = 1.0/86400.0 * bbio(i,j)**(cbio(i,j) * ztem)

         enddo  ! i
       enddo  ! k

      
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Begin nested time-stepping for ecosystem model !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do tn = 1,ts_npzd

        ! Initialise some arrays to zero
        phy_chlc(:,:) = 0.0      ! chlorophyll to carbon ratio of phytoplankton
        dia_chlc(:,:) = 0.0      ! chlorophyll to carbon ratio of diatom
        phy_FeC(:,:) = 0.0       ! iron to carbon ratio of phytoplankton
        dia_FeC(:,:) = 0.0       ! iron to carbon ratio of diatom
        zoo_FeC(:,:) = 0.0       ! iron to carbon ratio of zooplankton
        mes_FeC(:,:) = 0.0       ! iron to carbon ratio of mesozooplankton
        det_FeC(:,:) = 0.0       ! iron to carbon ratio of slow sinking detritus
        poc_FeC(:,:) = 0.0       ! iron to carbon ratio of fast sinking detritus
        phy_Feupt(:,:) = 0.0     ! iron uptake rate of phytoplankton
        dia_Feupt(:,:) = 0.0     ! iron uptake rate of diatom
        phy_limfe(:,:) = 0.0     ! iron limitation of phytoplankton
        dia_limfe(:,:) = 0.0     ! iron limitation of diatom
        phy_lpar(:,:) = 0.0      ! light limitation of phytoplankton growth
        dia_lpar(:,:) = 0.0      ! light limitation of diatom growth
        mu_phy(:,:) = 0.0        ! phytoplankton realised growth rate (light limited)
        mu_dia(:,:) = 0.0        ! diatom realised growth rate (light limited)
        phy_chl_lpar(:,:) = 0.0  ! light limitation of chlorophyll production by phytoplankton
        dia_chl_lpar(:,:) = 0.0  ! light limitation of chlorodiall production by diatom
        phy_chl_prod(:,:) = 0.0  ! phytoplankton chlorophyll realised growth rate
        dia_chl_prod(:,:) = 0.0  ! diatom chlorodiall realised growth rate
        mu_pchl(:,:) = 0.0       ! Growth of chlorophyll in phytoplankton
        mu_dchl(:,:) = 0.0       ! Growth of chlorophyll in diatoms
        nitrif(:,:) = 0.0        ! nitrification rate (NH4 --> NO3)

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
          if (id_nh4.ne.0)   bionh4 = max(0.0,biotr(i,k,id_nh4))
          if (id_poc.ne.0)   biopoc = max(0.0,biotr(i,k,id_poc))
          if (id_dia.ne.0)   biodia = max(0.0,biotr(i,k,id_dia))
          if (id_mes.ne.0)   biomes = max(0.0,biotr(i,k,id_mes))
          if (id_caco3.ne.0) biocaco3 = max(0.0,biotr(i,k,id_caco3))
          if (id_fe.ne.0)    biofe = max(0.0,biotr(i,k,id_fe))
          if (id_pchl.ne.0)  biopchl = max(0.0,biotr(i,k,id_pchl))
          if (id_dchl.ne.0)  biodchl = max(0.0,biotr(i,k,id_dchl))
          if (id_phyfe.ne.0) biophyfe = max(0.0,biotr(i,k,id_phyfe))
          if (id_diafe.ne.0) biodiafe = max(0.0,biotr(i,k,id_diafe))
          if (id_zoofe.ne.0) biozoofe = max(0.0,biotr(i,k,id_zoofe))
          if (id_mesfe.ne.0) biomesfe = max(0.0,biotr(i,k,id_mesfe))
          if (id_detfe.ne.0) biodetfe = max(0.0,biotr(i,k,id_detfe))
          if (id_pocfe.ne.0) biopocfe = max(0.0,biotr(i,k,id_pocfe))

          ! Get important ratios
          if (biophy*biopchl.gt.0.0) phy_chlc(i,k) = biopchl / (biophy * 12)
          if (biodia*biodchl.gt.0.0) dia_chlc(i,k) = biodchl / (biodia * 12)
          if (id_fe.ne.0) then 
            if (biophy*biophyfe.gt.0.0) phy_FeC(i,k) = biophyfe/biophy
            if (biodia*biodiafe.gt.0.0) dia_FeC(i,k) = biodiafe/biodia
            if (biozoo*biozoofe.gt.0.0) zoo_FeC(i,k) = biozoofe/biozoo
            if (biomes*biomesfe.gt.0.0) mes_FeC(i,k) = biomesfe/biomes
            if (biodet*biodetfe.gt.0.0) det_FeC(i,k) = biodetfe/biodet
            if (biopoc*biopocfe.gt.0.0) poc_FeC(i,k) = biopocfe/biopoc
          endif

          ! Temperature and salinity
          ztem = t_prog(index_temp)%field(i,j,k,time%tau)
          ztemk = ztem + 273.15
          zsal = t_prog(index_salt)%field(i,j,k,time%tau)


          !!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! Dissolved Fe chemistry !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! 
          !!! THINGS NEEDING TO BE DONE:
          !!!   - Variable ligand concentration rather than constant?
          !!!   - Variable pH value (currently set to constant 7.9)
          !!!

          ! 1. Determine equilibrium fractionation of total dFe into Fe` and L-Fe
          ligand = 0.7  ! nM ligand concentration
          fe_keq = 10**( 17.27 - 1565.7 / ztemk )
          fe_III = ( -( 1. + fe_keq * ligand*1e-9 - fe_keq * biofe*1e-9 )                          &
                     + SQRT( ( 1. + fe_keq * ligand*1e-9 - fe_keq * biofe*1e-9 )**2                &
                             + 4. * biofe*1e-9 * fe_keq) ) / ( 2. * fe_keq ) * 1e9
          fe_lig = max(0.0, biofe - fe_III)

          ! 2. Precipitation of Fe` (creation of nanoparticles)
          zval = 19.924 * zsal / ( 1000. - 1.005 * zsal)
          fesol1 = 10**(-13.486 - 0.1856*zval**0.5 + 0.3073*zval + 5254.0/max(ztemk, 278.15) )
          fesol2 = 10**(2.517 - 0.8885*zval**0.5 + 0.2139*zval - 1320.0/max(ztemk, 278.15) )
          fesol3 = 10**(0.4511 - 0.3305*zval**0.5 - 1996.0/max(ztemk, 278.15) )
          fesol4 = 10**(-0.2965 - 0.7881*zval**0.5 - 4086.0/max(ztemk, 278.15) )
          fesol5 = 10**(4.4466 - 0.8505*zval**0.5 - 7980.0/max(ztemk, 278.15) )
          hp = 10**(-7.9)
          fe3sol = fesol1 * ( hp**3 + fesol2 * hp**2 + fesol3 * hp + fesol4 + fesol5 / hp ) * 1e9
          precip = max(0.0, ( fe_III - fe3sol ) ) * 0.01/86400.0

          ! 3. Scavenging of Fe` (Fe` --> DET & POC)
          partic = (biodet + biopoc + biocaco3)
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

          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Limitation terms for resource availability !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! Calculate the variable half-saturation Ks (change with population biomass)
          phy_conc1 = max(0.0, biophy - 1.0)  ! pjb replace 1.0 with something from bgc_params.nc
          phy_conc2 = biophy - phy_conc1
          phy_k_nit = max( k1bio(i,j)*0.1,                                                         &
                           ( phy_conc2 * k1bio(i,j)*0.1 + phy_conc1 * 3.0 * k1bio(i,j)*0.1 )       &
                            / ( biophy + epsi ) )
          phy_k_fer = max( 0.1,                                                                    &
                           ( phy_conc2 * 0.1 + phy_conc1 * 3.0 * 0.1 )                             &
                           / ( biophy + epsi ) )
          phy_k_pho = max( 0.07,                                                                   &
                           ( phy_conc2 * 0.07 + phy_conc1 * 3.0 * 0.07 )                           &
                           / ( biophy + epsi ) )
          dia_conc1 = max(0.0, biodia - 1.0)  ! pjb replace 1.0 with something from bgc_params.nc
          dia_conc2 = biodia - dia_conc1
          dia_k_nit = max( k1bio(i,j),                                                             &
                           ( dia_conc2 * k1bio(i,j) + dia_conc1 * 3.0 * k1bio(i,j) )               &
                           / ( biodia + epsi ) )
          dia_k_fer = max( 0.1,                                                                    &
                           ( dia_conc2 * 0.1 + dia_conc1 * 3.0 * 0.1 )                             &
                           / ( biodia + epsi ) )
          dia_k_pho = max( 0.07,                                                                   &
                           ( dia_conc2 * 0.07 + dia_conc1 * 3.0 * 0.07 )                           &
                           / ( biodia + epsi ) )


          ! Nitrogen (NH4 & NO3)
          phy_zlimnh4 = bionh4 / (phy_k_nit + bionh4)
          phy_zlimno3 = biono3 / (phy_k_nit + biono3)
          phy_zlimnit = ( bionh4 + biono3 ) / ( phy_k_nit + bionh4 + biono3 )
          phy_limnh4 = phy_zlimnit * 5.0 * phy_zlimnh4 / ( phy_zlimno3 + 5.0 * phy_zlimnh4 + epsi )
          phy_limno3 = phy_zlimnit * phy_zlimno3       / ( phy_zlimno3 + 5.0 * phy_zlimnh4 + epsi )
          phy_limnit(i,k) = phy_limnh4 + phy_limno3

          dia_zlimnh4 = bionh4 / (dia_k_nit + bionh4)
          dia_zlimno3 = biono3 / (dia_k_nit + biono3)
          dia_zlimnit = ( bionh4 + biono3 ) / ( dia_k_nit + bionh4 + biono3 )
          dia_limnh4 = dia_zlimnit * 5.0 * dia_zlimnh4 / ( dia_zlimno3 + 5.0 * dia_zlimnh4 + epsi )
          dia_limno3 = dia_zlimnit * dia_zlimno3       / ( dia_zlimno3 + 5.0 * dia_zlimnh4 + epsi )
          dia_limnit(i,k) = dia_limnh4 + dia_limno3


          ! Iron
          if (id_fe.ne.0) then 
            if (biopchl*biophy.gt.0.0) then
              phy_minFeC = 0.0016 / 55.85 * biopchl/biophy +                                       &
                           1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * phy_limnit(i,k) +          &
                           1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * phy_limno3
            else
              phy_minFeC = 0.0016 / 55.85 * 0.004*12 +                                             &
                           1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * phy_limnit(i,k) +          &
                           1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * phy_limno3
            endif
            if (biodchl*biodia.gt.0.0) then
              dia_minFeC = 0.0016 / 55.85 * biodchl/biodia +                                       &
                           1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * dia_limnit(i,k) +          &
                           1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * dia_limno3
            else
              dia_minFeC = 0.0016 / 55.85 * 0.004*12 +                                             &
                           1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * dia_limnit(i,k) +          &
                           1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * dia_limno3
            endif
            phy_limfe(i,k) = min(1.0, max(0.0, (phy_FeC(i,k) - phy_minFeC)/7e-6 ))
            dia_limfe(i,k) = min(1.0, max(0.0, (dia_FeC(i,k) - dia_minFeC)/7e-6 ))
          else
            phy_limfe(i,k) = 1.0
            dia_limfe(i,k) = 1.0
          endif


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Limitation term for light availability !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          !   1. Get initial slope of Photosynthesis-Irradiance curve
          phy_pisl(i,k) = MAX(alphabio(i,j) * phy_chlc(i,k), alphabio(i,j) * 0.004)
          dia_pisl(i,k) = MAX(alphabio(i,j) * dia_chlc(i,k), alphabio(i,j) * 0.004)
          !   2. Alter the slope to account for respiration and daylength limitation
          !phy_pisl2= phy_pisl(i,k) / ( (1.0 + muepbio(i,j)*86400.0) * phy_lday(i,k) )
          !dia_pisl2= dia_pisl(i,k) / ( (1.0 + muepbio(i,j)*86400.0) * dia_lday(i,k) )
          phy_pisl2= phy_pisl(i,k) / ( (1.0 + 0.025) * phy_lday(i,k) )
          dia_pisl2= dia_pisl(i,k) / ( (1.0 + 0.025) * dia_lday(i,k) )
          !   3. Calculate light limitation
          phy_lpar(i,k) = 1.0 - EXP(-phy_pisl2 * par_phy(i,k))
          dia_lpar(i,k) = 1.0 - EXP(-dia_pisl2 * par_dia(i,k))
          !   4. Apply to maximum growth rate
          mu_phy(i,k) = mumax_phy(i,k) * phy_lpar(i,k)
          mu_dia(i,k) = mumax_dia(i,k) * dia_lpar(i,k)



          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Realised growth rate given resource availability !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          u_npz = min( phy_limnit(i,k), phy_limfe(i,k) ) * mu_phy(i,k)
          u_ndz = min( dia_limnit(i,k), dia_limfe(i,k) ) * mu_dia(i,k)


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Calculate growth in chlorophyll (mg Chl / m3) !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          !     1. Light-limitation of chlorophyll production
          phy_chl_pisl = phy_pisl(i,k) / (mumax_phy(i,k) * 86400.0 * chl_lday(i,k)) 
          phy_chl_lpar(i,k) = 1.0 - EXP( -phy_chl_pisl * par_phymld(i,k))
          phy_chl_prod(i,k) = mumax_phy(i,k) * phy_chl_lpar(i,k)
          dia_chl_pisl = dia_pisl(i,k) / (mumax_dia(i,k) * 86400.0 * chl_lday(i,k)) 
          dia_chl_lpar(i,k) = 1.0 - EXP( -dia_chl_pisl * par_diamld(i,k))
          dia_chl_prod(i,k) = mumax_dia(i,k) * dia_chl_lpar(i,k)
          !     2. Minimum rate of chlorophyll growth
          mumin_pchl = 0.004 * u_npz * biophy * 12. 
          mumin_dchl = 0.004 * u_ndz * biodia * 12. 
          !     3. Maximum chlorophyll:C ratio based on temperature
          max_pchl = MIN(0.033, (0.033 / (1.0 - 1.14 / 43.4 * ztem)) * (1.0 - 1.14 / 43.4 * 20.0))
          max_dchl = MIN(0.050, (0.050 / (1.0 - 1.14 / 43.4 * ztem)) * (1.0 - 1.14 / 43.4 * 20.0))
          !     4. Calculate mg Chl produced per mg C Phy 
          mu_pchl(i,k) = u_npz * biophy * 12. * 86400.0 *                                          &
                         (phy_chl_prod(i,k) * min(phy_limnit(i,k), phy_limfe(i,k))) 
          mu_dchl(i,k) = u_ndz * biodia * 12. * 86400.0 *                                          &
                         (dia_chl_prod(i,k) * min(dia_limnit(i,k), dia_limfe(i,k))) 
          if (mu_pchl(i,k).gt.0.0) then
            mu_pchl(i,k) = mumin_pchl + ( (max_pchl - 0.004) * mu_pchl(i,k) )                      &
                                      / ( phy_pisl(i,k) * par_phymld(i,k)/chl_lday(i,k) )
          endif
          if (mu_dchl(i,k).gt.0.0) then
            mu_dchl(i,k) = mumin_dchl + ( (max_dchl - 0.004) * mu_dchl(i,k) )                      &
                                      / ( dia_pisl(i,k) * par_diamld(i,k)/chl_lday(i,k) )
          endif


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Calculate phytoplankton uptake of Fe (nmol Fe / m3) !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          if (id_fe.ne.0) then 
            if (phy_FeC(i,k).gt.0.0) then 
              !     1. Find the maximum iron content of phytoplankton biomass
              phy_maxQFe = biophy * 40e-6
              !     2. Find iron half-saturation coefficient, increase as phytoplankton biomass increases based
              !        the understanding that more productive systems have larger cells
              phy_KFemin = 1.0      ! replace with input from bgc_params.nc
              if (biophy.gt.0.0) then
                phy_xsize = (min(1.0,biophy)*phy_KFemin + phy_KFemin*3*max(0.0,biophy-1.0))/biophy
                phy_KFe = max(phy_KFemin, phy_xsize) 
              else
                phy_xsize = 1.0
                phy_KFe = phy_KFemin
              endif
              !     3. Ensure that Fe uptake increases or decreases in response to standing quotas
              !        This is like the cells upregulating or downregulating investment in Fe transporters
              phy_Feupt_upreg   = (4.0 - 4.5 * phy_limfe(i,k) / (phy_limfe(i,k) + 0.5)) 
              phy_Feupt_downreg = max(0.0, (1.0 - biophyfe/phy_maxQFe) / abs(1.05 - biophyfe/phy_maxQFe) )
              !     4. Calculate the iron uptake by phytoplankton
              phy_Feupt(i,k) = biophy * mumax_phy(i,k) * 40e-6 *                                   &
                               biofe / (biofe + phy_KFe) * phy_Feupt_downreg * phy_Feupt_upreg
            endif
            if (dia_FeC(i,k).gt.0.0) then 
              !     1. Find the maximum iron content of diatom biomass
              dia_maxQFe = biodia * 40e-6
              !     2. Find iron half-saturation coefficient, increase as diatom biomass increases based
              !        the understanding that more productive systems have larger cells
              dia_KFemin = 1.0      ! replace with input from bgc_params.nc
              if (biodia.gt.0.0) then
                dia_xsize = (min(1.0,biodia)*dia_KFemin + dia_KFemin*3*max(0.0,biodia-1.0))/biodia
                dia_KFe = max(dia_KFemin, dia_xsize) 
              else
                dia_xsize = 1.0
                dia_KFe = dia_KFemin
              endif
              !     3. Ensure that Fe uptake increases or decreases in response to standing quotas
              !        This is like the cells upregulating or downregulating investment in Fe transporters
              dia_Feupt_upreg   = (4.0 - 4.5 * dia_limfe(i,k) / (dia_limfe(i,k) + 0.5)) 
              dia_Feupt_downreg = max(0.0, (1.0 - biodiafe/dia_maxQFe) / abs(1.05 - biodiafe/dia_maxQFe) )
              !     4. Calculate the iron uptake by diatom
              dia_Feupt(i,k) = biodia * mumax_dia(i,k) * 40e-6 *                                   &
                               biofe / (biofe + dia_KFe) * dia_Feupt_downreg * dia_Feupt_upreg
            endif
          endif
          

          !!!!!!!!!!!!!!!!!!!!
          ! Grazing function !
          !!!!!!!!!!!!!!!!!!!!

          biophy_min = MAX(0.0, biophy - 0.001)  ! Minimum Phy beneath which grazing doesn't happen
          biodia_min = MAX(0.0, biodia - 0.001)  ! Minimum dia beneath which grazing doesn't happen
          biodet_min = MAX(0.0, biodet - 0.001)  ! Minimum Det beneath which grazing doesn't happen
          biopoc_min = MAX(0.0, biopoc - 0.001)  ! Minimum poc beneath which grazing doesn't happen
          biozoo_min = MAX(0.0, biozoo - 0.001)  ! Minimum mes beneath which grazing doesn't happen

          ! preferences and prey availability
          zprefphy = 1.00
          zprefdia = 0.50
          zprefdet = 0.25
          zprefpoc = 0.10
          mprefphy = 0.25
          mprefdia = 1.00
          mprefdet = 0.35
          mprefpoc = 0.35
          mprefzoo = 1.00
          zooprey = zprefphy*biophy_min + zprefdia*biodia_min + zprefdia*biodet_min                &
                    + zprefpoc*biopoc_min
          mesprey = mprefphy*biophy_min + mprefdia*biodia_min + mprefdet*biodet_min                &
                    + mprefpoc*biopoc_min + mprefzoo*biozoo_min

!          gmax = gbio(i,j) * bbio(i,j)**(cbio(i,j) * ztem)
          gmax = gbio(i,j) * 1.075**(cbio(i,j) * ztem)
!          g_npz = gmax * epsbio(i,j) * zooprey*zooprey / ( gmax + epsbio(i,j) * zooprey*zooprey )
          g_npz = gmax*0.50 * 1.0/86400.0 * zooprey*zooprey / ( gmax*0.50 + 1.0/86400.0 * zooprey*zooprey )
          g_npm = gmax*0.25 * 0.5/86400.0 * mesprey*mesprey / ( gmax*0.25 + 0.5/86400.0 * mesprey*mesprey )


          !!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Temperature dependency !
          !!!!!!!!!!!!!!!!!!!!!!!!!!

!          fbc = bbio(i,j)**(cbio(i,j) * ztem)
          fbc = 1.075**(cbio(i,j) * ztem)


          !!!!!!!!!!!!!!!!!!!!!
          ! Mortality scalers !
          !!!!!!!!!!!!!!!!!!!!!

          ! Depress mortality in nutrient limited regime. Assumes investment in cell maintenance over growth.
          phy_slmor = MIN(1.0, MAX(0.0, min(phy_limnit(i,k), phy_limfe(i,k))/0.3))
          dia_slmor = MIN(1.0, MAX(0.0, min(dia_limnit(i,k), dia_limfe(i,k))/0.3))
          ! Depress mortality as phytoplankton population dilutes. Assumes density-dependent transfer
!          phy_sqmor = biophy_min / (muepsbio(i,j)*86400 + biophy)
!          dia_sqmor = biodia_min / (muepsbio(i,j)*86400 + biodia)


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Collect biomass gain and loss terms !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! PRIMARY PRODUCTION
          f11phy = u_npz * biophy                              * biomask(i,k)    ! DIC --> PHY
          f11dia = u_ndz * biodia                              * biomask(i,k)    ! DIC --> DIA

          ! GRAZING
          if (zooprey.gt.0.0) then
            z21phy = g_npz * biozoo * (zprefphy*biophy_min/zooprey)* biomask(i,k)    ! PHY --> ZOO & DET & POC
            z21dia = g_npz * biozoo * (zprefdia*biodia_min/zooprey)* biomask(i,k)    ! DIA --> ZOO & DET & POC
            z21det = g_npz * biozoo * (zprefdet*biodet_min/zooprey)* biomask(i,k)    ! DET --> ZOO & DET
            z21poc = g_npz * biozoo * (zprefpoc*biopoc_min/zooprey)* biomask(i,k)    ! POC --> ZOO & POC
            m21phy = g_npm * biomes * (mprefphy*biophy_min/mesprey)* biomask(i,k)    ! PHY --> MES & DET & POC
            m21dia = g_npm * biomes * (mprefdia*biodia_min/mesprey)* biomask(i,k)    ! DIA --> MES & DET & POC
            m21det = g_npm * biomes * (mprefdet*biodet_min/mesprey)* biomask(i,k)    ! DET --> MES & DET
            m21poc = g_npm * biomes * (mprefpoc*biopoc_min/mesprey)* biomask(i,k)    ! POC --> MES & POC
            m21zoo = g_npm * biomes * (mprefzoo*biozoo_min/mesprey)* biomask(i,k)    ! ZOO --> MES & POC
          else
            z21phy = 0.0 
            z21dia = 0.0 
            z21det = 0.0
            z21poc = 0.0
            m21phy = 0.0 
            m21dia = 0.0 
            m21det = 0.0
            m21poc = 0.0
            m21zoo = 0.0
          endif

          ! MORTALITY
          f22phy = 0.025/86400.0 * fbc * phy_slmor * biophy    * biomask(i,k)    ! PHY --> DIC
          f22dia = 0.025/86400.0 * fbc * dia_slmor * biodia    * biomask(i,k)    ! DIA --> DIC
          f23phy = 0.1/86400.0 * biophy * biophy               * biomask(i,k)    ! PHY --> DET & POC
          f23dia = 0.1/86400.0 * biodia * biodia               * biomask(i,k)    ! DIA --> DET & POC
          z31 = gam2bio(i,j)*0.5 * biozoo * fbc                    * biomask(i,k)    ! ZOO --> DIC
          z32 = muezbio(i,j)*0.1 * biozoo*biozoo               * biomask(i,k)    ! ZOO --> DET & POC
          m31 = gam2bio(i,j)*0.5 * biomes * fbc                    * biomask(i,k)    ! MES --> DIC
          m32 = muezbio(i,j) * biomes*biomes                   * biomask(i,k)    ! MES --> DET & POC
          f41det = muedbio(i,j) * biodet * fbc                 * biomask(i,k)    ! DET --> DIC
          f41poc = muedbio(i,j) * biopoc * fbc                 * biomask(i,k)    ! POC --> DIC

          ! Decrease the remineralisation rate beneath a certain depth
          !     ( This will be irrelevant once explicit bacteria is in )
          if (grid%zw(k) .ge. 180) f41det = f41det*.5   ! reduce decay below 180m
          if (grid%zw(k) .ge. 180) f41poc = f41poc*.5   ! reduce decay below 180m

          if (id_caco3.ne.0) f51 = muecaco3(i,j)* biocaco3 * biomask(i,k)


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Chemoautotrophic metabolisms !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! Nitrification (set at maximum of 5% of NH4 per day)
          nitrif(i,k) = 0.05/86400.0 * bionh4


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Apply biomass gain and loss terms !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! oxygen equation
          biotr(i,k,id_o2)   = biotr(i,k,id_o2)  - biomask(i,k) * 172./122. * dtsb * (             &
                                                   f41det + f41poc + z31 + m31 + f22phy + f22dia   &
                                                   - f11phy - f11dia )

          ! ammonium equation
          biotr(i,k,id_nh4)  = biotr(i,k,id_nh4) - dtsb * nitrif(i,k) + dtsb * 16./122. * (        &
                                                   f41det + f41poc + z31 + m31 + f22phy + f22dia   &
                                                 - f11phy *(phy_limnh4+epsi)/(phy_limnit(i,k)+epsi)&
                                                 - f11dia *(dia_limnh4+epsi)/(dia_limnit(i,k)+epsi) )


          ! nitrate equation
          biotr(i,k,id_no3)  = biotr(i,k,id_no3) + dtsb * nitrif(i,k) - dtsb * 16./122. * (        &
                                                   f11phy *(phy_limno3+epsi)/(phy_limnit(i,k)+epsi)&
                                                 + f11dia *(dia_limno3+epsi)/(dia_limnit(i,k)+epsi) )

          if (id_fe.ne.0)                                                                          &
          biotr(i,k,id_fe)   = biotr(i,k,id_fe)  + dtsb * (f41det*det_FeC(i,k) +                   &
                                                           f41poc*poc_FeC(i,k) +                   &
                                                           z31*zoo_FeC(i,k) +                      &
                                                           m31*mes_FeC(i,k) +                      &
                                                           f22phy*phy_FeC(i,k) +                   &
                                                           f22dia*dia_FeC(i,k) -                   &
                                                           phy_Feupt(i,k) - dia_Feupt(i,k))*1e3    &
                                             - dtsb * (precip + scaven + feloss + fe2det + fe2poc)

          ! phytoplankton equation
          biotr(i,k,id_phy)  = biotr(i,k,id_phy) + dtsb * (f11phy - z21phy - m21phy - f22phy - f23phy)
          if (id_phyfe.ne.0)                                                                       &
          biotr(i,k,id_phyfe)= biotr(i,k,id_phyfe)+dtsb * (phy_Feupt(i,k) -                        &
                                                           (z21phy+m21phy+f22phy+f23phy)*phy_FeC(i,k))

          ! diatom equation
          biotr(i,k,id_dia)  = biotr(i,k,id_dia) + dtsb * (f11dia - z21dia - m21dia - f22dia - f23dia)
          if (id_diafe.ne.0)                                                                       &
          biotr(i,k,id_diafe)= biotr(i,k,id_diafe)+dtsb * (dia_Feupt(i,k) -                        &
                                                           (z21dia+m21dia+f22dia+f23dia)*dia_FeC(i,k))


          ! zooplankton equation
          biotr(i,k,id_zoo)  = biotr(i,k,id_zoo) + dtsb * (gam1bio(i,j) *                          &
                                                          (z21phy + z21dia + z21det + z21poc) -    &
                                                           z31 - z32 - m21zoo)
          if (id_zoofe.ne.0)                                                                       &
          biotr(i,k,id_zoofe)= biotr(i,k,id_zoofe)+dtsb * (gam1bio(i,j)*(z21phy*phy_FeC(i,k) +     &
                                                                         z21dia*dia_FeC(i,k) +     &
                                                                         z21det*det_FeC(i,k) +     &
                                                                         z21poc*poc_FeC(i,k)) -    &
                                                           (z31+z32+m21zoo)*zoo_FeC(i,k) )
          
          ! mesozooplankton equation
          biotr(i,k,id_mes)  = biotr(i,k,id_mes) + dtsb * (gam1bio(i,j) *                          &
                                                          (m21phy+m21dia+m21det+m21poc+m21zoo) -   &
                                                           m31 - m32)
          if (id_mesfe.ne.0)                                                                       &
          biotr(i,k,id_mesfe)= biotr(i,k,id_mesfe)+dtsb * (gam1bio(i,j)*(m21phy*phy_FeC(i,k) +     &
                                                                         m21dia*dia_FeC(i,k) +     &
                                                                         m21det*det_FeC(i,k) +     &
                                                                         m21zoo*zoo_FeC(i,k) +     &
                                                                         m21poc*poc_FeC(i,k)) -    &
                                                           (m31+m32)*mes_FeC(i,k) )
          
          if (id_poc.ne.0) then
            phy_det2poc = min(0.80, 1.0/(phy_xsize + epsi))
            dia_det2poc = min(0.66, 1.0/(dia_xsize + epsi))
            zoo_det2poc = max(0.0, min( 0.7, (phy_det2poc*biophy + dia_det2poc*biodia)             &
                                             / (biophy + biodia + epsi) ))
            mes_det2poc = max(0.0, min( 0.5,                                                      &
                          (phy_det2poc*biophy + dia_det2poc*biodia + zoo_det2poc*biozoo)           &
                                             / (biophy + biodia + biozoo + epsi) ))
          else
            phy_det2poc = 1.0
            dia_det2poc = 1.0
            zoo_det2poc = 1.0
            mes_det2poc = 1.0
          endif

          ! detritus equation
          biotr(i,k,id_det)  = biotr(i,k,id_det) + dtsb * ( (1.-gam1bio(i,j)) * z21det             &
                                                           +(1.-gam1bio(i,j)) * z21phy*phy_det2poc &
                                                           +(1.-gam1bio(i,j)) * z21dia*dia_det2poc &
                                                           +(1.-gam1bio(i,j)) * m21det             &
                                                           +(1.-gam1bio(i,j)) * m21phy*phy_det2poc &
                                                           +(1.-gam1bio(i,j)) * m21dia*dia_det2poc &
                                                           +(1.-gam1bio(i,j)) * m21zoo*zoo_det2poc &
                                                           + f23phy * phy_det2poc                  &
                                                           + f23dia * dia_det2poc                  &
                                                           + z32 * zoo_det2poc                     &
                                                           + m32 * mes_det2poc                     &
                                                           - (z21det + m21det + f41det) )
          if (id_detfe.ne.0)                                                                       &
          biotr(i,k,id_detfe)= biotr(i,k,id_detfe)+dtsb * ((1.-gam1bio(i,j)) * z21det*det_FeC(i,k) &
                                              +(1.-gam1bio(i,j)) * z21phy*phy_det2poc*phy_FeC(i,k) &
                                              +(1.-gam1bio(i,j)) * z21dia*dia_det2poc*dia_FeC(i,k) &
                                              +            (1.-gam1bio(i,j)) * m21det*det_FeC(i,k) &
                                              +(1.-gam1bio(i,j)) * m21phy*phy_det2poc*phy_FeC(i,k) &
                                              +(1.-gam1bio(i,j)) * m21dia*dia_det2poc*dia_FeC(i,k) &
                                              +(1.-gam1bio(i,j)) * m21zoo*zoo_det2poc*zoo_FeC(i,k) &
                                              + f23phy * phy_FeC(i,k) * phy_det2poc                & 
                                              + f23dia * dia_FeC(i,k) * dia_det2poc                &
                                              + z32 * zoo_FeC(i,k) * zoo_det2poc                   &
                                              + m32 * mes_FeC(i,k) * mes_det2poc                   &
                                              - (z21det + m21det + f41det)*det_FeC(i,k)            &
                                              + (scadet + fe2det)*1e-3 )

          ! particulates equation
          if (id_poc.ne.0)                                                                         &
          biotr(i,k,id_poc)  = biotr(i,k,id_poc) + dtsb * ( (1.-gam1bio(i,j)) * z21poc             &
                                                     +(1.-gam1bio(i,j)) * z21phy*(1.-phy_det2poc)  &
                                                     +(1.-gam1bio(i,j)) * z21dia*(1.-dia_det2poc)  &
                                                     +      (1.-gam1bio(i,j)) * m21poc             &
                                                     +(1.-gam1bio(i,j)) * m21phy*(1.-phy_det2poc)  &
                                                     +(1.-gam1bio(i,j)) * m21dia*(1.-dia_det2poc)  &
                                                     +(1.-gam1bio(i,j)) * m21zoo*(1.-zoo_det2poc)  &
                                                     + f23phy * (1.-phy_det2poc)                   &
                                                     + f23dia * (1.-dia_det2poc)                   &
                                                     + z32 * (1.-zoo_det2poc)                      &
                                                     + m32 * (1.-mes_det2poc)                      &
                                                     - (z21poc + m21poc + f41poc) )
          if (id_pocfe.ne.0)                                                                       &
          biotr(i,k,id_pocfe)= biotr(i,k,id_pocfe)+dtsb * ( (1.-gam1bio(i,j))*z21poc*poc_FeC(i,k)  &
                                         +(1.-gam1bio(i,j)) * z21phy*(1.-phy_det2poc)*phy_FeC(i,k) &
                                         +(1.-gam1bio(i,j)) * z21dia*(1.-dia_det2poc)*dia_FeC(i,k) &
                                         +                 (1.-gam1bio(i,j)) * m21poc*poc_FeC(i,k) &
                                         +(1.-gam1bio(i,j)) * m21phy*(1.-phy_det2poc)*phy_FeC(i,k) &
                                         +(1.-gam1bio(i,j)) * m21dia*(1.-dia_det2poc)*dia_FeC(i,k) &
                                         +(1.-gam1bio(i,j)) * m21zoo*(1.-zoo_det2poc)*zoo_FeC(i,k) &
                                         + f23phy * phy_FeC(i,k) * (1.-phy_det2poc)                &
                                         + f23dia * dia_FeC(i,k) * (1.-dia_det2poc)                &
                                         + z32 * zoo_FeC(i,k) * (1.-zoo_det2poc)                   &
                                         + m32 * mes_FeC(i,k) * (1.-mes_det2poc)                   &
                                         - (z21poc + m21poc + f41poc)*poc_FeC(i,k)                 &
                                         + (scapoc + fe2poc)*1e-3 )
                     

          !if (i.eq.330) then
          !  if (j.eq.194) then
          !    print*, "Grid k =", k, biofe, biomask(i,k)
          !    print*, biophy, biozoo, biodet
          !    print*, phy_FeC(i,k), zoo_FeC(i,k), det_FeC(i,k)
          !    print*, phy_minFeC, phy_limfe(i,k), phy_Feupt(i,k)
          !    !print*, phy_Feupt(i,k), z21phy*phy_FeC(i,k), f22phy*phy_FeC(i,k), f23phy*phy_FeC(i,k)
          !    !print*, z21det*det_FeC(i,k), z31*zoo_FeC(i,k), z32*zoo_FeC(i,k), f41det*det_FeC(i,k)
          !  endif
          !endif

          !if (i.eq.160) then
          !  if (j.eq.80) then
          !    print*, "Grid =", grid%xt(i,j), grid%yt(i,j), grid%zt(k)
          !    print*, biono3, limno3 
          !    print*, biofe, phy_limfe(i,k)
          !    print*, biophy, biozoo, biodet
          !    print*, ztem, mumax_phy(i,k), mort_s, lmor_s
          !    !print*, phy_FeC(i,k), zoo_FeC(i,k), det_FeC(i,k)
          !    print*, f11phy, z21phy, z21det, f22phy, f23phy
          !  endif
          !endif

          ! caco3 equation
          if (id_caco3.ne.0)                                                                       &
          biotr(i,k,id_caco3) = biotr(i,k,id_caco3) + dtsb * ( (   &
                                (1.-gam1bio(i,j))*z21phy + f23phy + z32 ) * f_inorg(i,j) - f51)   

          ! chlorophyll equations
          if (id_pchl.ne.0)                                                                        &
          biotr(i,k,id_pchl) = biotr(i,k,id_pchl) + dtsb * ( mu_pchl(i,k) - phy_chlc(i,k) *        &
                                                             (z21phy + f22phy + f23phy) * 12.0 )
          if (id_dchl.ne.0)                                                                        &
          biotr(i,k,id_dchl) = biotr(i,k,id_dchl) + dtsb * ( mu_dchl(i,k) - dia_chlc(i,k) *        &
                                                             (z21dia + f22dia + f23dia) * 12.0 )


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Estimate primary and secondary production !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! Gross NPP
          pprod_gross(i,j,k) = pprod_gross(i,j,k) + dtsb*(f11phy+f11dia)

          ! Net NPP (gross minus linear mortality (respiration))
	      npp3d(i,j,k) = npp3d(i,j,k) + dtsb*(f11phy+f11dia - f22phy-f22dia)

          ! Secondary productivity
          zprod_gross(i,j,k) = zprod_gross(i,j,k) + dtsb*gam1bio(i,j)*(z21phy+z21dia+z21det+z21poc)
          mprod_gross(i,j,k) = mprod_gross(i,j,k) + dtsb*gam1bio(i,j)*(m21phy+m21dia+m21det+m21poc+m21zoo)

          ! Net Secondary production (gross minus respiration losses)
          nsp3d(i,j,k) = nsp3d(i,j,k) + dtsb*((z21phy+z21dia+z21det+z21poc +                       &
                                               m21phy+m21dia+m21det+m21poc+m21zoo)*gam1bio(i,j)    &
                                               - z31 - m31)

          ! Light limitation of phytoplankton
          phy_parlimit(i,j,k) = phy_parlimit(i,j,k) + phy_lpar(i,k)*dtsb/dtts
          dia_parlimit(i,j,k) = dia_parlimit(i,j,k) + dia_lpar(i,k)*dtsb/dtts
          ! Iron limitation of phytoplankton
          phy_Felimit(i,j,k) = phy_Felimit(i,j,k) + phy_limfe(i,k)*dtsb/dtts
          dia_Felimit(i,j,k) = dia_Felimit(i,j,k) + dia_limfe(i,k)*dtsb/dtts
          ! Nitrogen limitation of phytoplankton
          phy_Nlimit(i,j,k) = phy_Nlimit(i,j,k) + phy_limnit(i,k)*dtsb/dtts
          dia_Nlimit(i,j,k) = dia_Nlimit(i,j,k) + dia_limnit(i,k)*dtsb/dtts
          
          ! Chlorophyll : C ratio of phytoplankton
          phy_chl2c(i,j,k) = phy_chl2c(i,j,k) + phy_chlc(i,k)*dtsb/dtts
          dia_chl2c(i,j,k) = dia_chl2c(i,j,k) + dia_chlc(i,k)*dtsb/dtts
          
          ! Specific grazing pressure of zooplankton (µM Z per µM P per second)
          if (zooprey.gt.0.0) then
            zoo_grazpres(i,j,k) = zoo_grazpres(i,j,k) + dtsb *                                     &
                                 (z21phy+z21dia+z21det+z21poc) /                                   & 
                                 (biophy+biodia+biodet+biopoc) 
          endif
          if (mesprey.gt.0.0) then
            mes_grazpres(i,j,k) = mes_grazpres(i,j,k) + dtsb *                                     &
                                 (m21phy+m21dia+m21det+m21poc+m21zoo) /                            & 
                                 (biophy+biodia+biodet+biopoc+biozoo) 
          endif
          
          ! Nitrification
          nitrif1(i,j,k) = nitrif1(i,j,k) + dtsb * nitrif(i,k)


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
         
         nh4_bgc_change = grid%tmask(i,j,k) *                                                      &
             (biotr(i,k,id_nh4) - max(0.0,t_prog(ind_nh4)%field(i,j,k,time%taum1)))
         
         ! Collect total change in CaCO3, but add back CaCO3 if negative at tau-1
         caco3_bgc_change = grid%tmask(i,j,k) *                                                    &
            (biotr(i,k,id_caco3) - max(0.0,t_prog(ind_caco3)%field(i,j,k,time%taum1)))

          
         ! Add sources and sinks of BGC tracers to tracer arrays
         t_prog(ind_o2)%field(i,j,k,time%taum1)  = biotr(i,k,id_o2)
         t_prog(ind_no3)%field(i,j,k,time%taum1) = biotr(i,k,id_no3)
         t_prog(ind_phy)%field(i,j,k,time%taum1) = biotr(i,k,id_phy)
         t_prog(ind_zoo)%field(i,j,k,time%taum1) = biotr(i,k,id_zoo)
         t_prog(ind_det)%field(i,j,k,time%taum1) = biotr(i,k,id_det)
         if ( id_nh4.ne.0) t_prog(ind_nh4)%field(i,j,k,time%taum1) = biotr(i,k,id_nh4)
         if ( id_dia.ne.0) t_prog(ind_dia)%field(i,j,k,time%taum1) = biotr(i,k,id_dia)
         if ( id_mes.ne.0) t_prog(ind_mes)%field(i,j,k,time%taum1) = biotr(i,k,id_mes)
         if ( id_poc.ne.0)   t_prog(ind_poc)%field(i,j,k,time%taum1) = biotr(i,k,id_poc)
         if ( id_caco3.ne.0) t_prog(ind_caco3)%field(i,j,k,time%taum1) = biotr(i,k,id_caco3)
         if ( id_fe.ne.0) t_prog(ind_fe)%field(i,j,k,time%taum1) = biotr(i,k,id_fe)
!         t_prog(ind_fe)%field(i,j,k,time%taum1) = biotr(i,k,id_fe) - dtts * tscav_fe(i,j) *        &
!                                                  max(0.0, (biotr(i,k,id_fe) - fe_bkgnd(i,j)) ) 
         if (id_dic.ne.0) &
         t_prog(ind_dic)%field(i,j,k,time%taum1) = t_prog(ind_dic)%field(i,j,k,time%taum1) +       &
                                                   122./16. * (nh4_bgc_change + no3_bgc_change)    &
                                                   - caco3_bgc_change
         if (id_adic.ne.0) &
         t_prog(ind_adic)%field(i,j,k,time%taum1) = t_prog(ind_adic)%field(i,j,k,time%taum1) +     &
                                                    122./16. * (nh4_bgc_change + no3_bgc_change)   &
                                                    - caco3_bgc_change 
         if (id_alk.ne.0) &
         t_prog(ind_alk)%field(i,j,k,time%taum1) = t_prog(ind_alk)%field(i,j,k,time%taum1) +       &  
                                                   ( -2.0 * caco3_bgc_change                       &
                                                   + nh4_bgc_change - no3_bgc_change)
         if (id_pchl.ne.0)  t_prog(ind_pchl)%field(i,j,k,time%taum1) = biotr(i,k,id_pchl) 
         if (id_dchl.ne.0)  t_prog(ind_dchl)%field(i,j,k,time%taum1) = biotr(i,k,id_dchl) 
         if (id_phyfe.ne.0) t_prog(ind_phyfe)%field(i,j,k,time%taum1) = biotr(i,k,id_phyfe) 
         if (id_diafe.ne.0) t_prog(ind_diafe)%field(i,j,k,time%taum1) = biotr(i,k,id_diafe) 
         if (id_zoofe.ne.0) t_prog(ind_zoofe)%field(i,j,k,time%taum1) = biotr(i,k,id_zoofe) 
         if (id_mesfe.ne.0) t_prog(ind_mesfe)%field(i,j,k,time%taum1) = biotr(i,k,id_mesfe) 
         if (id_detfe.ne.0) t_prog(ind_detfe)%field(i,j,k,time%taum1) = biotr(i,k,id_detfe) 
         if (id_pocfe.ne.0) t_prog(ind_pocfe)%field(i,j,k,time%taum1) = biotr(i,k,id_pocfe) 
         
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



       ! initialise advective sinking arrays
       adv_fb1(:,:) = 0.0
       adv_fb2(:,:) = 0.0
       adv_fb3(:,:) = 0.0
       adv_fb4(:,:) = 0.0

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! RASF upstream sinking of detritus and particulates !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do k=2,grid%nk+1
        do i=isc,iec
         adv_fb1(i,k) = wdetbio(i,j)*0.1 * biotr(i,k-1,id_det)
         adv_fb2(i,k) = wdetbio(i,j)*0.1 * biotr(i,k-1,id_detfe)
         adv_fb3(i,k) = wdetbio(i,j) * biotr(i,k-1,id_poc)
         adv_fb4(i,k) = wdetbio(i,j) * biotr(i,k-1,id_pocfe)
        enddo  ! i
       enddo  ! k

       ! Deposit tracer to sediment as tracer sinks through base of column, 
       !    accounting for burial and denitrification
       do i = isc,iec 
        k = grid%kmt(i,j)
        if (k.gt.0) then 
         biotic(n)%det_sed_depst(i,j) = adv_fb1(i,k+1) + adv_fb3(i,k+1)
         biotic(n)%detfe_sed_depst(i,j) = adv_fb2(i,k+1) + adv_fb4(i,k+1)
        endif  ! k.gt.0
       enddo  ! i

       do k =1,grid%nk  !{
        do i =isc,iec  !{
          t_prog(ind_det)%field(i,j,k,time%taum1) = t_prog(ind_det)%field(i,j,k,time%taum1) +      &
              grid%tmask(i,j,k) * dtts * (-adv_fb1(i,k+1) + adv_fb1(i,k))/Thickness%dzt(i,j,k)
          t_prog(ind_detfe)%field(i,j,k,time%taum1) = t_prog(ind_detfe)%field(i,j,k,time%taum1) +  &
              grid%tmask(i,j,k) * dtts * (-adv_fb2(i,k+1) + adv_fb2(i,k))/Thickness%dzt(i,j,k)
          t_prog(ind_poc)%field(i,j,k,time%taum1) = t_prog(ind_poc)%field(i,j,k,time%taum1) +      &
              grid%tmask(i,j,k) * dtts * (-adv_fb3(i,k+1) + adv_fb3(i,k))/Thickness%dzt(i,j,k)
          t_prog(ind_pocfe)%field(i,j,k,time%taum1) = t_prog(ind_pocfe)%field(i,j,k,time%taum1) +  &
              grid%tmask(i,j,k) * dtts * (-adv_fb4(i,k+1) + adv_fb4(i,k))/Thickness%dzt(i,j,k)
        enddo  !} i
       enddo  !} k


       adv_fb1(:,:) = 0.0
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! RASF upstream sinking of caco3 !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
           t_prog(ind_caco3)%field(i,j,k,time%taum1)=t_prog(ind_caco3)%field(i,j,k,time%taum1) +&
               grid%tmask(i,j,k) * dtts * &
               (-adv_fb1(i,k+1) + adv_fb1(i,k))/Thickness%dzt(i,j,k)
         enddo  ! i
        enddo  ! k

       endif  ! end loop for caco3

  enddo  !} j

enddo  !} n

return
end subroutine bio_v4
