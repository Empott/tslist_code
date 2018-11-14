!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine prints out the current value of variables at all specified
!   time series locations that are within the current patch.
!
! Michael G. Duda -- 25 August 2005
! vertical profiles added by Torge Lorenz -- June 2012
! Emily Potter -- November 2017. Changed to three letter file endings, momentum terms added and outputs only for domain 4. Altered lines marked by!EP or altered sections by !begin EP and !end EP.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calc_ts_locations( grid )

   USE module_domain, ONLY : domain, get_ijk_from_grid
   USE module_configure, ONLY : model_config_rec, grid_config_rec_type, model_to_grid_config_rec
   USE module_dm, ONLY : wrf_dm_min_real
   USE module_llxy
   USE module_state_description

   IMPLICIT NONE

   ! Arguments
   TYPE (domain), INTENT(INOUT) :: grid

   ! Externals
   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
   INTEGER, EXTERNAL :: get_unused_unit

   ! Local variables
   INTEGER :: ntsloc_temp
   INTEGER :: i, j, k, iunit
   REAL :: ts_rx, ts_ry, ts_xlat, ts_xlong, ts_hgt
   REAL :: known_lat, known_lon
   CHARACTER (LEN=132) :: message
   CHARACTER (LEN=24) :: ts_profile_filename
   CHARACTER (LEN=3), DIMENSION(39) :: ts_file_endings = (/ 'UUU', 'VVV', 'PHH', 'THH', 'QVV','WWW', &
                                                           'URU', 'UPH', 'UPG', 'UAD', 'UHD', 'UND', 'UCO', 'UCU', 'URS', & !EP
                                                           'VRV', 'VPH', 'VPG', 'VAD', 'VHD', 'VND', 'VCO', 'VCU', 'VRS', & !EP
                                                           'WRW', 'WBU', 'WAD', 'WHD', 'WND', 'WCO', 'WCU',               & !EP
                                                           'TRT', 'TBL', 'TRA', 'TAD', 'THD', 'TND', 'TDI', 'TRS'   /)   !EP
   TYPE (PROJ_INFO) :: ts_proj
   TYPE (grid_config_rec_type) :: config_flags


   INTEGER :: ids, ide, jds, jde, kds, kde,        &
              ims, ime, jms, jme, kms, kme,        &
              ips, ipe, jps, jpe, kps, kpe,        &
              imsx, imex, jmsx, jmex, kmsx, kmex,  &
              ipsx, ipex, jpsx, jpex, kpsx, kpex,  &
              imsy, imey, jmsy, jmey, kmsy, kmey,  &
              ipsy, ipey, jpsy, jpey, kpsy, kpey

   IF ( grid%id == 4) THEN !EP only output tslists for domain 4

   IF ( grid%ntsloc .LE. 0 ) RETURN

#if ((EM_CORE == 1) && (DA_CORE != 1))
   IF ( grid%dfi_stage == DFI_FST ) THEN
#endif
      CALL get_ijk_from_grid ( grid ,                               &
                               ids, ide, jds, jde, kds, kde,        &
                               ims, ime, jms, jme, kms, kme,        &
                               ips, ipe, jps, jpe, kps, kpe,        &
                               imsx, imex, jmsx, jmex, kmsx, kmex,  &
                               ipsx, ipex, jpsx, jpex, kpsx, kpex,  &
                               imsy, imey, jmsy, jmey, kmsy, kmey,  &
                               ipsy, ipey, jpsy, jpey, kpsy, kpey )

      CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags )
   
      ! Set up map transformation structure
      CALL map_init(ts_proj)
   
      IF (ips <= 1 .AND. 1 <= ipe .AND. &
          jps <= 1 .AND. 1 <= jpe) THEN
#if(NMM_CORE==1)
         known_lat = grid%hlat(1,1)
         known_lon = grid%hlon(1,1)
#else
         known_lat = grid%xlat(1,1)
         known_lon = grid%xlong(1,1)
#endif
      ELSE
         known_lat = 9999.
         known_lon = 9999.
      END IF
      known_lat = wrf_dm_min_real(known_lat)
      known_lon = wrf_dm_min_real(known_lon)
   
      ! Mercator
      IF (config_flags%map_proj == PROJ_MERC) THEN
         CALL map_set(PROJ_MERC, ts_proj,               &
                      truelat1 = config_flags%truelat1, &
                      lat1     = known_lat,             &
                      lon1     = known_lon,             &
                      knowni   = 1.,                    &
                      knownj   = 1.,                    &
                      dx       = config_flags%dx)
   
      ! Lambert conformal
      ELSE IF (config_flags%map_proj == PROJ_LC) THEN
         CALL map_set(PROJ_LC, ts_proj,                  &
                      truelat1 = config_flags%truelat1,  &
                      truelat2 = config_flags%truelat2,  &
                      stdlon   = config_flags%stand_lon, &
                      lat1     = known_lat,              &
                      lon1     = known_lon,              &
                      knowni   = 1.,                     &
                      knownj   = 1.,                     &
                      dx       = config_flags%dx)
   
      ! Polar stereographic
      ELSE IF (config_flags%map_proj == PROJ_PS) THEN
         CALL map_set(PROJ_PS, ts_proj,                  &
                      truelat1 = config_flags%truelat1,  &
                      stdlon   = config_flags%stand_lon, &
                      lat1     = known_lat,              &
                      lon1     = known_lon,              &
                      knowni   = 1.,                     &
                      knownj   = 1.,                     &
                      dx       = config_flags%dx)
   
#if (EM_CORE == 1)
      ! Cassini (global ARW)
      ELSE IF (config_flags%map_proj == PROJ_CASSINI) THEN
         CALL map_set(PROJ_CASSINI, ts_proj,                            &
                      latinc   = grid%dy*360.0/(2.0*EARTH_RADIUS_M*PI), &
                      loninc   = grid%dx*360.0/(2.0*EARTH_RADIUS_M*PI), & 
                      lat1     = known_lat,                             &
                      lon1     = known_lon,                             &
                      lat0     = config_flags%pole_lat,                 &
                      lon0     = config_flags%pole_lon,                 &
                      knowni   = 1.,                                    &
                      knownj   = 1.,                                    &
                      stdlon   = config_flags%stand_lon)
#endif

      ! Rotated latitude-longitude
      ELSE IF (config_flags%map_proj == PROJ_ROTLL) THEN
         CALL map_set(PROJ_ROTLL, ts_proj,                      &
! I have no idea how this should work for NMM nested domains
                      ixdim    = grid%e_we-1,                   &
                      jydim    = grid%e_sn-1,                   &
                      phi      = real(grid%e_sn-2)*grid%dy/2.0, &
                      lambda   = real(grid%e_we-2)*grid%dx,     &
                      lat1     = config_flags%cen_lat,          &
                      lon1     = config_flags%cen_lon,          &
                      latinc   = grid%dy,                       &
                      loninc   = grid%dx,                       &
                      stagger  = HH)
   
      END IF
   
      ! Determine time series locations for domain
      IF (.NOT. grid%have_calculated_tslocs) THEN
         grid%have_calculated_tslocs = .TRUE.
         WRITE(message, '(A43,I3)') 'Computing time series locations for domain ', grid%id
         CALL wrf_message(message)
   
         ntsloc_temp = 0
         DO k=1,grid%ntsloc
         
            IF (config_flags%map_proj == 0) THEN  ! For idealized cases, no map transformation needed
               ts_rx = grid%lattsloc(k)           ! NB: (x,y) = (lat,lon) rather than (x,y) = (lon,lat)
               ts_ry = grid%lontsloc(k)
            ELSE
               CALL latlon_to_ij(ts_proj, grid%lattsloc(k), grid%lontsloc(k), ts_rx, ts_ry)
            END IF
            

            ntsloc_temp = ntsloc_temp + 1
            grid%itsloc(ntsloc_temp) = NINT(ts_rx)
            grid%jtsloc(ntsloc_temp) = NINT(ts_ry)
            grid%id_tsloc(ntsloc_temp) = k
   
            ! Is point outside of domain (or on the edge of domain)?
            IF (grid%itsloc(ntsloc_temp) < ids .OR. grid%itsloc(ntsloc_temp) > ide .OR. &
                grid%jtsloc(ntsloc_temp) < jds .OR. grid%jtsloc(ntsloc_temp) > jde) THEN
               ntsloc_temp = ntsloc_temp - 1
   
            END IF
   
         END DO
   
         grid%next_ts_time = 1
   
         grid%ntsloc_domain = ntsloc_temp
   
         DO k=1,grid%ntsloc_domain
   
            ! If location is outside of patch, we need to get lat/lon of TS grid cell from another patch
            IF (grid%itsloc(k) < ips .OR. grid%itsloc(k) > ipe .OR. &
                grid%jtsloc(k) < jps .OR. grid%jtsloc(k) > jpe) THEN
               ts_xlat  = 1.E30
               ts_xlong = 1.E30
               ts_hgt   = 1.E30
            ELSE
#if(NMM_CORE==1)
               ts_xlat  = grid%hlat(grid%itsloc(k),grid%jtsloc(k))
               ts_xlong = grid%hlon(grid%itsloc(k),grid%jtsloc(k))
#else
               ts_xlat  = grid%xlat(grid%itsloc(k),grid%jtsloc(k))
               ts_xlong = grid%xlong(grid%itsloc(k),grid%jtsloc(k))
#endif
#if (EM_CORE == 1)
               ts_hgt   = grid%ht(grid%itsloc(k),grid%jtsloc(k))
#endif
            END IF
#if DM_PARALLEL
            ts_xlat  = wrf_dm_min_real(ts_xlat)
            ts_xlong = wrf_dm_min_real(ts_xlong)
            ts_hgt   = wrf_dm_min_real(ts_hgt)
#endif
   
            IF ( wrf_dm_on_monitor() ) THEN
               iunit = get_unused_unit()
               IF ( iunit <= 0 ) THEN
                  CALL wrf_error_fatal('Error in calc_ts_locations: could not find a free Fortran unit.')
               END IF
               WRITE(grid%ts_filename(k),'(A)') TRIM(grid%nametsloc(grid%id_tsloc(k)))//'.d00.TSS' !EP 3 letter ending
               i = LEN_TRIM(grid%ts_filename(k))
               WRITE(grid%ts_filename(k)(i-5:i-4),'(I2.2)') grid%id
               OPEN(UNIT=iunit, FILE=TRIM(grid%ts_filename(k)), FORM='FORMATTED', STATUS='REPLACE')
#if (EM_CORE == 1)
               WRITE(UNIT=iunit, &
                     FMT='(A26,I2,I3,A6,A2,F7.3,A1,F8.3,A3,I4,A1,I4,A3,F7.3,A1,F8.3,A2,F6.1,A7)') &
                     grid%desctsloc(grid%id_tsloc(k))//' ', grid%id, grid%id_tsloc(k), &
                     ' '//grid%nametsloc(grid%id_tsloc(k)), &
                     ' (', grid%lattsloc(grid%id_tsloc(k)), ',', grid%lontsloc(grid%id_tsloc(k)), ') (', &
                     grid%itsloc(k), ',', grid%jtsloc(k), ') (', &
                     ts_xlat, ',', ts_xlong, ') ', &
                     ts_hgt,' meters'
#else
               WRITE(UNIT=iunit, &
                     FMT='(A26,I2,I3,A6,A2,F7.3,A1,F8.3,A3,I4,A1,I4,A3,F7.3,A1,F8.3,A2)') &
                     grid%desctsloc(grid%id_tsloc(k))//' ', grid%id, grid%id_tsloc(k), &
                     ' '//grid%nametsloc(grid%id_tsloc(k)), &
                     ' (', grid%lattsloc(grid%id_tsloc(k)), ',', grid%lontsloc(grid%id_tsloc(k)), ') (', &
                     grid%itsloc(k), ',', grid%jtsloc(k), ') (', &
                     ts_xlat, ',', ts_xlong, ') '
#endif
               CLOSE(UNIT=iunit)

               ts_profile_filename = grid%ts_filename(k)
               DO j=1,39  
                  ! Create the output files for the vertical profiles, one file for each variable
                  iunit = get_unused_unit()
                  IF ( iunit <= 0 ) THEN
                     CALL wrf_error_fatal('Error in calc_ts_locations: could not find a free Fortran unit.')
                  END IF
                  i = LEN_TRIM(ts_profile_filename)
                  WRITE(ts_profile_filename(i-2:i),'(A3)') ts_file_endings(j)  !EP changed to 3 letter ending
                  OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), FORM='FORMATTED', STATUS='REPLACE')
#if (EM_CORE == 1)
                  WRITE(UNIT=iunit, &
                        FMT='(A26,I2,I3,A6,A2,F7.3,A1,F8.3,A3,I4,A1,I4,A3,F7.3,A1,F8.3,A2,F6.1,A7)') &
                        grid%desctsloc(grid%id_tsloc(k))//' ', grid%id, grid%id_tsloc(k), &
                        ' '//grid%nametsloc(grid%id_tsloc(k)), &
                        ' (', grid%lattsloc(grid%id_tsloc(k)), ',', grid%lontsloc(grid%id_tsloc(k)), ') (', &
                        grid%itsloc(k), ',', grid%jtsloc(k), ') (', &
                        ts_xlat, ',', ts_xlong, ') ', &
                        ts_hgt,' meters'
#else
                  WRITE(UNIT=iunit, &
                        FMT='(A26,I2,I3,A6,A2,F7.3,A1,F8.3,A3,I4,A1,I4,A3,F7.3,A1,F8.3,A2)') &
                        grid%desctsloc(grid%id_tsloc(k))//' ', grid%id, grid%id_tsloc(k), &
                        ' '//grid%nametsloc(grid%id_tsloc(k)), &
                        ' (', grid%lattsloc(grid%id_tsloc(k)), ',', grid%lontsloc(grid%id_tsloc(k)), ') (', &
                        grid%itsloc(k), ',', grid%jtsloc(k), ') (', &
                        ts_xlat, ',', ts_xlong, ') '
#endif               
                  CLOSE(UNIT=iunit)
               END DO
            END IF
         END DO
   
      END IF
#if ((EM_CORE == 1) && (DA_CORE != 1))
   END IF
#endif

   END IF !EP grid%id ==4

END SUBROUTINE calc_ts_locations


SUBROUTINE calc_ts( grid )

   USE module_domain
   USE module_model_constants

   IMPLICIT NONE

   ! Arguments
   TYPE (domain), INTENT(INOUT) :: grid

   LOGICAL, EXTERNAL :: wrf_dm_on_monitor

   ! Local variables
   INTEGER :: i, k, mm, n, ix, iy, rc
   REAL :: earth_u, earth_v,                       &
           output_t, output_q, clw, xtime_minutes, output_w, &
           output_trt, output_tbl, output_tra, output_tad, output_thd, output_tnd, output_trs, output_tdi,   &
           output_wrw, output_wnd, output_wad, output_wbu, output_wco, output_wcu, output_whd,                &
           earth_uru, earth_uph, earth_upg, earth_uad, earth_uhd, earth_und, earth_uco, earth_ucu, earth_urs, & !EP
           earth_vrv, earth_vph, earth_vpg, earth_vad, earth_vhd, earth_vnd, earth_vco, earth_vcu, earth_vrs, & !EP
           earth_wrw, earth_wbu, earth_wad, earth_whd, earth_wnd, earth_wco, earth_wcu,  & !EP
           earth_trt, earth_tbl, earth_tra, earth_tad, earth_thd, earth_tnd, earth_tdi, earth_trs !EP

   REAL, ALLOCATABLE, DIMENSION(:) :: p8w
   REAL, ALLOCATABLE, DIMENSION(:) :: earth_u_profile, earth_v_profile, &
                                      earth_uru_profile, earth_uph_profile, earth_upg_profile, earth_uad_profile, & !EP
                                      earth_uhd_profile, earth_und_profile, earth_uco_profile, earth_ucu_profile, earth_urs_profile, & !EP
                                      earth_vrv_profile, earth_vph_profile, earth_vpg_profile, earth_vad_profile, & !EP
                                      earth_vhd_profile, earth_vnd_profile, earth_vco_profile, earth_vcu_profile, earth_vrs_profile !EP


   ! Parameter ts_model_level:  
       ! TRUE to output T, Q, and wind at lowest model level
       ! FALSE to output T and Q at 2-m and wind at 10-m diagnostic levels:
   LOGICAL, PARAMETER :: ts_model_level = .FALSE. 
 
   IF ( grid%id == 4 )  THEN !EP only if d04

   !Allocate the arrays for wind components 
#if ( EM_CORE == 1 )   
   ALLOCATE ( earth_u_profile(grid%max_ts_level), earth_v_profile(grid%max_ts_level) )
!begin EP
   ALLOCATE ( earth_uru_profile(grid%max_ts_level), earth_uph_profile(grid%max_ts_level), earth_upg_profile(grid%max_ts_level) )
   ALLOCATE ( earth_uad_profile(grid%max_ts_level), earth_uhd_profile(grid%max_ts_level), earth_und_profile(grid%max_ts_level) )
   ALLOCATE ( earth_uco_profile(grid%max_ts_level), earth_ucu_profile(grid%max_ts_level), earth_urs_profile(grid%max_ts_level) )
   ALLOCATE ( earth_vrv_profile(grid%max_ts_level), earth_vph_profile(grid%max_ts_level), earth_vpg_profile(grid%max_ts_level) )
   ALLOCATE ( earth_vad_profile(grid%max_ts_level), earth_vhd_profile(grid%max_ts_level), earth_vnd_profile(grid%max_ts_level) )
   ALLOCATE ( earth_vco_profile(grid%max_ts_level), earth_vcu_profile(grid%max_ts_level), earth_vrs_profile(grid%max_ts_level) )
#endif
!end EP
   IF ( grid%ntsloc_domain .LE. 0 ) RETURN

#if ((EM_CORE == 1) && (DA_CORE != 1))
   IF ( grid%dfi_opt /= DFI_NODFI .AND. grid%dfi_stage /= DFI_FST ) RETURN
#endif

   n = grid%next_ts_time

   ALLOCATE(p8w(grid%sm32:grid%em32))

   DO i=1,grid%ntsloc_domain

      ix = grid%itsloc(i)
      iy = grid%jtsloc(i)
  
      IF (grid%sp31 <= ix .AND. ix <= grid%ep31 .AND. &
          grid%sp33 <= iy .AND. iy <= grid%ep33) THEN
       
         IF (ts_model_level) THEN
   
            !
            ! Output from the lowest model computational level: EP additional comment:cosa is the local cosine of map rotation.
#if (EM_CORE == 1)
            earth_u = grid%u_2(ix,1,iy)*grid%cosa(ix,iy)-grid%v_2(ix,1,iy)*grid%sina(ix,iy)
            earth_v = grid%v_2(ix,1,iy)*grid%cosa(ix,iy)+grid%u_2(ix,1,iy)*grid%sina(ix,iy)
             !begin EP
            earth_uru = grid%ru_tend(ix,1,iy)*grid%cosa(ix,iy)-grid%rv_tend(ix,1,iy)*grid%sina(ix,iy)
            earth_uph = grid%ru_tend_phys(ix,1,iy)*grid%cosa(ix,iy)-grid%rv_tend_phys(ix,1,iy)*grid%sina(ix,iy)
            earth_upg = grid%ru_tend_pgf(ix,1,iy)*grid%cosa(ix,iy)-grid%rv_tend_pgf(ix,1,iy)*grid%sina(ix,iy)
            earth_uad = grid%ru_tend_adv(ix,1,iy)*grid%cosa(ix,iy)-grid%rv_tend_adv(ix,1,iy)*grid%sina(ix,iy)
            earth_uhd = grid%ru_tend_hdiff(ix,1,iy)*grid%cosa(ix,iy)-grid%rv_tend_hdiff(ix,1,iy)*grid%sina(ix,iy)
            earth_und = grid%ru_tend_6thdiff(ix,1,iy)*grid%cosa(ix,iy)-grid%rv_tend_6thdiff(ix,1,iy)*grid%sina(ix,iy)
            earth_uco = grid%ru_tend_cor(ix,1,iy)*grid%cosa(ix,iy)-grid%rv_tend_cor(ix,1,iy)*grid%sina(ix,iy)
            earth_ucu = grid%ru_tend_curv(ix,1,iy)*grid%cosa(ix,iy)-grid%rv_tend_curv(ix,1,iy)*grid%sina(ix,iy)
            earth_urs = grid%utend_sound(ix,1,iy)*grid%cosa(ix,iy)-grid%vtend_sound(ix,1,iy)*grid%sina(ix,iy)

            earth_vrv = grid%rv_tend(ix,1,iy)*grid%cosa(ix,iy)+grid%ru_tend(ix,1,iy)*grid%sina(ix,iy)
            earth_vph =grid%rv_tend_phys(ix,1,iy)*grid%cosa(ix,iy)+grid%ru_tend_phys(ix,1,iy)*grid%sina(ix,iy)
            earth_vpg =grid%rv_tend_pgf(ix,1,iy)*grid%cosa(ix,iy)+grid%ru_tend_pgf(ix,1,iy)*grid%sina(ix,iy)
            earth_vad =grid%rv_tend_adv(ix,1,iy)*grid%cosa(ix,iy)+grid%ru_tend_adv(ix,1,iy)*grid%sina(ix,iy)
            earth_vhd =grid%rv_tend_hdiff(ix,1,iy)*grid%cosa(ix,iy)+grid%ru_tend_hdiff(ix,1,iy)*grid%sina(ix,iy)
            earth_vnd =grid%rv_tend_6thdiff(ix,1,iy)*grid%cosa(ix,iy)+grid%ru_tend_6thdiff(ix,1,iy)*grid%sina(ix,iy)
            earth_vco =grid%rv_tend_cor(ix,1,iy)*grid%cosa(ix,iy)+grid%ru_tend_cor(ix,1,iy)*grid%sina(ix,iy)
            earth_vcu =grid%rv_tend_curv(ix,1,iy)*grid%cosa(ix,iy)+grid%ru_tend_curv(ix,1,iy)*grid%sina(ix,iy)
            earth_vrs =grid%vtend_sound(ix,1,iy)*grid%cosa(ix,iy)+grid%utend_sound(ix,1,iy)*grid%sina(ix,iy)
            !end EP
            output_t = grid%t_2(ix,1,iy)

            output_w = grid%w_2(ix,1,iy)
            
#else
            earth_u = grid%u(ix,1,iy)
            earth_v = grid%v(ix,1,iy)
            !begin EP
            earth_uru = grid%ru_tend(ix,1,iy)
            earth_uph = grid%ru_tend_phys(ix,1,iy)
            earth_upg = grid%ru_tend_pgf(ix,1,iy)
            earth_uad = grid%ru_tend_adv(ix,1,iy)
            earth_uhd = grid%ru_tend_hdiff(ix,1,iy)
            earth_und = grid%ru_tend_6thdiff(ix,1,iy)
            earth_uco = grid%ru_tend_cor(ix,1,iy)
            earth_ucu = grid%ru_tend_curv(ix,1,iy)
            earth_urs = grid%utend_sound(ix,1,iy)

            earth_vrv = grid%rv_tend(ix,1,iy)
            earth_vph =grid%rv_tend_phys(ix,1,iy)
            earth_vpg =grid%rv_tend_pgf(ix,1,iy)
            earth_vad =grid%rv_tend_adv(ix,1,iy)
            earth_vhd =grid%rv_tend_hdiff(ix,1,iy)
            earth_vnd =grid%rv_tend_6thdiff(ix,1,iy)
            earth_vco =grid%rv_tend_cor(ix,1,iy)
            earth_vcu =grid%rv_tend_curv(ix,1,iy)
            earth_vrs =grid%vtend_sound(ix,1,iy)
            !end EP
            output_t = grid%t(ix,1,iy)
            output_w = grid%w(ix,1,iy)

#endif
            !begin EP
            output_trt = grid%t_tend(ix,1,iy)
            output_tbl = grid%t_tend_phys_bl(ix,1,iy)
            output_tra = grid%t_tend_phys_rad(ix,1,iy)
            output_tad = grid%t_tend_adv(ix,1,iy)
            output_thd = grid%t_tend_hdiff(ix,1,iy)
            output_tnd = grid%t_tend_6thdiff(ix,1,iy)
            output_tdi = grid%h_diabatic(ix,1,iy)
            output_trs = grid%ttend_sound(ix,1,iy)

            output_wrw = grid%rw_tend(ix,1,iy)
            output_wbu = grid%rw_tend_buoy(ix,1,iy)
            output_wad = grid%rw_tend_adv(ix,1,iy)
            output_whd = grid%rw_tend_hdiff(ix,1,iy)
            output_wnd = grid%rw_tend_6thdiff(ix,1,iy)
            output_wco = grid%rw_tend_cor(ix,1,iy)
            output_wcu = grid%rw_tend_curv(ix,1,iy)
            !end EP
            output_q = grid%moist(ix,1,iy,P_QV)
   
         ELSE
   
            !
            ! Output at 2-m and 10-m diagnostic levels:
#if (EM_CORE == 1)
            DO k=1,grid%max_ts_level
            earth_u_profile(k) = grid%u_2(ix,k,iy)*grid%cosa(ix,iy)-grid%v_2(ix,k,iy)*grid%sina(ix,iy)
            earth_v_profile(k) = grid%v_2(ix,k,iy)*grid%cosa(ix,iy)+grid%u_2(ix,k,iy)*grid%sina(ix,iy)
            !begin EP
            earth_uru_profile(k) = grid%ru_tend(ix,k,iy)*grid%cosa(ix,iy)-grid%rv_tend(ix,k,iy)*grid%sina(ix,iy)
            earth_vrv_profile(k) = grid%rv_tend(ix,k,iy)*grid%cosa(ix,iy)+grid%ru_tend(ix,k,iy)*grid%sina(ix,iy)

            earth_uru_profile(k) = grid%ru_tend(ix,k,iy)*grid%cosa(ix,iy)-grid%rv_tend(ix,k,iy)*grid%sina(ix,iy)
            earth_uph_profile(k) = grid%ru_tend_phys(ix,k,iy)*grid%cosa(ix,iy)-grid%rv_tend_phys(ix,k,iy)*grid%sina(ix,iy)
            earth_upg_profile(k) = grid%ru_tend_pgf(ix,k,iy)*grid%cosa(ix,iy)-grid%rv_tend_pgf(ix,k,iy)*grid%sina(ix,iy)
            earth_uad_profile(k) = grid%ru_tend_adv(ix,k,iy)*grid%cosa(ix,iy)-grid%rv_tend_adv(ix,k,iy)*grid%sina(ix,iy)
            earth_uhd_profile(k) = grid%ru_tend_hdiff(ix,k,iy)*grid%cosa(ix,iy)-grid%rv_tend_hdiff(ix,k,iy)*grid%sina(ix,iy)
            earth_und_profile(k) = grid%ru_tend_6thdiff(ix,k,iy)*grid%cosa(ix,iy)-grid%rv_tend_6thdiff(ix,k,iy)*grid%sina(ix,iy)
            earth_uco_profile(k) = grid%ru_tend_cor(ix,k,iy)*grid%cosa(ix,iy)-grid%rv_tend_cor(ix,k,iy)*grid%sina(ix,iy)
            earth_ucu_profile(k) = grid%ru_tend_curv(ix,k,iy)*grid%cosa(ix,iy)-grid%rv_tend_curv(ix,k,iy)*grid%sina(ix,iy)
            earth_urs_profile(k) = grid%utend_sound(ix,k,iy)*grid%cosa(ix,iy)-grid%vtend_sound(ix,k,iy)*grid%sina(ix,iy)

            earth_vrv_profile(k) = grid%rv_tend(ix,k,iy)*grid%cosa(ix,iy)+grid%ru_tend(ix,k,iy)*grid%sina(ix,iy)
            earth_vph_profile(k) =grid%rv_tend_phys(ix,k,iy)*grid%cosa(ix,iy)+grid%ru_tend_phys(ix,k,iy)*grid%sina(ix,iy)
            earth_vpg_profile(k) =grid%rv_tend_pgf(ix,k,iy)*grid%cosa(ix,iy)+grid%ru_tend_pgf(ix,k,iy)*grid%sina(ix,iy)
            earth_vad_profile(k) =grid%rv_tend_adv(ix,k,iy)*grid%cosa(ix,iy)+grid%ru_tend_adv(ix,k,iy)*grid%sina(ix,iy)
            earth_vhd_profile(k) =grid%rv_tend_hdiff(ix,k,iy)*grid%cosa(ix,iy)+grid%ru_tend_hdiff(ix,k,iy)*grid%sina(ix,iy)
            earth_vnd_profile(k) =grid%rv_tend_6thdiff(ix,k,iy)*grid%cosa(ix,iy)+grid%ru_tend_6thdiff(ix,k,iy)*grid%sina(ix,iy)
            earth_vco_profile(k) =grid%rv_tend_cor(ix,k,iy)*grid%cosa(ix,iy)+grid%ru_tend_cor(ix,k,iy)*grid%sina(ix,iy)
            earth_vcu_profile(k) =grid%rv_tend_curv(ix,k,iy)*grid%cosa(ix,iy)+grid%ru_tend_curv(ix,k,iy)*grid%sina(ix,iy)
            earth_vrs_profile(k) =grid%vtend_sound(ix,k,iy)*grid%cosa(ix,iy)+grid%utend_sound(ix,k,iy)*grid%sina(ix,iy)
            !end EP
            END DO
            earth_u = grid%u10(ix,iy)*grid%cosa(ix,iy)-grid%v10(ix,iy)*grid%sina(ix,iy)
            earth_v = grid%v10(ix,iy)*grid%cosa(ix,iy)+grid%u10(ix,iy)*grid%sina(ix,iy)
            output_q = grid%q2(ix,iy)
#else
            earth_u = grid%u10(ix,iy)
            earth_v = grid%v10(ix,iy)
            output_q = grid%qsfc(ix,iy)
#endif
            output_t = grid%t2(ix,iy)
   
         END IF
   
#if (EM_CORE == 1)
         ! Calculate column-integrated liquid/ice  (kg/m^2 or mm)
         CALL calc_p8w(grid, ix, iy, p8w, grid%sm32, grid%em32)
         clw=0.
         DO mm = 1, num_moist
            IF ( (mm == P_QC) .OR. (mm == P_QR) .OR. (mm == P_QI) .OR. &
                 (mm == P_QS) .OR. (mm == P_QG) ) THEN
               DO k=grid%sm32,grid%em32-1
                  clw=clw+grid%moist(ix,k,iy,mm)*(p8w(k)-p8w(k+1))
               END DO
           END IF
         END DO
         clw = clw / g
#endif
   
         CALL domain_clock_get( grid, minutesSinceSimulationStart=xtime_minutes )
         grid%ts_hour(n,i) = xtime_minutes / 60.
#if (EM_CORE == 1)
            !Create vertical profiles, from lowest model level up to desired level max_ts_level
            DO k=1,grid%max_ts_level
            grid%ts_u_profile(n,i,k)   = earth_u_profile(k)
            grid%ts_v_profile(n,i,k)   = earth_v_profile(k)
            !begin EP
            grid%ts_uru_profile(n,i,k) = earth_uru_profile(k)
            grid%ts_uph_profile(n,i,k) = earth_uph_profile(k)
            grid%ts_upg_profile(n,i,k) = earth_upg_profile(k)
            grid%ts_uad_profile(n,i,k) = earth_uad_profile(k) 
            grid%ts_uhd_profile(n,i,k) = earth_uhd_profile(k) 
            grid%ts_und_profile(n,i,k) = earth_und_profile(k) 
            grid%ts_uco_profile(n,i,k) = earth_uco_profile(k) 
            grid%ts_ucu_profile(n,i,k) = earth_ucu_profile(k) 
            grid%ts_urs_profile(n,i,k) = earth_urs_profile(k) 

            grid%ts_vrv_profile(n,i,k) = earth_vrv_profile(k)
            grid%ts_vph_profile(n,i,k) = earth_vph_profile(k) 
            grid%ts_vpg_profile(n,i,k) = earth_vpg_profile(k) 
            grid%ts_vad_profile(n,i,k) = earth_vad_profile(k) 
            grid%ts_vhd_profile(n,i,k) = earth_vhd_profile(k) 
            grid%ts_vnd_profile(n,i,k) = earth_vnd_profile(k) 
            grid%ts_vco_profile(n,i,k) = earth_vco_profile(k) 
            grid%ts_vcu_profile(n,i,k) = earth_vcu_profile(k) 
            grid%ts_vrs_profile(n,i,k) = earth_vrs_profile(k)
            !end EP
            grid%ts_gph_profile(n,i,k) = (grid%phb(ix,k,iy)+grid%ph_2(ix,k,iy))/9.81 
            grid%ts_th_profile(n,i,k)  = grid%t_2(ix,k,iy) + 300
            grid%ts_w_profile(n,i,k)   = grid%w_2(ix,k,iy)
            !begin EP
            grid%ts_trt_profile(n,i,k) = grid%t_tend(ix,k,iy)
            grid%ts_tbl_profile(n,i,k) = grid%t_tend_phys_bl(ix,k,iy)
            grid%ts_tra_profile(n,i,k) = grid%t_tend_phys_rad(ix,k,iy)
            grid%ts_tad_profile(n,i,k) = grid%t_tend_adv(ix,k,iy)
            grid%ts_thd_profile(n,i,k) = grid%t_tend_hdiff(ix,k,iy)
            grid%ts_tnd_profile(n,i,k) = grid%t_tend_6thdiff(ix,k,iy)
            grid%ts_tdi_profile(n,i,k) = grid%h_diabatic(ix,k,iy)
            grid%ts_trs_profile(n,i,k) = grid%ttend_sound(ix,k,iy)

            grid%ts_wrw_profile(n,i,k) = grid%rw_tend(ix,k,iy)
            grid%ts_wbu_profile(n,i,k) = grid%rw_tend_buoy(ix,k,iy)
            grid%ts_wad_profile(n,i,k) = grid%rw_tend_adv(ix,k,iy)
            grid%ts_whd_profile(n,i,k) = grid%rw_tend_hdiff(ix,k,iy)
            grid%ts_wnd_profile(n,i,k) = grid%rw_tend_6thdiff(ix,k,iy)
            grid%ts_wco_profile(n,i,k) = grid%rw_tend_cor(ix,k,iy)
            grid%ts_wcu_profile(n,i,k) = grid%rw_tend_curv(ix,k,iy)
            !end EP

            grid%ts_qv_profile(n,i,k)  = grid%moist(ix,k,iy,P_QV)
            END DO
#endif
         grid%ts_u(n,i)    = earth_u
         grid%ts_v(n,i)    = earth_v
         !begin EP
         grid%ts_uru(n,i) = earth_uru 
         grid%ts_uph(n,i) = earth_uph 
         grid%ts_upg(n,i) = earth_upg 
         grid%ts_uad(n,i) = earth_uad 
         grid%ts_uhd(n,i) = earth_uhd 
         grid%ts_und(n,i) = earth_und 
         grid%ts_uco(n,i) = earth_uco 
         grid%ts_ucu(n,i) = earth_ucu 
         grid%ts_urs(n,i) = earth_urs 

         grid%ts_vrv(n,i) = earth_vrv 
         grid%ts_vph(n,i) = earth_vph 
         grid%ts_vpg(n,i) = earth_vpg 
         grid%ts_vad(n,i) = earth_vad 
         grid%ts_vhd(n,i) = earth_vhd 
         grid%ts_vnd(n,i) = earth_vnd 
         grid%ts_vco(n,i) = earth_vco 
         grid%ts_vcu(n,i) = earth_vcu 
         grid%ts_vrs(n,i) = earth_vrs 
         !end EP

         grid%ts_t(n,i)    = output_t
         grid%ts_w(n,i)    = output_w
         !begin EP
         grid%ts_trt(n,i) = output_trt
         grid%ts_tbl(n,i) = output_tbl
         grid%ts_tra(n,i) = output_tra
         grid%ts_tad(n,i) = output_tad
         grid%ts_thd(n,i) = output_thd
         grid%ts_tnd(n,i) = output_tnd
         grid%ts_tdi(n,i) = output_tdi
         grid%ts_trs(n,i) = output_trs

         grid%ts_wrw(n,i) = output_wrw
         grid%ts_wbu(n,i) = output_wbu
         grid%ts_wad(n,i) = output_wad
         grid%ts_whd(n,i) = output_whd
         grid%ts_wnd(n,i) = output_wnd
         grid%ts_wco(n,i) = output_wco
         grid%ts_wcu(n,i) = output_wcu
         !end EP

         grid%ts_q(n,i)    = output_q
         grid%ts_psfc(n,i) = grid%psfc(ix,iy)
#if (EM_CORE == 1)
!EP mut added
         grid%ts_mut(n,i)  = grid%mut(ix,iy)
         grid%ts_muu(n,i)  = grid%muu(ix,iy)
         grid%ts_muv(n,i)  = grid%muv(ix,iy)
         grid%ts_glw(n,i)  = grid%glw(ix,iy)
         grid%ts_gsw(n,i)  = grid%gsw(ix,iy)
         grid%ts_hfx(n,i)  = grid%hfx(ix,iy)
         grid%ts_lh(n,i)   = grid%lh(ix,iy)
         grid%ts_clw(n,i)  = clw
         grid%ts_rainc(n,i)  = grid%rainc(ix,iy)
         grid%ts_rainnc(n,i) = grid%rainnc(ix,iy)
         grid%ts_tsk(n,i)  = grid%tsk(ix,iy)
#else
         grid%ts_tsk(n,i)  = grid%nmm_tsk(ix,iy)
#endif
         grid%ts_tslb(n,i) = grid%tslb(ix,1,iy)
   
      ELSE
#if (EM_CORE == 1 )   
         DO k=1,grid%max_ts_level
         grid%ts_u_profile(n,i,k)     = 1.E30
         grid%ts_v_profile(n,i,k)     = 1.E30
         grid%ts_gph_profile(n,i,k)   = 1.E30
         grid%ts_th_profile(n,i,k)    = 1.E30
         grid%ts_qv_profile(n,i,k)    = 1.E30
         !begin EP
         grid%ts_uru_profile(n,i,k)   = 1.E30
         grid%ts_uph_profile(n,i,k)   = 1.E30
         grid%ts_upg_profile(n,i,k)   = 1.E30
         grid%ts_uad_profile(n,i,k)   = 1.E30 
         grid%ts_uhd_profile(n,i,k)   = 1.E30 
         grid%ts_und_profile(n,i,k)   = 1.E30 
         grid%ts_uco_profile(n,i,k)   = 1.E30 
         grid%ts_ucu_profile(n,i,k)   = 1.E30 
         grid%ts_urs_profile(n,i,k)   = 1.E30 

         grid%ts_vrv_profile(n,i,k)   = 1.E30
         grid%ts_vph_profile(n,i,k)   = 1.E30 
         grid%ts_vpg_profile(n,i,k)   = 1.E30 
         grid%ts_vad_profile(n,i,k)   = 1.E30 
         grid%ts_vhd_profile(n,i,k)   = 1.E30 
         grid%ts_vnd_profile(n,i,k)   = 1.E30 
         grid%ts_vco_profile(n,i,k)   = 1.E30 
         grid%ts_vcu_profile(n,i,k)   = 1.E30 
         grid%ts_vrs_profile(n,i,k)   = 1.E30  

         grid%ts_trt_profile(n,i,k)   = 1.E30
         grid%ts_tbl_profile(n,i,k)   = 1.E30
         grid%ts_tra_profile(n,i,k)   = 1.E30
         grid%ts_tad_profile(n,i,k)   = 1.E30
         grid%ts_thd_profile(n,i,k)   = 1.E30
         grid%ts_tnd_profile(n,i,k)   = 1.E30
         grid%ts_tdi_profile(n,i,k)   = 1.E30
         grid%ts_trs_profile(n,i,k)   = 1.E30

         grid%ts_wrw_profile(n,i,k)   = 1.E30
         grid%ts_wbu_profile(n,i,k)   = 1.E30
         grid%ts_wad_profile(n,i,k)   = 1.E30
         grid%ts_whd_profile(n,i,k)   = 1.E30
         grid%ts_wnd_profile(n,i,k)   = 1.E30
         grid%ts_wco_profile(n,i,k)   = 1.E30
         grid%ts_wcu_profile(n,i,k)   = 1.E30
         !end EP
         END DO
#endif   
         grid%ts_hour(n,i) = 1.E30
         grid%ts_u(n,i)    = 1.E30
         grid%ts_v(n,i)    = 1.E30
         grid%ts_t(n,i)    = 1.E30
         grid%ts_q(n,i)    = 1.E30
         grid%ts_psfc(n,i) = 1.E30
         !begin EP
         grid%ts_uru(n,i)   = 1.E30
         grid%ts_uph(n,i)   = 1.E30
         grid%ts_upg(n,i)   = 1.E30
         grid%ts_uad(n,i)   = 1.E30
         grid%ts_uhd(n,i)   = 1.E30
         grid%ts_und(n,i)   = 1.E30
         grid%ts_uco(n,i)   = 1.E30
         grid%ts_ucu(n,i)   = 1.E30
         grid%ts_urs(n,i)   = 1.E30

         grid%ts_vrv(n,i)   = 1.E30
         grid%ts_vph(n,i)   = 1.E30
         grid%ts_vpg(n,i)   = 1.E30
         grid%ts_vad(n,i)   = 1.E30
         grid%ts_vhd(n,i)   = 1.E30
         grid%ts_vnd(n,i)   = 1.E30
         grid%ts_vco(n,i)   = 1.E30
         grid%ts_vcu(n,i)   = 1.E30
         grid%ts_vrs(n,i)   = 1.E30

         grid%ts_trt(n,i)   = 1.E30
         grid%ts_tbl(n,i)   = 1.E30
         grid%ts_tra(n,i)   = 1.E30
         grid%ts_tad(n,i)   = 1.E30
         grid%ts_thd(n,i)   = 1.E30
         grid%ts_tnd(n,i)   = 1.E30
         grid%ts_tdi(n,i)   = 1.E30
         grid%ts_trs(n,i)   = 1.E30

         grid%ts_wrw(n,i)   = 1.E30
         grid%ts_wbu(n,i)   = 1.E30
         grid%ts_wad(n,i)   = 1.E30
         grid%ts_whd(n,i)   = 1.E30
         grid%ts_wnd(n,i)   = 1.E30
         grid%ts_wco(n,i)   = 1.E30
         grid%ts_wcu(n,i)   = 1.E30
         !end EP



#if (EM_CORE == 1)
         grid%ts_mut(n,i)  = 1.E30  !EP
         grid%ts_muu(n,i)  = 1.E30  !EP
         grid%ts_muv(n,i)  = 1.E30  !EP
         grid%ts_glw(n,i)  = 1.E30
         grid%ts_gsw(n,i)  = 1.E30
         grid%ts_hfx(n,i)  = 1.E30
         grid%ts_lh(n,i)   = 1.E30
         grid%ts_clw(n,i)  = 1.E30
         grid%ts_rainc(n,i)  = 1.E30
         grid%ts_rainnc(n,i) = 1.E30
#endif
         grid%ts_tsk(n,i)  = 1.E30
         grid%ts_tslb(n,i) = 1.E30
   
      END IF
   END DO

#if (EM_CORE == 1) 
   DEALLOCATE(p8w, earth_u_profile, earth_v_profile)
   !begin EP
   DEALLOCATE(earth_uru_profile, earth_uph_profile, earth_upg_profile, earth_uad_profile)
   DEALLOCATE(earth_uhd_profile, earth_und_profile, earth_uco_profile, earth_ucu_profile)
   DEALLOCATE(earth_urs_profile, earth_vrv_profile, earth_vph_profile, earth_vpg_profile)
   DEALLOCATE(earth_vad_profile, earth_vhd_profile, earth_vnd_profile, earth_vco_profile)
   DEALLOCATE(earth_vcu_profile, earth_vrs_profile)
   !end EP
#endif 
 
   grid%next_ts_time = grid%next_ts_time + 1  

   IF ( grid%next_ts_time > grid%ts_buf_size ) CALL write_ts(grid)
   
   END IF !EP grid%id ==4

END SUBROUTINE calc_ts


SUBROUTINE write_ts( grid )

   USE module_domain, ONLY : domain
   USE module_dm, ONLY : wrf_dm_min_reals
   USE module_state_description

   IMPLICIT NONE

   ! Arguments
   TYPE (domain), INTENT(INOUT) :: grid

   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
   INTEGER, EXTERNAL :: get_unused_unit

   ! Local variables
   INTEGER :: i, n, ix, iy, iunit, k
   REAL, ALLOCATABLE, DIMENSION(:,:) :: ts_buf
   CHARACTER (LEN=24) :: ts_profile_filename
   CHARACTER (LEN=26) :: profile_format

   IF ( grid%ntsloc_domain .LE. 0 ) RETURN

#if ((EM_CORE == 1) && (DA_CORE != 1))
   IF ( grid%dfi_opt /= DFI_NODFI .AND. grid%dfi_stage /= DFI_FST ) RETURN
#endif

#ifdef DM_PARALLEL
   ALLOCATE(ts_buf(grid%ts_buf_size,grid%max_ts_locs))

   ts_buf(:,:) = grid%ts_hour(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_hour(:,:),grid%ts_buf_size*grid%max_ts_locs)

#if (EM_CORE == 1) 
   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_u_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_u_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO
 
   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_v_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_v_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO


!begin EP

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_uru_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_uru_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_uph_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_uph_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_upg_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_upg_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_uad_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_uad_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_uhd_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_uhd_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_und_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_und_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_uco_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_uco_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_ucu_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_ucu_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO


   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_urs_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_urs_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO


   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_vrv_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vrv_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_vph_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vph_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_vpg_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vpg_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_vad_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vad_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_vhd_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vhd_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_vnd_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vnd_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_vco_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vco_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_vcu_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vcu_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO


   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_vrs_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vrs_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

!end EP


   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_gph_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_gph_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_th_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_th_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

!begin EP

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_trt_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_trt_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_tbl_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_tbl_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_tra_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_tra_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_tad_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_tad_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_thd_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_thd_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_tnd_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_tnd_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_trs_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_trs_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_tdi_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_tdi_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

!end EP
   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_w_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_w_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO
!begin EP
   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_wrw_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_wrw_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_wbu_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_wbu_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_wad_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_wad_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_whd_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_whd_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_wnd_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_wnd_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_wco_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_wco_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO

   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_wcu_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_wcu_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO
!end EP


   DO k=1,grid%max_ts_level
   ts_buf(:,:) = grid%ts_qv_profile(:,:,k)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_qv_profile(:,:,k),grid%ts_buf_size*grid%max_ts_locs)
   END DO
#endif         

   ts_buf(:,:) = grid%ts_u(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_u(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_v(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_v(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_t(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_t(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_q(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_q(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_psfc(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_psfc(:,:),grid%ts_buf_size*grid%max_ts_locs)

!begin EP


   ts_buf(:,:) = grid%ts_uru(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_uru(:,:),grid%ts_buf_size*grid%max_ts_locs)
   
   ts_buf(:,:) = grid%ts_uph(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_uph(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_upg(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_upg(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_uad(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_uad(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_uhd(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_uhd(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_und(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_und(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_uco(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_uco(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_ucu(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_ucu(:,:),grid%ts_buf_size*grid%max_ts_locs)
   
     ts_buf(:,:) = grid%ts_urs(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_urs(:,:),grid%ts_buf_size*grid%max_ts_locs)
   
   
   ts_buf(:,:) = grid%ts_vrv(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vrv(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_vph(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vph(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_vpg(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vpg(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_vad(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vad(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_vhd(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vhd(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_vnd(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vnd(:,:),grid%ts_buf_size*grid%max_ts_locs)
     
   ts_buf(:,:) = grid%ts_vco(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vco(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_vcu(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vcu(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_vrs(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_vrs(:,:),grid%ts_buf_size*grid%max_ts_locs)
   
   
   ts_buf(:,:) = grid%ts_trt(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_trt(:,:),grid%ts_buf_size*grid%max_ts_locs)
     
   ts_buf(:,:) = grid%ts_tbl(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_tbl(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_tra(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_tra(:,:),grid%ts_buf_size*grid%max_ts_locs)
   
   ts_buf(:,:) = grid%ts_tad(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_tad(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_thd(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_thd(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_tnd(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_tnd(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_trs(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_trs(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_tdi(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_tdi(:,:),grid%ts_buf_size*grid%max_ts_locs)
     
   ts_buf(:,:) = grid%ts_w(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_w(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_wrw(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_wrw(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_wbu(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_wbu(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_wad(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_wad(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_whd(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_whd(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_wnd(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_wnd(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_wco(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_wco(:,:),grid%ts_buf_size*grid%max_ts_locs)
      
   ts_buf(:,:) = grid%ts_wcu(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_wcu(:,:),grid%ts_buf_size*grid%max_ts_locs)
   
!end EP

#if (EM_CORE == 1)

   ts_buf(:,:) = grid%ts_mut(:,:) !EP
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_mut(:,:),grid%ts_buf_size*grid%max_ts_locs) !EP

   ts_buf(:,:) = grid%ts_muu(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_muu(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_muv(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_muv(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_glw(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_glw(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_gsw(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_gsw(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_hfx(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_hfx(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_lh(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_lh(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_clw(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_clw(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_rainc(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_rainc(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_rainnc(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_rainnc(:,:),grid%ts_buf_size*grid%max_ts_locs)
#endif

   ts_buf(:,:) = grid%ts_tsk(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_tsk(:,:),grid%ts_buf_size*grid%max_ts_locs)

   ts_buf(:,:) = grid%ts_tslb(:,:)
   CALL wrf_dm_min_reals(ts_buf(:,:),grid%ts_tslb(:,:),grid%ts_buf_size*grid%max_ts_locs)

   DEALLOCATE(ts_buf)
#endif

   IF ( wrf_dm_on_monitor() ) THEN

      iunit = get_unused_unit()
      IF ( iunit <= 0 ) THEN
         CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
      END IF

      DO i=1,grid%ntsloc_domain

         ix = grid%itsloc(i)
         iy = grid%jtsloc(i)

         OPEN(UNIT=iunit, FILE=TRIM(grid%ts_filename(i)), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

#if (EM_CORE == 1)
            WRITE(UNIT=iunit,FMT='(i2,f13.6,i5,i5,i5,1x,50(f13.5,1x))')  & !EP changed 14(f13.5....) to 47(f13.5...
                              grid%id, grid%ts_hour(n,i),                &
                              grid%id_tsloc(i), ix, iy,                  &
                              grid%ts_t(n,i),                            &
                              grid%ts_q(n,i),                            &
                              grid%ts_u(n,i),                            &
                              grid%ts_v(n,i),                            &
                              grid%ts_psfc(n,i),                         &
                              grid%ts_mut(n,i),                          & !EP
                              grid%ts_muu(n,i),                          &
                              grid%ts_muv(n,i),                          &
                              grid%ts_glw(n,i),                          &
                              grid%ts_gsw(n,i),                          &
                              grid%ts_hfx(n,i),                          &
                              grid%ts_lh(n,i),                           &
                              grid%ts_tsk(n,i),                          &
                              grid%ts_tslb(n,i),                         &
                              grid%ts_rainc(n,i),                        &
                              grid%ts_rainnc(n,i),                       &
                              grid%ts_clw(n,i),                          &
!begin EP

                              grid%ts_uru(n,i), &
                              grid%ts_uph(n,i), &
                              grid%ts_upg(n,i), &
                              grid%ts_uad(n,i), &
                              grid%ts_uhd(n,i), &
                              grid%ts_und(n,i), &
                              grid%ts_uco(n,i), &
                              grid%ts_ucu(n,i), &
                              grid%ts_urs(n,i), &

                              grid%ts_vrv(n,i), &
                              grid%ts_vph(n,i), &
                              grid%ts_vpg(n,i), &
                              grid%ts_vad(n,i), &
                              grid%ts_vhd(n,i), &
                              grid%ts_vnd(n,i), &
                              grid%ts_vco(n,i), &
                              grid%ts_vcu(n,i), &
                              grid%ts_vrs(n,i), &

                              grid%ts_trt(n,i), &
                              grid%ts_tbl(n,i), &
                              grid%ts_tra(n,i), &
                              grid%ts_tad(n,i), &
                              grid%ts_thd(n,i), &
                              grid%ts_tnd(n,i), &
                              grid%ts_tdi(n,i), &
                              grid%ts_trs(n,i), &

                              grid%ts_wrw(n,i), &
                              grid%ts_wbu(n,i), &
                              grid%ts_wad(n,i), &
                              grid%ts_whd(n,i), &
                              grid%ts_wnd(n,i), &
                              grid%ts_wco(n,i), &
                              grid%ts_wcu(n,i)

!end EP


#else
            WRITE(UNIT=iunit,FMT='(i2,f13.6,i5,i5,i5,1x,40(f13.5,1x))')   & !EP changed from 7(f to 40(f
                              grid%id, grid%ts_hour(n,i),                &
                              grid%id_tsloc(i), ix, iy,                  &
                              grid%ts_t(n,i),                            &
                              grid%ts_q(n,i),                            &
                              grid%ts_u(n,i),                            &
                              grid%ts_v(n,i),                            &
                              grid%ts_psfc(n,i),                         &
                              grid%ts_tsk(n,i),                          &
                              grid%ts_tslb(n,i),                         &
!begin EP

                              grid%ts_uru(n,i), &
                              grid%ts_uph(n,i), &
                              grid%ts_upg(n,i), &
                              grid%ts_uad(n,i), &
                              grid%ts_uhd(n,i), &
                              grid%ts_und(n,i), &
                              grid%ts_uco(n,i), &
                              grid%ts_ucu(n,i), &
                              grid%ts_urs(n,i), &

                              grid%ts_vrv(n,i), &
                              grid%ts_vph(n,i), &
                              grid%ts_vpg(n,i), &
                              grid%ts_vad(n,i), &
                              grid%ts_vhd(n,i), &
                              grid%ts_vnd(n,i), &
                              grid%ts_vco(n,i), &
                              grid%ts_vcu(n,i), &
                              grid%ts_vrs(n,i), &

                              grid%ts_trt(n,i), &
                              grid%ts_tbl(n,i), &
                              grid%ts_tra(n,i), &
                              grid%ts_tad(n,i), &
                              grid%ts_thd(n,i), &
                              grid%ts_tnd(n,i), &
                              grid%ts_tdi(n,i), &
                              grid%ts_trs(n,i), &

                              grid%ts_wrw(n,i), &
                              grid%ts_wbu(n,i), &
                              grid%ts_wad(n,i), &
                              grid%ts_whd(n,i), &
                              grid%ts_wnd(n,i), &
                              grid%ts_wco(n,i), &
                              grid%ts_wcu(n,i)

!end EP

#endif
         END DO
         CLOSE(UNIT=iunit)
         
         !Set write format for vertical profiles, depending on the highest model level of interest
#if (EM_CORE == 1)   
         profile_format = '(f13.6,1x,000(f13.5,1x))'
         k= LEN_TRIM(profile_format)
         WRITE(profile_format(12:14),'(I3.3)') grid%max_ts_level

         !Write u wind component profile to separate file
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for u wind component profiles
         WRITE(ts_profile_filename,'(A)') TRIM(grid%nametsloc(grid%id_tsloc(i)))//'.d00.TSS'
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-5:k-4),'(I2.2)') grid%id
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'UUU'
         
         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)           &
                              grid%ts_hour(n,i),                      &
                              grid%ts_u_profile(n,i,1:grid%max_ts_level)                                       
         END DO
         CLOSE(UNIT=iunit)


!begin EP U terms

         !Write uru wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for uru wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'URU'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_uru_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write uph wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for uph wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'UPH'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_uph_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write upg wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for upg wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'UPG'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_upg_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write uad wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for uad wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'UAD'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_uad_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write uhd wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for uhd wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'UHD'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_uhd_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write und wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for und wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'UND'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_und_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write uco wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for uco wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'UCO'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_uco_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write ucu wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for ucu wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'UCU'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_ucu_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write urs wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for urs wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'URS'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_urs_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


!end EP


         !Write v wind component profile to separate file
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for v wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'VVV'
         
         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             & 
                              grid%ts_v_profile(n,i,1:grid%max_ts_level)                                       
         END DO
         CLOSE(UNIT=iunit)


!begin EP V terms

         !Write vrv wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for vrv wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'VRV'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_vrv_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write vph wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for vph wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'VPH'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_vph_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write vpg wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for vpg wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'VPG'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_vpg_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write vad wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for vad wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'VAD'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_vad_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write vhd wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for vhd wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'VHD'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_vhd_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write vnd wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for vnd wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'VND'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_vnd_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write vco wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for vco wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'VCO'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_vco_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write vcu wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for vcu wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'VCU'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_vcu_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write vrs wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for vrs wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'VRS'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_vrs_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)

!end EP


         !Write geopotential height profile to separate file
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for geopotential height profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'PHH'
         
         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_gph_profile(n,i,1:grid%max_ts_level)                                       
         END DO
         CLOSE(UNIT=iunit)
        
 
         !Write potential temperature profile to separate file
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for potential temperature profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'THH'
         
         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)   &
                              grid%ts_hour(n,i),              &
                              grid%ts_th_profile(n,i,1:grid%max_ts_level)                                       
         END DO
         CLOSE(UNIT=iunit) 
       

!begin EP theta tendency terms 


         !Write trt profile to separate file
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for trt profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'TRT'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)   &
                              grid%ts_hour(n,i),              &
                              grid%ts_trt_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)



         !Write tbl profile to separate file
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for tbl profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'TBL'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')


         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)   &
                              grid%ts_hour(n,i),              &
                              grid%ts_tbl_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write tra profile to separate file
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for tra profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'TRA'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)   &
                              grid%ts_hour(n,i),              &
                              grid%ts_tra_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)



         !Write tad profile to separate file
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for tad profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'TAD'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)   &
                              grid%ts_hour(n,i),              &
                              grid%ts_tad_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)



         !Write thd profile to separate file
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for thd profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'THD'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)   &
                              grid%ts_hour(n,i),              &
                              grid%ts_thd_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)



         !Write tnd profile to separate file
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for tnd profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'TND'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)   &
                              grid%ts_hour(n,i),              &
                              grid%ts_tnd_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)



         !Write trs profile to separate file
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for trs profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'TRS'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)   &
                              grid%ts_hour(n,i),              &
                              grid%ts_trs_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write tdi profile to separate file
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for tdi profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'TDI'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)   &
                              grid%ts_hour(n,i),              &
                              grid%ts_tdi_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)

!end EP
!begin EP w terms

         !Write w wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for w wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'WWW'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_w_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)



         !Write wrw wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for wrw wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'WRW'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_wrw_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)



         !Write wbu wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for wbu wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'WBU'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_wbu_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write wad wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for wad wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'WAD'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_wad_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write whd wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for whd wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'WHD'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_whd_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write wnd wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for wnd wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'WND'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_wnd_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write wco wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for wco wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'WCO'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_wco_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)


         !Write wcu wind component profile to separate file EP
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for wcu wind component profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'WCU'

         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_wcu_profile(n,i,1:grid%max_ts_level)
         END DO
         CLOSE(UNIT=iunit)

!end EP




         !Write water vapor mixing ratio profile to separate file
         iunit = get_unused_unit()
            IF ( iunit <= 0 ) THEN
            CALL wrf_error_fatal('Error in write_ts: could not find a free Fortran unit.')
            END IF
         !Recreate filename for water vapor mixing ratio profiles
         k = LEN_TRIM(ts_profile_filename)
         WRITE(ts_profile_filename(k-2:k),'(A3)') 'QVV'
         
         OPEN(UNIT=iunit, FILE=TRIM(ts_profile_filename), STATUS='unknown', POSITION='append', FORM='formatted')

         DO n=1,grid%next_ts_time - 1

            WRITE(UNIT=iunit,FMT=profile_format)  &
                              grid%ts_hour(n,i),             &
                              grid%ts_qv_profile(n,i,1:grid%max_ts_level)                                       
         END DO
         CLOSE(UNIT=iunit)
#endif   

      END DO

   END IF

   grid%next_ts_time = 1

END SUBROUTINE write_ts


#if (EM_CORE == 1)
SUBROUTINE calc_p8w(grid, ix, iy, p8w, k_start, k_end)

   USE module_domain
   USE module_model_constants

   IMPLICIT NONE

   ! Arguments
   TYPE (domain), INTENT(IN) :: grid
   INTEGER, INTENT(IN) :: ix, iy, k_start, k_end
   REAL, DIMENSION(k_start:k_end), INTENT(OUT) :: p8w

   ! Local variables
   INTEGER :: k
   REAL    :: z0, z1, z2, w1, w2 
   REAL, DIMENSION(k_start:k_end)   :: z_at_w
   REAL, DIMENSION(k_start:k_end-1) :: z


   DO k = k_start, k_end
      z_at_w(k) = (grid%phb(ix,k,iy)+grid%ph_2(ix,k,iy))/g
   END DO

   DO k = k_start, k_end-1
      z(k) = 0.5*(z_at_w(k) + z_at_w(k+1))
   END DO

   DO k = k_start+1, k_end-1
      p8w(k) = grid%fnm(k)*(grid%p(ix,k,iy)+grid%pb(ix,k,iy)) + &
               grid%fnp(k)*(grid%p(ix,k-1,iy)+grid%pb(ix,k-1,iy))
   END DO

   z0 = z_at_w(k_start)
   z1 = z(k_start)
   z2 = z(k_start+1)
   w1 = (z0 - z2)/(z1 - z2)
   w2 = 1. - w1
   p8w(k_start) = w1*(grid%p(ix,k_start,iy)+grid%pb(ix,k_start,iy)) + &
                  w2*(grid%p(ix,k_start+1,iy)+grid%pb(ix,k_start+1,iy))

   z0 = z_at_w(k_end)
   z1 = z(k_end-1)
   z2 = z(k_end-2)
   w1 = (z0 - z2)/(z1 - z2)
   w2 = 1. - w1
   p8w(k_end) = exp(w1*log(grid%p(ix,k_end-1,iy)+grid%pb(ix,k_end-1,iy)) + &
                    w2*log(grid%p(ix,k_end-2,iy)+grid%pb(ix,k_end-2,iy)))

END SUBROUTINE calc_p8w
#endif
