subroutine BASGRA( PARAMS, MATRIX_WEATHER, &
                   CALENDAR_FERT, CALENDAR_NDEP, DAYS_HARVEST, &
                   NDAYS, NYEARS, SIZE_WEATHER, soilcn_option, NOUT, y)
!-------------------------------------------------------------------------------
! This is the BASic GRAss model originally written in MATLAB/Simulink by Marcel
! van Oijen, Mats Hoglind, Stig Morten Thorsen and Ad Schapendonk.
! 2011-07-13: Translation to FORTRAN by David Cameron and Marcel van Oijen.
! 2014-03-17: Extra category of tillers added
! 2014-04-03: Vernalization added
! 2014-04-03: Lower limit of temperature-driven leaf senescence no longer zero
! 2015-03-26: Introducing N, following example of BASFOR for the soil part.
!             Plant N is restricted to LV, ST, RT, not present in STUB and RES.
! 2015-09-24: More N-processes added
! 2016-12-09: Digestibility and fibre content added
!-------------------------------------------------------------------------------

use parameters_site
use parameters_plant
use environment
use resources
use soil
use plant
use basgra2yasso
use nasso
use set_params_mod

implicit none

!#integer, parameter :: nyears = 100
integer, dimension(100,2) :: DAYS_HARVEST
real                      :: PARAMS(120)
#ifdef weathergen  
  integer, parameter      :: NWEATHER =  7
#else
  integer, parameter      :: NWEATHER =  8
#endif
real                      :: MATRIX_WEATHER(SIZE_WEATHER,NWEATHER)
real   , dimension(100,3) :: CALENDAR_FERT, CALENDAR_NDEP
integer, dimension(100,2) :: DAYS_FERT    , DAYS_NDEP
real   , dimension(100)   :: NFERTV       , NDEPV

integer                   :: day, doy, i, NDAYS, NOUT, year, nyears, size_weather
integer                   :: soilcn_option
real                      :: y(NDAYS*NYEARS,NOUT)

! State variables plants
real    :: CLV, CLVD, CRES, CRT, CST, CSTUB, LAI, LT50, PHEN
real    :: ROOTD, TILG1, TILG2, TILV
integer :: VERN
real    :: YIELD, YIELD_LAST, YIELD_TOT
real    :: NRT, NSH

! Output variables constructed from plant state variables
real    :: DM, DMLV, DMRES, DMSH, DMST, DMSTUB, DM_MAX, TILTOT
real    :: NSH_DMSH
real    :: ENERGY_DM, F_ASH, F_PROTEIN, PROTEIN
real    :: NPP, NEE

! State variables soil
real    :: CLITT, CSOMF, CSOMS, DRYSTOR, Fdepth
real    :: NLITT, NSOMF, NSOMS, NMIN, O2, Sdepth
real    :: TANAER, WAL, WAPL, WAPS, WAS, WETSTOR
real    :: Nfert_TOT

! Intermediate and rate variables
real :: DeHardRate, DLAI, DLV, DPHEN, DRT, DSTUB, dTANAER, DTILV, EVAP, EXPLOR
real :: Frate, FREEZEL, FREEZEPL, GLAI, GLV, GPHEN, GRES, GRT, GST, GSTUB, GTILV, HardRate
real :: HARVLA, HARVLV, HARVPH, HARVRE, HARVST, HARVTILG2, INFIL, IRRIG, O2IN
real :: O2OUT, PackMelt, poolDrain, poolInfil, Psnow, reFreeze, RESMOB
real :: RGRTVG1, RROOTD, SnowMelt, THAWPS, THAWS, TILVG1, TILG1G2, TRAN, Wremain
real :: NCSHI, NCGSH, NCDSH, NCHARVSH, GNSH, DNSH, HARVNSH, GNRT, DNRT
real :: NSHmob, NSHmobsoil, Nupt
! yasso
real :: norg_runoff, soilc_runoff(num_c_pools)

real :: Ndep, Nfert

real :: F_DIGEST_DM, F_DIGEST_DMSH, F_DIGEST_LV, F_DIGEST_ST, F_DIGEST_WALL
real :: F_WALL_DM  , F_WALL_DMSH  , F_WALL_LV  , F_WALL_ST
real :: GRESSI, ALLOTOT
integer :: cycyear, total_day

! yasso state variables
real :: yasso_clim(num_clim_par)
real :: cumlittc, cumsoilr, cumphot, nmin_trunc
type(yasso_t) :: yasso_inst
character(len=10) :: soilcn_model
real :: ndemand_plant(1), nalloc_plant(1), nalloc_soil(2), ndemand_soil(2)

print *, 'Allocate for', size_weather, ndays, nyears
call alloc_environment(size_weather)

if (soilcn_option == 1) then
   soilcn_model = 'basgra'
else if (soilcn_option == 2) then
   soilcn_model = 'yasso'
   cumlittc = 0.0
   cumsoilr = 0.0
   cumphot = 0.0
   nmin_trunc = 0.0
else if (soilcn_option == 3) then
   ! basgra with yasso initial condition
   soilcn_model = 'basgra'
else
   print *, 'Bad soilcn_option', soilcn_option, ndays
   return
end if

! Parameters
call set_params(PARAMS)

! Calendar & weather
YEARI  = MATRIX_WEATHER(:,1)
DOYI   = MATRIX_WEATHER(:,2)
GRI    = MATRIX_WEATHER(:,3)
TMMNI  = MATRIX_WEATHER(:,4)
TMMXI  = MATRIX_WEATHER(:,5)
#ifdef weathergen  
  RAINI = MATRIX_WEATHER(:,6)
  PETI  = MATRIX_WEATHER(:,7)
#else
  VPI   = MATRIX_WEATHER(:,6)
  RAINI = MATRIX_WEATHER(:,7)
  WNI   = MATRIX_WEATHER(:,8)
#endif

! Calendars
DAYS_FERT  = CALENDAR_FERT (:,1:2)
DAYS_NDEP  = CALENDAR_NDEP (:,1:2)
NFERTV     = CALENDAR_FERT (:,3) * NFERTMULT
NDEPV      = CALENDAR_NDEP (:,3)

! Initial constants for plant state variables
CLV        = CLVI
CLVD       = CLVDI
CRES       = CRESI
CRT        = CRTI
CST        = CSTI
CSTUB      = CSTUBI ! set in parameters_plant
LAI        = LAII
LT50       = LT50I
NRT        = NCR * CRTI
  NCSHI    = NCSHMAX * (1-EXP(-K*LAII)) / (K*LAII)
NSH        = NCSHI * (CLVI+CSTI)
PHEN       = PHENI
ROOTD      = ROOTDM
TILG1      = TILTOTI *       FRTILGI *    FRTILGG1I
TILG2      = TILTOTI *       FRTILGI * (1-FRTILGG1I)
TILV       = TILTOTI * (1. - FRTILGI)
VERN       = 1
YIELD      = YIELDI
YIELD_LAST = YIELDI
YIELD_TOT  = YIELDI

Nfert_TOT  = 0
DM_MAX     = 0
NEE        = 0
NPP        = 0

! Initial constants for soil state variables
CLITT      = CLITT0
CSOMF      = CSOM0 * FCSOMF0
CSOMS      = CSOM0 * (1-FCSOMF0)
DRYSTOR    = DRYSTORI
Fdepth     = FdepthI
NLITT      = CLITT0 / CNLITT0
NSOMF      = (CSOM0 *    FCSOMF0)  / CNSOMF0
NSOMS      = (CSOM0 * (1-FCSOMF0)) / CNSOMS0
NMIN       = NMIN0
O2         = FGAS * ROOTDM * FO2MX * 1000./22.4
Sdepth     = SDEPTHI
TANAER     = TANAERI
WAL        = 1000. * ROOTDM * WCI
WAPL       = WAPLI
WAPS       = WAPSI
WAS        = WASI
WETSTOR    = WETSTORI

if (soilcn_model == 'yasso') then
   call init_yasso_from_files('parameters/yasso_parameters', 'initialisation/yasso_init', yasso_inst)
   call load_yasso_clim('weather/yasso_clim', yasso_clim)
   call map_yasso_to_basgra(yasso_inst, CLITT, CSOMF, CSOMS, NLITT, NSOMF, NSOMS)
else if (soilcn_option == 3) then
   call yasso_initcn2basgra('parameters/yasso_parameters', 'initialisation/yasso_init', &
        CLITT, CSOMF, CSOMS, NLITT, NSOMF, NSOMS)
end if


total_day = 0
do cycyear = 1, nyears
   do day = 1, NDAYS
      total_day = total_day + 1
      if (day == 1 .and. soilcn_model == 'basgra') then
         print *, 'Cycle:', cycyear
         print *, 'C:N, LITT', clitt / nlitt
         print *, 'C:N, SOMF', csomf / nsomf
         print *, 'C:N, SOMS', csoms / nsoms
      else if (day > 0 .and. soilcn_model == 'yasso') then
         print *, 'Year:', year, cycyear, day
         print *, 'C:N, SOM', get_totc(yasso_inst) / get_norg(yasso_inst)
         print *, 'C:N, AWEN', get_cn_awen(yasso_inst)
         print '(A, 5F10.2)', 'AWENH', yasso_inst%state(1:5)
         print *, 'Cumulative negative NMIN', nmin_trunc
         print *, 'NOHARV:', NOHARV, PHOT
      end if
      ! Environment
      call DDAYL          (doy)
      call set_weather_day(day,DRYSTOR,                    year,doy)
      call SoilWaterContent(Fdepth,ROOTD,WAL)
      call Physics        (DAVTMP,Fdepth,ROOTD,Sdepth,WAS, Frate)
      call MicroClimate   (doy,DRYSTOR,Fdepth,Frate,LAI,Sdepth,Tsurf,WAPL,WAPS,WETSTOR, &
           FREEZEPL,INFIL,PackMelt,poolDrain,poolInfil, &
           pSnow,reFreeze,SnowMelt,THAWPS,wRemain)
#ifdef weathergen  
      call PEVAPINPUT     (LAI)
#else
      call PENMAN         (LAI)
#endif
      ! Resources
      call Light          (DAYL,DTR,LAI,PAR)
      call EVAPTRTRF      (Fdepth,PEVAP,PTRAN,ROOTD,WAL,   EVAP,TRAN)
      call ROOTDG         (Fdepth,ROOTD,WAL,               EXPLOR,RROOTD)
      ! Soil
      call FRDRUNIR       (EVAP,Fdepth,Frate,INFIL,poolDRAIN,ROOTD,TRAN,WAL,WAS, &
           FREEZEL,IRRIG,THAWS)
      call O2status       (O2,ROOTD)
      ! Plant
      call Harvest        (CLV,CRES,CST,year,doy,DAYS_HARVEST,LAI,PHEN,TILG1,TILG2,TILV, &
           GSTUB,HARVLA,HARVLV,HARVPH,HARVRE,HARVST,HARVTILG2)
      call Biomass        (CLV,CRES,CST)
      call Phenology      (DAYL,PHEN,                      DPHEN,GPHEN)
      call Foliage1
      call LUECO2TM       (PARAV)
      call HardeningSink  (CLV,DAYL,doy,LT50,Tsurf)
      ! call Growth         (LAI,NSH,NMIN,CLV,CRES,CST,PARINT,TILG1,TILG2,TILV,TRANRF, &
      !                                                      GLV,GRES,GRT,GST,RESMOB,NSHmob)

      ! 
      call growth_demand(LAI,NSH,NMIN,CLV,CRES,CST,PARINT,TILG1,TILG2,TILV,TRANRF, &
           GLV,GRES,GRT,GST,RESMOB,NSHmob, ALLOTOT, GRESSI, ndemand_plant(1))
      if (soilcn_model == 'yasso') then
         call cn_decomp_demand(yasso_clim, delt, yasso_inst)
         ndemand_soil = get_n_demand(yasso_inst)
         call resolve_ndemand('demand_based', NMIN / TCNUPT, &
              ndemand_soil, ndemand_plant, &
              nalloc_soil, nalloc_plant)
         call cn_decomp_final(yasso_inst, nalloc_soil)
      else
         ndemand_soil = 0.0
         call resolve_ndemand('demand_based', NMIN / TCNUPT, &
              ndemand_soil, ndemand_plant, &
              nalloc_soil, nalloc_plant)
      end if
      call growth_final(nalloc_plant(1), GRES, GRT, GLV, GST, ALLOTOT, GRESSI, NSHmob)

      call PlantRespiration(FO2,RESPHARD)
      call Senescence     (CLV,CRT,CSTUB,LAI,LT50,PERMgas,TANAER,TILV,Tsurf, &
           DeHardRate,DLAI,DLV,DRT,DSTUB,dTANAER,DTILV,HardRate)
      call Foliage2       (DAYL,GLV,LAI,TILV,TILG1,TRANRF,Tsurf,VERN, &
           GLAI,GTILV,TILVG1,TILG1G2)
      ! Soil 2
      call O2fluxes       (O2,PERMgas,ROOTD,RplantAer,     O2IN,O2OUT)
      call N_fert         (year,doy,DAYS_FERT,NFERTV,      Nfert)
      call N_dep          (year,doy,DAYS_NDEP,NDEPV,       Ndep)

      if (soilcn_model == 'basgra') then
         call CNsoil         (ROOTD,RWA,WFPS,WAL,GRT,CLITT,CSOMF,NLITT,NSOMF,NSOMS,NMIN,CSOMS)
      else if (soilcn_model == 'yasso') then
         call cn_set_tend(yasso_clim, delt, yasso_inst)
         call cn_otherflux(ROOTD, RWA, WFPS, WAL, GRT, yasso_inst, nmin, &
              soilc_runoff, norg_runoff, Nleaching, NemissionNO, NemissionN2O, Nfixation)
         Nemission = NemissionNO + NemissionN2O
         ! abuse the basgra state variables with a rough mapping to Yasso pools
         call map_yasso_to_basgra(yasso_inst, CLITT, CSOMF, CSOMS, NLITT, NSOMF, NSOMS)
         RSOIL = get_soilresp(yasso_inst)
      else
         print *, 'Bad soilcn_model', soilcn_model
         return
      end if

      call Nplant         (NSHmob,CLV,CRT,CST,DLAI,DLV,DRT,GLAI,GLV,GRT,GST, &
           HARVLA,HARVLV,HARVST,LAI,NRT,NSH, &
           DNRT,DNSH,GNRT,GNSH,HARVNSH, &
           NCDSH,NCGSH,NCHARVSH,NSHmobsoil,Nupt)

      ! Extra variables
      DMLV      = CLV   / 0.45           ! Leaf dry matter; g m-2
      DMST      = CST   / 0.45           ! Stem dry matter; g m-2
      DMSTUB    = CSTUB / 0.45
      DMRES     = CRES  / 0.40
      DMSH      = DMLV + DMST + DMRES
      DM        = DMSH + DMSTUB
      TILTOT    = TILG1 + TILG2 + TILV

      NSH_DMSH  = NSH / DMSH             ! N content in shoot DM; g N g-1 DM

      PROTEIN   = NSH * 6.25             ! Crude protein; g m-2
      F_PROTEIN = PROTEIN / DMSH         ! Crude protein in shoot dry matter; g g-1   
      F_ASH	    = 0.069 + 0.14*F_PROTEIN

      call Digestibility  (DM,DMLV,DMRES,DMSH,DMST,DMSTUB,PHEN, &
           F_WALL_DM,F_WALL_DMSH,F_WALL_LV,F_WALL_ST, &
           F_DIGEST_DM,F_DIGEST_DMSH,F_DIGEST_LV,F_DIGEST_ST,F_DIGEST_WALL)

      !================
      ! Outputs
      !================
      y(total_day, 1) = year + (doy-0.5)/366 ! "Time" = Decimal year (approximation)
      y(total_day, 2) = year
      y(total_day, 3) = doy
      y(total_day, 4) = DAVTMP

      y(total_day, 5) = CLV
      y(total_day, 6) = CLVD
      y(total_day, 7) = YIELD_LAST
      y(total_day, 8) = CRES
      y(total_day, 9) = CRT
      y(total_day,10) = CST
      y(total_day,11) = CSTUB
      y(total_day,12) = DRYSTOR
      y(total_day,13) = Fdepth
      y(total_day,14) = LAI
      y(total_day,15) = LT50
      y(total_day,16) = O2
      y(total_day,17) = PHEN
      y(total_day,18) = ROOTD
      y(total_day,19) = Sdepth
      y(total_day,20) = TANAER
      y(total_day,21) = TILG1 + TILG2          ! "TILG"
      y(total_day,22) = TILV
      y(total_day,23) = WAL
      y(total_day,24) = WAPL
      y(total_day,25) = WAPS
      y(total_day,26) = WAS
      y(total_day,27) = WETSTOR

      y(total_day,28) = DM                     ! "DM"      = Aboveground dry matter in g m-2
      y(total_day,29) = DMRES / DM             ! "RES"     = Reserves in g g-1 aboveground dry matter
      y(total_day,30) = LERG                   !
      y(total_day,31) = NELLVG                 !
      y(total_day,32) = RLEAF                  !
      y(total_day,33) = LAI / DMLV             ! "SLA"     = m2 leaf area g-1 dry matter leaves
      y(total_day,34) = TILTOT                 ! "TILTOT"  = Total tiller number in # m-2
      y(total_day,35) = (TILG1+TILG2) / TILTOT ! "FRTILG"  = Fraction of tillers that is generative
      y(total_day,36) =  TILG1        / TILTOT ! "FRTILG1" = Fraction of tillers that is in TILG1
      y(total_day,37) =        TILG2  / TILTOT ! "FRTILG2" = Fraction of tillers that is in TILG2
      y(total_day,38) = RDRT
      y(total_day,39) = VERN

      y(total_day,40) = CLITT           ! g C m-2
      y(total_day,41) = CSOMF           ! g C m-2
      y(total_day,42) = CSOMS           ! g C m-2
      y(total_day,43) = NLITT           ! g N m-2
      y(total_day,44) = NSOMF           ! g N m-2
      y(total_day,45) = NSOMS           ! g N m-2
      y(total_day,46) = NMIN            ! g N m-2
      y(total_day,47) = Rsoil           ! g C m-2 d-1
      y(total_day,48) = NemissionN2O    ! g N m-2 d-1
      y(total_day,49) = NemissionNO     ! g N m-2 d-1
      y(total_day,50) = Nfert           ! g N m-2 d-1
      y(total_day,51) = Ndep            ! g N m-2 d-1
      y(total_day,52) = RWA             ! -
      y(total_day,53) = NSH             ! g N m-2
      y(total_day,54) = GNSH
      y(total_day,55) = DNSH
      y(total_day,56) = HARVNSH 
      y(total_day,57) = NSH / (CLV+CST) ! - "NCSH"
      y(total_day,58) = NCGSH
      y(total_day,59) = NCDSH
      y(total_day,60) = NCHARVSH
      y(total_day,61) = fNgrowth
      y(total_day,62) = RGRTV
      y(total_day,63) = FSPOT
      y(total_day,64) = RESNOR
      y(total_day,65) = TV2TIL
      y(total_day,66) = NSHNOR
      y(total_day,67) = KNMAX
      y(total_day,68) = KN

      y(total_day,69) = DMLV
      y(total_day,70) = DMST
      y(total_day,71) = NSH_DMSH
      y(total_day,72) = Nfert_TOT
      y(total_day,73) = YIELD_TOT
      y(total_day,74) = DM_MAX

      y(total_day,75) = F_PROTEIN
      y(total_day,76) = F_ASH

      y(total_day,77) = F_WALL_DM
      y(total_day,78) = F_WALL_DMSH
      y(total_day,79) = F_WALL_LV
      y(total_day,80) = F_WALL_ST
      y(total_day,81) = F_DIGEST_DM
      y(total_day,82) = F_DIGEST_DMSH
      y(total_day,83) = F_DIGEST_LV
      y(total_day,84) = F_DIGEST_ST
      y(total_day,85) = F_DIGEST_WALL

      y(total_day,86) = RDRS
      y(total_day,87) = RAIN
      y(total_day,88) = Nleaching

      y(total_day,89) = NSHmob
      y(total_day,90) = NSHmobsoil
      y(total_day,91) = Nfixation
      y(total_day,92) = Nupt
      y(total_day,93) = Nmineralisation

      y(total_day,94) = NSOURCE
      y(total_day,95) = NSINK

      y(total_day,96) = NRT
      y(total_day,97) = NRT / CRT

      y(total_day,98) = rNLITT
      y(total_day,99) = rNSOMF

      y(total_day,100) = DAYL
      y(total_day,101) = GSTUB+GRT+GST+GLV ! GTOT
      y(total_day,102) = DNRT
      y(total_day,103) = norg_runoff
      y(total_day,104) = GRT
      y(total_day,105) = NPP
      y(total_day,106) = NEE
      y(total_day,107) = RESMOB
      y(total_day,108) = DLV+DSTUB ! FLITTC_LEAF
      y(total_day,109) = DRT       ! FLITTC_ROOT
      y(total_day,110) = DNSH ! FLITTN_LEAF
      y(total_day,111) = DNRT ! FLITTN_ROOT
      

      ! State equations plants
      CLV     = CLV     + GLV   - DLV    - HARVLV
      CLVD    = CLVD            + DLV
      CRES    = CRES    + GRES  - RESMOB - HARVRE
      CRT     = CRT     + GRT   - DRT
      CST     = CST     + GST           - HARVST
      CSTUB   = CSTUB   + GSTUB - DSTUB
      LAI     = LAI     + GLAI - DLAI   - HARVLA
      LT50    = LT50    + DeHardRate - HardRate
      PHEN    = min(1., PHEN + GPHEN - DPHEN - HARVPH)
      ROOTD   = ROOTD   + RROOTD
      TILG1   = TILG1           + TILVG1 - TILG1G2
      TILG2   = TILG2                    + TILG1G2 - HARVTILG2
      TILV    = TILV    + GTILV - TILVG1           - DTILV
      if((LAT>0).AND.(doy==305)) VERN = 0  
      if((LAT<0).AND.(doy==122)) VERN = 0  
      if(DAVTMP<TVERN)           VERN = 1
      YIELD     = (HARVLV + HARVST*HAGERE)/0.45 + HARVRE/0.40
      if(YIELD>0) YIELD_LAST = YIELD
      YIELD_TOT = YIELD_TOT + YIELD

      NRT       = NRT   + GNRT - DNRT
      NSH       = NSH   + GNSH - DNSH - HARVNSH - NSHmob

      Nfert_TOT = Nfert_TOT + Nfert
      DM_MAX    = max( DM, DM_MAX )


      ! State equations soil

      ! Default BASGRA C & N
      if (soilcn_model == 'basgra') then
         CLITT   = CLITT + DLV + DSTUB              - rCLITT - dCLITT
         CSOMF   = CSOMF + DRT         + dCLITTsomf - rCSOMF - dCSOMF
         CSOMS   = CSOMS               + dCSOMFsoms          - dCSOMS
         NLITT   = NLITT + DNSH             - rNLITT - dNLITT
         NSOMF   = NSOMF + DNRT + NLITTsomf - rNSOMF - dNSOMF 
         NSOMS   = NSOMS        + NSOMFsoms          - dNSOMS
      else if (soilcn_model == 'yasso') then
         call cn_state_update(yasso_inst, &
              get_littc_awenh(DLV+DSTUB, DRT) - soilc_runoff, &
              DNSH + DNRT-norg_runoff)
         !call cn_state_update(yasso_inst, &
         !     get_littc_awenh(0.0, DRT) - soilc_runoff, &
         !     DNRT-norg_runoff)

         nMineralisation = get_netmin_act(yasso_inst)
         cumlittc = cumlittc + sum(get_littc_awenh(DLV+DSTUB, DRT))
         cumsoilr = cumsoilr + get_soilresp(yasso_inst)
         cumphot = cumphot + phot
      end if

      ! Mineral N the same in basgra and yasso
      NMIN = NMIN  + Ndep + Nfert + Nmineralisation + Nfixation + NSHmobsoil &
           - Nupt - Nleaching - Nemission
      if (NMIN < 0) then
         nmin_trunc = nmin_trunc - NMIN
      end if
      NMIN    = max(0.,NMIN)

      NPP = GLV + GRES - RESMOB + GRT + GST + GSTUB
      NEE = NPP - RSOIL

      DRYSTOR = DRYSTOR + reFreeze + Psnow - SnowMelt
      Fdepth  = Fdepth  + Frate
      O2      = O2      + O2IN - O2OUT
      Sdepth  = Sdepth  + Psnow/RHOnewSnow - PackMelt
      TANAER  = TANAER  + dTANAER
      WAL     = WAL  + THAWS  - FREEZEL  + poolDrain + INFIL +EXPLOR+IRRIG-DRAIN-RUNOFF-EVAP-TRAN
      WAPL    = WAPL + THAWPS - FREEZEPL + poolInfil - poolDrain
      WAPS    = WAPS - THAWPS + FREEZEPL
      WAS     = WAS  - THAWS  + FREEZEL
      WETSTOR = WETSTOR + Wremain - WETSTOR

   enddo
end do

if (soilcn_model == 'yasso') then
   call store_yasso_state(yasso_inst, 'yasso.state')
end if

call dealloc_environment()

end subroutine BASGRA

