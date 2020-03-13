subroutine BASGRA( PARAMS, MATRIX_WEATHER, weather_yasso, yasso_init, yasso_params, &
                   CALENDAR_FERT, CALENDAR_NDEP, DAYS_HARVEST, soilcn_option, if_weathergen, &
                   NDAYS, NYEARS, NWEATHER, SIZE_WEATHER, size_calendar, size_yasso_init, size_yasso_params, NOUT, y)
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
use yassobase
use nmodel
use cascademodel
use set_params_mod
use yassocore, only : num_c_pools, num_clim_par, num_yasso_param => num_parameters

implicit none

real, intent(in)          :: PARAMS(113) ! Parameters for BASGRA
real, intent(in)          :: MATRIX_WEATHER(SIZE_WEATHER,NWEATHER) ! (number of days, number of parameters)
real, intent(in)          :: weather_yasso(5, SIZE_WEATHER)        ! (number of parameters, number of days)
real, intent(in)          :: yasso_init(size_yasso_init)           ! (yasso initial state)
real, intent(in)          :: yasso_params(size_yasso_params)       ! yasso parameters
! CALENDAR_FERT: (ind_fert, year - 0 means any; doy; amount in gN/m2; type (1.0=mineral, 2.0=soluble organic); C:N ratio)
real, intent(in), dimension(size_calendar,5) :: CALENDAR_FERT
real, intent(in), dimension(size_calendar,3) :: CALENDAR_NDEP      ! (ind_entry, year-doy-rate)
integer, intent(in), dimension(size_calendar,3) :: DAYS_HARVEST    ! (ind_harvest, year - 0 means any; doy; flag(0=harvest, >0 = cut only))

! Possible values for soilcn_option:
integer, parameter :: soilcn_basgra = 1, &
     soilcn_nmodel1 = 2, &
     soilcn_nmodel2cls = 3, &
     soilcn_basgra_nm1 = 4, &
     soilcn_basgra_nm2c = 5, &
     soilcn_nmodel1_soil = 6

integer, intent(in)                   :: soilcn_option             ! use the Basgra default soil C & N vs. use Yasso. See above.
logical, intent(in)                   :: if_weathergen             ! "weather handling mode" for BASGRA
! number of days for each cycled year (365 for full years when nyears > 1, > 365 when running several real years):
integer, intent(in)                   :: NDAYS                     
integer, intent(in)                   :: nyears                    ! number of years to cycle over, set to 1 when running "real" calendar years. 
integer, intent(in)                   :: size_weather              ! number of weather records
integer, intent(in)                   :: size_calendar             ! number of fert/ndep/harvest calendar records
integer, intent(in)                   :: NWEATHER                  ! number of weather parameters (in each record)
integer, intent(in)                   :: size_yasso_init, size_yasso_params
integer, intent(in)                   :: NOUT                      ! number of outputs
real, intent(out)                     :: y(NDAYS*NYEARS,NOUT)      ! output

integer :: day, i, year, doy
integer, dimension(size_calendar,2) :: DAYS_FERT, DAYS_NDEP
real, dimension(size_calendar)   :: NDEPV

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
real :: harv_c_to_litt, harv_n_to_litt


real :: Ndep, Nfert_min, Nfert_org, Cfert, fertflag

real :: F_DIGEST_DM, F_DIGEST_DMSH, F_DIGEST_LV, F_DIGEST_ST, F_DIGEST_WALL
real :: F_WALL_DM  , F_WALL_DMSH  , F_WALL_LV  , F_WALL_ST
real :: GRESSI, ALLOTOT
integer :: cycyear, total_day
logical :: if_cut_only, wx_found

! yasso state variables
real :: yasso_clim(num_clim_par)     ! The yasso drivers. Content depends on the driver option.
real :: cumlittc, cumsoilr, cumphot, nmin_trunc
class(base_yasso_t), allocatable :: yasso_inst

real :: ndemand_plant(1)             ! The plant N demand evaluated from the potential photolysis
real :: nalloc_plant(1)              ! The N allocated for plant depending on N availability and soil demand.
! The soil N demand/alloc may have several components (= several processes that release/fix N).
real, allocatable :: ndemand_soil(:) ! The N demand of litter decomposition. Negative means net mineralization. Alway negative for the Basgra C/N.
real, allocatable :: nalloc_soil(:)  ! The N allocated for litter decomposition (if required).

logical :: get_ic_from_yasso, need_yasso

! yasso flux variables
real, allocatable :: norg_runoff(:), soilc_runoff(:), nflux(:), cflux(:)

! Parameters
call set_params(PARAMS)

! Calendar & weather
YEARI  = MATRIX_WEATHER(:,1)
DOYI   = MATRIX_WEATHER(:,2)
GRI    = MATRIX_WEATHER(:,3)
TMMNI  = MATRIX_WEATHER(:,4)
TMMXI  = MATRIX_WEATHER(:,5)
if (NWEATHER == 7) then
   RAINI = MATRIX_WEATHER(:,6)
   PETI  = MATRIX_WEATHER(:,7)
else
   VPI   = MATRIX_WEATHER(:,6)
   RAINI = MATRIX_WEATHER(:,7)
   WNI   = MATRIX_WEATHER(:,8)
end if

! Calendars
DAYS_FERT  = CALENDAR_FERT (:,1:2)
DAYS_NDEP  = CALENDAR_NDEP (:,1:2)
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




select case(soilcn_option)
case (soilcn_basgra)
   get_ic_from_yasso = .false.
   allocate(soilc_runoff(1), norg_runoff(1), ndemand_soil(1), nalloc_soil(1))
   soilc_runoff = 0.0
   norg_runoff = 0.0
   
case (soilcn_nmodel1, soilcn_nmodel2cls, soilcn_nmodel1_soil)
   cumlittc = 0.0
   cumsoilr = 0.0
   cumphot = 0.0
   nmin_trunc = 0.0
   call yasso_factory(yasso_inst, soilcn_option, yasso_init, yasso_params)
   allocate(soilc_runoff(yasso_inst%get_csize()), norg_runoff(yasso_inst%get_nsize()))
   allocate(cflux(yasso_inst%get_csize()), nflux(yasso_inst%get_nsize()), &
        ndemand_soil(yasso_inst%get_demand_size()), nalloc_soil(yasso_inst%get_demand_size()))
   
case (soilcn_basgra_nm1)
   call yasso_factory(yasso_inst, soilcn_nmodel1, yasso_init, yasso_params)
   call yasso_inst%map_to_basgra(CLITT, CSOMF, CSOMS, NLITT, NSOMF, NSOMS)
   deallocate(yasso_inst)
   allocate(soilc_runoff(1), norg_runoff(1))
   
case default
   print *, 'Bad soilcn_option'
   return
end select

total_day = 0
do cycyear = 1, nyears
   do day = 1, NDAYS
      total_day = total_day + 1
      if (day == 1 .and. .not. allocated(yasso_inst)) then
         print *, 'Cycle:', cycyear
      else if (day > 0 .and. allocated(yasso_inst)) then
         print *, 'Year:', year, cycyear, day
         if (yasso_inst%get_totc() / yasso_inst%get_norg() < 0 .or. yasso_inst%get_totc() / yasso_inst%get_norg() > 1000) then
            print *, 'Inconsistent yasso state'
            call yasso_inst%report()
            stop
         end if
      end if
      
      ! Environment
      call DDAYL          (doy)
      call set_weather_day(day, DRYSTOR, year, doy, if_weathergen)
      call SoilWaterContent(Fdepth,ROOTD,WAL)
      call Physics        (DAVTMP,Fdepth,ROOTD,Sdepth,WAS, Frate)
      call MicroClimate   (doy,DRYSTOR,Fdepth,Frate,LAI,Sdepth,Tsurf,WAPL,WAPS,WETSTOR, &
           FREEZEPL,INFIL,PackMelt,poolDrain,poolInfil, &
           pSnow,reFreeze,SnowMelt,THAWPS,wRemain)

      if (if_weathergen) then
         call PEVAPINPUT     (LAI)
      else
         call PENMAN         (LAI)
      end if

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
           GSTUB,HARVLA,HARVLV,HARVPH,HARVRE,HARVST,HARVTILG2,if_cut_only)
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
      if (allocated(yasso_inst)) then
         if (soilcn_option == soilcn_nmodel1_soil) then
            call get_soil_clim(yasso_clim)
         else
            call get_daily_clim(year, doy, weather_yasso, yasso_clim, wx_found)
            if (.not. wx_found) then
               print *, 'Failed to find weather for yasso: year, doy', year, doy
               return
            end if
         end if
         call yasso_inst%decomp_demand(yasso_clim, delt)
         ndemand_soil = yasso_inst%get_n_demand()
         call resolve_ndemand('demand_based', NMIN / TCNUPT, &
              ndemand_soil, ndemand_plant, &
              nalloc_soil, nalloc_plant)
         call yasso_inst%decomp_final(yasso_clim, delt, nalloc_soil)
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
      call N_fert         (year, doy, CALENDAR_FERT, Nfert_min, Nfert_org, Cfert, fertflag)
      call N_dep          (year,doy,DAYS_NDEP,NDEPV,       Ndep)
      
      if (allocated(yasso_inst)) then
         call eval_otherfluxes(ROOTD, RWA, WFPS, WAL, GRT, yasso_inst, nmin, &
              soilc_runoff, norg_runoff, Nleaching, NemissionNO, NemissionN2O, Nfixation)
         Nemission = NemissionNO + NemissionN2O
         ! abuse the basgra state variables with a rough mapping to Yasso pools
         call yasso_inst%map_to_basgra(CLITT, CSOMF, CSOMS, NLITT, NSOMF, NSOMS)
         RSOIL = yasso_inst%get_soilresp()
      else ! basgra
         call CNsoil         (ROOTD,RWA,WFPS,WAL,GRT,CLITT,CSOMF,NLITT,NSOMF,NSOMS,NMIN,CSOMS)
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
      y(total_day,50) = Nfert_min       ! g N m-2 d-1
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
      y(total_day,103) = sum(norg_runoff)
      y(total_day,104) = GRT
      y(total_day,105) = NPP
      y(total_day,106) = NEE
      y(total_day,107) = RESMOB
      y(total_day,108) = DLV+DSTUB ! FLITTC_LEAF
      y(total_day,109) = DRT       ! FLITTC_ROOT
      y(total_day,110) = DNSH ! FLITTN_LEAF
      y(total_day,111) = DNRT ! FLITTN_ROOT
      y(total_day,112) = VP 

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

      ! Yield computed only if harvest is removed
      if (.not. if_cut_only) then
         ! why HAGERE here?
         YIELD     = (HARVLV + HARVST*HAGERE)/0.45 + HARVRE/0.40
         if(YIELD>0) YIELD_LAST = YIELD
         YIELD_TOT = YIELD_TOT + YIELD
         harv_c_to_litt = 0.0
         harv_n_to_litt = 0.0
      else
         YIELD = 0.0
         harv_c_to_litt = HARVLV + HARVST*HAGERE + HARVRE
         harv_n_to_litt = HARVNSH 
      end if
         
      NRT       = NRT   + GNRT - DNRT
      NSH       = NSH   + GNSH - DNSH - HARVNSH - NSHmob

      Nfert_TOT = Nfert_TOT + Nfert_min + Nfert_org
      DM_MAX    = max( DM, DM_MAX )

      ! State equations soil

      if (allocated(yasso_inst)) then
         nflux = yasso_inst%litt_n_to_nflux(DNSH + DNRT + harv_n_to_litt + nfert_org)
         cflux = yasso_inst%litt_awenh_to_cflux(get_littc_awenh(DLV+DSTUB+harv_c_to_litt, DRT) &
                                                + get_fertc_awenh(fertflag, Cfert))
         call yasso_inst%update_state(cflux - soilc_runoff, nflux-norg_runoff)
         nMineralisation = yasso_inst%get_netmin_act()
         cumlittc = cumlittc + sum(get_littc_awenh(DLV+DSTUB+harv_c_to_litt, DRT))
         cumsoilr = cumsoilr + yasso_inst%get_soilresp()
         cumphot = cumphot + phot
         call yasso_inst%report()
      else
         ! Default BASGRA C & N
         CLITT   = CLITT + DLV + DSTUB + harv_c_to_litt + Cfert - rCLITT - dCLITT
         CSOMF   = CSOMF + DRT         + dCLITTsomf - rCSOMF - dCSOMF
         CSOMS   = CSOMS               + dCSOMFsoms          - dCSOMS
         NLITT   = NLITT + DNSH        + harv_n_to_litt + Nfert_org      - rNLITT - dNLITT
         NSOMF   = NSOMF + DNRT + NLITTsomf - rNSOMF - dNSOMF 
         NSOMS   = NSOMS        + NSOMFsoms          - dNSOMS
      end if

      ! Mineral N the same in basgra and yasso
      NMIN = NMIN  + Ndep + Nfert_min + Nmineralisation + Nfixation + NSHmobsoil &
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

print *, 'Done time loop'

if (allocated(yasso_inst)) then
   print *, 'store state'
   call yasso_inst%store_state('yasso.state')
   print *, 'done'
end if
print *, 'Dealloc'
call dealloc_environment()

contains

  subroutine yasso_factory(inst, soilcn_model, state_init, yasso_params)
    ! Instantiate a Yasso model depending on the soilcn_model option.
    use yassocore, only : num_parameters
    class(base_yasso_t), allocatable :: inst
    integer, intent(in) :: soilcn_model
    real, intent(in) :: state_init(:)
    real, intent(in) :: yasso_params(:)
    
    real :: param(num_parameters)

    if (size(yasso_params) < num_parameters) then
       print *, 'yasso_params too small'
       stop
    end if

    select case(soilcn_model)
    case(soilcn_nmodel1)
       allocate(nmodel1::inst)
    case(soilcn_nmodel1_soil)
       allocate(nmodel1_soil::inst)
    case(soilcn_nmodel2cls)
       allocate(nmodel2cls::inst)
    case default
       print *, 'Bad soilcn_model', soilcn_model
       stop
    end select

    call inst%init(yasso_params)
    call inst%set_state(state_init)

    !call inst%init_from_files('parameters/yasso_parameters', filename_inicn)
       
  end subroutine yasso_factory
  
  subroutine eval_otherfluxes(ROOTD, RWA, WFPS, WAL, GCR, yasso_inst, nmin, &
       soilc_runoff, norg_runoff, nmin_leach, n_emis_no, n_emis_n2o, nfixation)
    ! Evaluate leaching and gaseous N fluxes when using Yasso for soil
    !
    use parameters_site, only : KNFIX, RRUNBULK, KNEMIT, RFN2O, WFPS50N2O, RNLEACH
    use soil, only : DRAIN, RUNOFF
    real, intent(in) :: ROOTD, RWA, WFPS, WAL, GCR
    class(base_yasso_t), intent(in) :: yasso_inst
    real, intent(in) :: nmin
    
    real, intent(out) :: soilc_runoff(:)    ! C runoff from yasso C state
    real, intent(out) :: norg_runoff(:)     ! N runoff from yasso N state
    real, intent(out) :: nmin_leach         ! N leaching from NMIN
    real, intent(out) :: n_emis_no          ! NOx emission
    real, intent(out) :: n_emis_n2o         ! N2O emission
    real, intent(out) :: nfixation          ! N fixation

    real :: fN2O, Nemission
    real :: cstate(yasso_inst%get_csize())
    real :: nstate(yasso_inst%get_nsize())
    
    call yasso_inst%cstate(get=cstate)
    soilc_runoff = (cstate*runoff / rootd)  * RRUNBULK * 0.001
    call yasso_inst%nstate(get=nstate)
    norg_runoff = (nstate*runoff / rootd) * RRUNBULK * 0.001

    nfixation       = gCR * KNFIX
    ! Nleaching       = (NMIN*RNLEACH*DRAIN) / WAL
    if ((WAL > 0.) .and. (NMIN > 0.)) then
       nmin_leach      = (NMIN*RNLEACH*DRAIN) / WAL
    else
       nmin_leach       = 0.0
    end if

    Nemission       = NMIN * KNEMIT * RWA
    fN2O            = 1. / (1. + exp(-RFN2O*(WFPS-WFPS50N2O)))
    n_emis_n2o    = Nemission *     fN2O
    n_emis_no     = Nemission * (1.-fN2O)
  end subroutine eval_otherfluxes

  subroutine get_soil_clim(clim)
    real, intent(out) :: clim(3)

    clim(1) = Tsurf
    clim(2) = RAIN*365.0
    clim(3) = 0.001 * WAL / (ROOTD*WCST) ! WAL is in mm and ROOTD in m
    
  end subroutine get_soil_clim
  
end subroutine BASGRA

