module yassobase
  use nasso
  implicit none
  private

  public resolve_ndemand
  public get_fertc_awenh
  public get_littc_awenh
  public load_yasso_clim
  public load_params
  public init_from_files
  public get_daily_clim
  
  integer, parameter :: fileunit_def = 30
  
  real, parameter :: litterfract_leaf(num_c_pools) = (/0.23, 0.68, 0.0, 0.09, 0.0/)
  real, parameter :: litterfract_fineroot(num_c_pools) = (/0.44, 0.34, 0.0, 0.22, 0.0/)

  public base_yasso_t
  type, abstract :: base_yasso_t
     integer  :: demand_size
     integer :: cstate_size
     integer :: nstate_size
     real :: param(num_parameters)

   contains
     
     procedure, pass(this) :: get_demand_size
     procedure, pass(this) :: get_csize
     procedure, pass(this) :: get_nsize
     
     procedure(get_n_demand), pass(this), deferred :: get_n_demand
     procedure(get_norg), pass(this), deferred :: get_norg
     procedure(get_netmin_act), pass(this), deferred :: get_netmin_act
     procedure(get_totc), pass(this), deferred :: get_totc
     procedure(get_cn_awen), pass(this), deferred :: get_cn_awen
     procedure(get_soilresp), pass(this), deferred :: get_soilresp
     procedure(map_to_basgra), pass(this), deferred :: map_to_basgra
     procedure(store_state), pass(this), deferred :: store_state
     procedure(load_state), pass(this), deferred :: load_state     
     procedure(decomp_demand), pass(this), deferred :: decomp_demand
     procedure(decomp_final), pass(this), deferred :: decomp_final
     !procedure(eval_tend), pass(this), deferred :: eval_tend
     procedure(update_state), pass(this), deferred :: update_state
     procedure(force_n), pass(this), deferred :: force_n
     procedure(cmatrix), pass(this), deferred :: cmatrix
     procedure(cstate), pass(this), deferred :: cstate
     procedure(nstate), pass(this), deferred :: nstate
     procedure(init), pass(this), deferred :: init
     procedure(litt_awenh_to_cflux), pass(this), deferred :: litt_awenh_to_cflux
     procedure(litt_n_to_nflux), pass(this), deferred :: litt_n_to_nflux
     procedure, pass(this) :: init_from_files
     procedure, pass(this) :: get_params
     procedure, pass(this) :: set_params
     procedure, pass(this) :: set_sizes
     procedure, pass(this) :: delete
     procedure, pass(this) :: test
     procedure(report), pass(this), deferred :: report
     !procedure, pass(this) :: init

  end type base_yasso_t

  interface load_yasso_clim
     module procedure load_yasso_clim_1d
     module procedure load_yasso_clim_2d
  end interface load_yasso_clim
  
  abstract interface

     !**********************************************************************************
     !
     ! Abstract methods for base_yasso_t
     !
     !**********************************************************************************
     
     subroutine report(this)
       import :: base_yasso_t
       class(base_yasso_t), intent(in) :: this
     end subroutine report
     
     function get_n_demand(this) result(demand)
       ! Return the N required to realize the potential decomposition. The model may have
       ! several N-consuming steps, which are here considered separately.
       import :: base_yasso_t
       class(base_yasso_t), intent(in) :: this
       real :: demand(this%demand_size)
     end function get_n_demand

     function get_norg(this) result(norg)
       ! Return the total organic N in SOM.
       import :: base_yasso_t
       class(base_yasso_t), intent(in) :: this
       real :: norg
     end function get_norg
     
     function get_netmin_act(this) result(net)
       ! Return the actual net mineralization flux after adjusting for N limitation.
       import :: base_yasso_t
       class(base_yasso_t), intent(in) :: this
       real :: net
     end function get_netmin_act

     function get_totc(this) result(totc)
       ! Return the total soil C.
       import :: base_yasso_t
       class(base_yasso_t), intent(in) :: this
       real :: totc
     end function get_totc

     function get_cn_awen(this) result(ratio)
       ! Return the C:N ratio in the AWEN pools, excluding H.
       import :: base_yasso_t
       class(base_yasso_t), intent(in) :: this
       real :: ratio
     end function get_cn_awen

     function get_soilresp(this) result(resp)
       ! Return the heterotrophic respiration flux.
       import :: base_yasso_t
       class(base_yasso_t), intent(in) :: this
       real :: resp
     end function get_soilresp

     subroutine map_to_basgra(this, clitt, csomf, csoms, nlitt, nsomf, nsoms)
       ! Map the C and N in the model pools to the BASGRA-N litter and SOM C and N pools.
       import :: base_yasso_t
       class(base_yasso_t), intent(in) :: this
       real, intent(out) :: clitt, csomf, csoms, nlitt, nsomf, nsoms
     end subroutine map_to_basgra

     subroutine store_state(this, filename, fileunit)
       ! Store the model state into a file.
       import :: base_yasso_t
       class(base_yasso_t), intent(in) :: this
       character(len=*), intent(in) :: filename
       integer, intent(in), optional :: fileunit
     end subroutine store_state

     subroutine decomp_demand(this, climate, timestep_days)
       ! Evaluate the potential decomposition and store it internally.
       import :: base_yasso_t, num_clim_par
       class(base_yasso_t), intent(inout) :: this
       real, intent(in) :: climate(:)
       real, intent(in) :: timestep_days
     end subroutine decomp_demand

     subroutine decomp_final(this, climate, timestep_days, n_alloc_soil)
       ! Evaluate the final decomposition given the N allocation.
       import :: base_yasso_t
       class(base_yasso_t), intent(inout) :: this
       real, intent(in) :: climate(:)
       real, intent(in) :: timestep_days
       real, intent(in) :: n_alloc_soil(:)
     end subroutine decomp_final

     subroutine eval_tend(this, climate, timestep_days)
       import :: base_yasso_t
       class(base_yasso_t), intent(inout) :: this
       real, intent(in) :: climate(:)
       real, intent(in) :: timestep_days
     end subroutine eval_tend

     subroutine update_state(this, cflux, nflux)
       ! Update the state vector according to the internal tendencies (calculated earlier)
       ! and forced tendencies (cflux and nflux).
       import :: base_yasso_t
       class(base_yasso_t), intent(inout) :: this
       real, intent(in) :: cflux(:), nflux(:) ! C and N flux into the model. These shall
                                              ! have the dimensions of csize and nsize.

     end subroutine update_state

     subroutine force_n(this, cn)
       ! force the N pool(s) given the current C pools, and a C:N ratio.
       import :: base_yasso_t
       class(base_yasso_t), intent(inout) :: this
       real, intent(in) :: cn
     end subroutine force_n

     subroutine cmatrix(this, clim, matr)
       ! if possible return matrix such that matr*cstate = ctend
       import :: base_yasso_t
       class(base_yasso_t), intent(in) :: this
       real, intent(in) :: clim(:)
       real, intent(out) :: matr(:,:)
     end subroutine cmatrix

     subroutine cstate(this, get, set)
       ! Set/get the vector (size given by get_csize()) of carbon pools.
       import :: base_yasso_t
       class(base_yasso_t) :: this
       real, intent(out), optional :: get(:)
       real, intent(in), optional :: set(:)
     end subroutine cstate

     subroutine nstate(this, get, set)
       ! Set/get the vector (size given by get_nsize()) of nitrogen pools.
       import :: base_yasso_t
       class(base_yasso_t) :: this
       real, intent(out), optional :: get(:)
       real, intent(in), optional :: set(:)
     end subroutine nstate

     subroutine load_state(this, filename, fileunit)
       ! Load the model state from a file as written by store_state.
       import :: base_yasso_t
       class(base_yasso_t), intent(inout) :: this
       character(len=*), intent(in) :: filename
       integer, intent(in), optional :: fileunit
     end subroutine load_state

     subroutine init(this, params, aux_int, aux_real)
       ! Initialize a generic yasso-N instance. The optional arguments allow passing
       ! additional integer or real parameters as needed by subclasses.
       import :: base_yasso_t
       class(base_yasso_t), intent(out) :: this
       real, intent(in) :: params(:)
       integer, intent(in), optional :: aux_int(:)
       real, intent(in), optional :: aux_real(:)
       ! No point in creating a default implementation for base_yasso_t since it seems
       ! that subclasses cannot call the methods of an abstract super.
     end subroutine init
     
     function litt_awenh_to_cflux(this, litt_awenh) result(cflux)
       ! Attribute the litter C input AWENH into the model pools as reqruied for
       ! update_state.
       import :: base_yasso_t
       class(base_yasso_t) :: this
       real, intent(in) :: litt_awenh(:)
       real :: cflux(this%get_csize())
     end function litt_awenh_to_cflux

     function litt_n_to_nflux(this, litt_n) result(nflux)
       ! Attribute the litter C input AWENH into the model pools as reqruied for
       ! update_state.
       import :: base_yasso_t
       class(base_yasso_t) :: this
       real, intent(in) :: litt_n
       real :: nflux(this%get_nsize())
     end function litt_n_to_nflux

  end interface

contains

  !**********************************************************************************
  !
  ! Non-abstract methods for base_yasso_t.
  !
  !**********************************************************************************
  
  pure integer function get_demand_size(this) result(s)
    ! Return the number of N demand components, as use in get_n_demand().
    class(base_yasso_t), intent(in) :: this
    s = this%demand_size
  end function get_demand_size

  pure integer function get_csize(this) result(s)
    ! Return the number of C pools, as required for the cstate subroutine and cflux
    ! argument in update_state.
    class(base_yasso_t), intent(in) :: this
    s = this%cstate_size
  end function get_csize

  pure integer function get_nsize(this) result(s)
    ! Return the number of N pools, as required for the cstate subroutine and cflux
    ! argument in update_state.
    class(base_yasso_t), intent(in) :: this
    s = this%nstate_size
  end function get_nsize

  subroutine delete(this)
    class(base_yasso_t), intent(in) :: this
    ! Subclasses that allocate memory should override this and connect it to the
    ! destructor.
  end subroutine delete

  logical function test(this) result(success)
    class(base_yasso_t), intent(in) :: this
    ! Subclasses may use this to run type-specific test code.
    success = .true.
  end function test
  
  subroutine init_from_files(this, filename_param, filename_initcn, fileunit)
    ! Initialize this instance from parameter and initial condition files.
    class(base_yasso_t), intent(out) :: this
    character(len=*), intent(in) :: filename_param
    character(len=*), intent(in) :: filename_initcn
    integer, intent(in), optional :: fileunit
    
    real :: param(num_parameters)
    integer :: unit

    if (present(fileunit)) then
       unit = fileunit
    else
       unit = fileunit_def
    end if
    call load_params(filename_param, param, unit)
    call this%init(param)
    call this%load_state(filename_initcn, unit)
    
  end subroutine init_from_files
  
  subroutine set_sizes(this, demand_size, csize, nsize)
    class(base_yasso_t), intent(inout) :: this
    integer, intent(in) :: demand_size, csize, nsize

    this%demand_size = demand_size
    this%cstate_size = csize
    this%nstate_size = nsize
    
  end subroutine set_sizes

  subroutine get_params(this, params)
    class(base_yasso_t), intent(in) :: this
    real, intent(out) :: params(num_parameters)

    params = this%param
    
  end subroutine get_params

  subroutine set_params(this, params)
    class(base_yasso_t), intent(inout) :: this
    real, intent(in) :: params(num_parameters)

    this%param = params
    
  end subroutine set_params


  !**********************************************************************************
  !
  ! Common subroutines for the YASSO models
  !
  !**********************************************************************************
  
  subroutine load_yasso_clim_1d(filename, clim, fileunit)
    ! Load the constant climate forcing, one line in file.
    character(len=*), intent(in) :: filename
    real, intent(out) :: clim(:) ! climpar
    integer, intent(in), optional :: fileunit
    
    integer :: unit
    if (present(fileunit)) then
       unit = fileunit
    else
       unit = fileunit_def
    end if
    open(unit, file=filename, action='read', status='old')
    read(unit, *) clim(1:num_clim_par)
    close(unit)
    
  end subroutine load_yasso_clim_1d

  subroutine load_yasso_clim_2d(filename, clim, fileunit)
    ! Load the daily climate forcing, multi-line file.
    use iso_fortran_env, only : iostat_end
    character(len=*), intent(in) :: filename
    real, intent(out) :: clim(:,:) ! climpar, day
    integer, intent(in), optional :: fileunit
    
    integer :: unit, iostat, ind_row
    real :: tmp(num_clim_par+2)
    
    if (present(fileunit)) then
       unit = fileunit
    else
       unit = fileunit_def
    end if
    open(unit, file=filename, action='read', status='old')
    ind_row = 0
    do
       read(unit, *, iostat=iostat) tmp
       if (iostat == 0) then
          ind_row = ind_row + 1
          if (ind_row > size(clim,2)) then
             print *, 'clim argument too small'
             stop
          end if
          clim(:, ind_row) = tmp
       else if (iostat == iostat_end) then
          exit
       else
          print *, 'Bad iostat:', iostat
          stop
       end if
    end do
    close(unit)
    clim(:,ind_row+1:) = 0.0
    
  end subroutine load_yasso_clim_2d
  
  subroutine get_daily_clim(year, doy, clim_full, clim_day, found)
    integer, intent(in) :: year, doy
    real, intent(in) :: clim_full(:,:)
    real, intent(out) :: clim_day(:)
    logical, intent(out) :: found
    
    integer :: ind_last = 1
    integer :: ind_row
    
    if (clim_full(1,ind_last) > year .or. (abs(clim_full(1,ind_last)-year) < 1e-5 .and. clim_full(2,ind_last) > doy)) then
       ! start over
       ind_last = 1
    end if
    found = .false.
    do ind_row = ind_last, size(clim_full, 2)
       if (abs(clim_full(1,ind_row)-year) < 1e-5 .and. abs(clim_full(2,ind_row)-doy) < 1e-5) then
          clim_day = clim_full(3:, ind_row)
          ind_last = ind_row
          found = .true.
          exit
       end if
    end do
    
  end subroutine get_daily_clim
  
  function get_littc_awenh(leaf_littc, fineroot_littc) result(litter_awenh)
    real, intent(in) :: leaf_littc, fineroot_littc
    real :: litter_awenh(5)
    litter_awenh = leaf_littc*litterfract_leaf + fineroot_littc*litterfract_fineroot
  end function get_littc_awenh

  function get_fertc_awenh(typeflag, fertc) result(fertc_awenh)
    real, intent(in) :: typeflag
    real, intent(in) :: fertc

    real :: fertc_awenh(5)

    if (typeflag == 0.0) then
       ! mineral
       fertc_awenh = 0.0
    else if (0.9 <  typeflag .and. typeflag < 1.1) then
       fertc_awenh = (/0.0, fertc, 0.0, 0.0, 0.0/)
    else
       print *, 'Bad typeflag to get_fertc_awenh'
       stop
    end if
    
  end function get_fertc_awenh

  subroutine load_params(filename_params, params, fileunit)
    character(len=*), intent(in) :: filename_params
    real, intent(out) :: params(:)
    integer, intent(in), optional :: fileunit

    integer :: fun

    if (present(fileunit)) then
       fun = fileunit
    else
       fun = fileunit_def
    end if
    open(fun, file=filename_params, action='read', status='old')
    read(fun, *) params(1:num_parameters)
    close(fun)
    
  end subroutine load_params
    
  subroutine resolve_ndemand(method, nmin_avail, &
       demand_soil, demand_plant, alloc_soil, alloc_plant)
    character(len=*), intent(in) :: method
    real, intent(in) :: nmin_avail
    real, intent(in) :: demand_soil(:)
    real, intent(in) :: demand_plant(:)
    real, intent(out) :: alloc_soil(:)
    real, intent(out) :: alloc_plant(:)

    real :: total_demand
    real :: inhib
    
    if (method == 'demand_based') then
       ! distribute available N proportional to demand. AWEN decomposition is assumed to
       ! acquire N from the H pool only in mineral form (note difference from get_actfluxes).
       total_demand = sum(demand_soil, mask=demand_soil > 0) + sum(demand_plant, mask=demand_plant > 0)
       if (total_demand < nmin_avail) then
          alloc_soil = demand_soil
          alloc_plant = demand_plant
          inhib = 1.0
       else
          inhib = nmin_avail / total_demand
          where (demand_soil > 0)
             alloc_soil = demand_soil*inhib
          elsewhere
             alloc_soil = demand_soil
          end where
          where (demand_plant > 0)
             alloc_plant = demand_plant*inhib
          elsewhere
             alloc_plant = demand_plant
          end where
       end if
    else
       print *, 'Bad method'
       stop
    end if

!!$    print *, 'Resolve demand:'
!!$    print *, 'N available:', nmin_avail
!!$    print *, 'Soil demand:', demand_soil
!!$    print *, 'Plant demand:', demand_plant
!!$    print *, 'Soil alloc:', alloc_soil
!!$    print *, 'Plant alloc:', alloc_plant
!!$    print *, 'Inhib:', inhib
    
  end subroutine resolve_ndemand

  
  
end module yassobase
