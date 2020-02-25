module cascademodel
  use nasso
  use yassobase
  implicit none
  
  private
  public casc_from_files
  public get_age_tend
  
  integer, parameter :: fileunit_def = 30
    
  public cascademodel_t
  type, extends(base_yasso_t) :: cascademodel_t
     !private
     real, allocatable :: state(:,:) ! (awenhN, age)
     real, allocatable :: ctend(:,:)
     real, allocatable :: ntend(:)

     real :: matr(num_c_pools, num_c_pools)
     real :: clim(num_clim_par)

     real, allocatable :: fluxes_awen_pot(:,:) ! (flux, age)
     real, allocatable :: fluxes_h(:,:)
     real, allocatable :: fluxes_awen_act(:,:)

     real, allocatable :: ranges(:) ! in days (num_cls)
     
     integer :: num_cls
          
     real, allocatable :: inhib(:)

   contains
     
     procedure, pass(this) :: get_n_demand
     procedure, pass(this) :: get_norg
     procedure, pass(this) :: get_netmin_act
     procedure, pass(this) :: get_totc
     procedure, pass(this) :: get_cn_awen
     procedure, pass(this) :: get_soilresp
     procedure, pass(this) :: map_to_basgra
     procedure, pass(this) :: store_state
     procedure, pass(this) :: decomp_demand
     procedure, pass(this) :: decomp_final
     procedure, pass(this) :: update_state
     procedure, pass(this) :: force_n
     procedure, pass(this) :: cstate
     procedure, pass(this) :: nstate
     procedure, pass(this) :: cmatrix
     procedure, pass(this) :: litt_awenh_to_cflux
     procedure, pass(this) :: litt_n_to_nflux
          
     procedure, pass(this), private :: setup
     procedure, pass(this), private :: eval_tend
     
     procedure, pass(this) :: init
     procedure, pass(this) :: load_state
     procedure, pass(this) :: cstate_cls
     procedure, pass(this) :: get_ranges
     procedure, pass(this) :: report
     
     final :: delete
  end type cascademodel_t

  public nmodel1cls
  type, extends(cascademodel_t) :: nmodel1cls
   contains
     procedure, pass(this) :: init => init_1cls
  end type nmodel1cls

  public nmodel2cls
  type, extends(cascademodel_t) :: nmodel2cls
   contains
     procedure, pass(this) :: init => init_2cls
  end type nmodel2cls
  
contains
  
  subroutine delete(this)
    type(cascademodel_t), intent(inout) :: this
    if (allocated(this%state)) then
       deallocate(this%state, this%ctend, this%ntend, this%fluxes_awen_pot, this%fluxes_h, this%fluxes_awen_act, this%ranges)
    end if
  end subroutine delete
  
  subroutine init(this, params, aux_int, aux_real)
    class(cascademodel_t), intent(out) :: this
    real, intent(in) :: params(:)
    integer, intent(in), optional :: aux_int(:)
    real, intent(in), optional :: aux_real(:)

    if (.not. present(aux_int) .or. .not. present(aux_real)) then
       print *, 'ERROR: aux_int(1) = num_cls and aux_real(1:num_cls) = ranges must be given for cascademodel_t'
       stop
    end if
    if (size(aux_real) < aux_int(1)) then
       print *, 'ERROR: aux_real too small'
       stop
    end if
    if (aux_int(1) > 1) then
       if (.not. all(aux_real(1:aux_int(1)-1) > 0)) then
          print *, 'ERROR: not all ranges(1:num_cls-1) positive'
          stop
       end if
    end if
    call this%set_params(params)
    call this%setup(aux_int(1), aux_real(1:aux_int(1)))
    
  end subroutine init

 function get_n_demand(this) result(demand)
    class(cascademodel_t), intent(in) :: this
    real :: demand(this%demand_size)

    integer :: ind_cls, ind_demand

    ind_demand = 1
    do ind_cls = 1, this%num_cls
       demand(ind_demand) = this%fluxes_awen_pot(ind_dn,ind_cls)
       demand(ind_demand+1) = this%fluxes_h(ind_dn,ind_cls)
       ind_demand = ind_demand + 2
    end do
    
  end function get_n_demand

  real function get_norg(this) result(norg)
    class(cascademodel_t), intent(in) :: this
    norg = sum(this%state(ind_norg,:))
  end function get_norg

  real function get_netmin_act(this) result(netmin)
    class(cascademodel_t), intent(in) :: this
    netmin = - (sum(this%fluxes_awen_act(ind_dn,:)) + sum(this%fluxes_h(ind_dn,:)))
  end function get_netmin_act

  real function get_totc(this) result(totc)
    class(cascademodel_t), intent(in) :: this
    totc = sum(this%state(1:num_c_pools,:))
  end function get_totc

  
  function get_cn_awen(this) result(ratio)
    class(cascademodel_t), intent(in) :: this
    real :: ratio
    
    real :: n_hum
    
    n_hum = sum(this%state(5,:)) * this%param(ind_nc_humus)
    ratio = sum(this%state(1:4,:)) / (get_norg(this) - n_hum)
  end function get_cn_awen

  function get_soilresp(this) result(resp)
    class(cascademodel_t), intent(in) :: this
    real :: resp
    resp = -sum(this%ctend)
  end function get_soilresp

  subroutine map_to_basgra(this, clitt, csomf, csoms, nlitt, nsomf, nsoms)
    class(cascademodel_t), intent(in) :: this
    real, intent(out) :: clitt, csomf, csoms, nlitt, nsomf, nsoms

    clitt = sum(this%state(1:3,:))
    csomf = sum(this%state(4,:))
    csoms = sum(this%state(5,:))

    nlitt = clitt / this%get_cn_awen()
    nsomf = csomf / this%get_cn_awen()
    nsoms = csoms * this%param(ind_nc_humus)
    
  end subroutine map_to_basgra

  subroutine setup(this, num_cls, ranges)
    class(cascademodel_t), intent(inout) :: this
    integer, intent(in) :: num_cls
    real, intent(in) :: ranges(:)

    this%num_cls = num_cls
    call this%set_sizes(demand_size=2*this%num_cls, csize=num_c_pools*this%num_cls, nsize=this%num_cls)
    allocate(this%state(6,num_cls), this%ctend(5,num_cls), this%ntend(num_cls), &
         this%fluxes_awen_pot(3,num_cls), this%fluxes_h(3,num_cls), this%fluxes_awen_act(3,num_cls), &
         this%inhib(num_cls), this%ranges(num_cls))

    this%ranges = ranges
    
    this%state = 0.0
    this%ctend = 0.0
    this%ntend = 0.0
    this%fluxes_awen_pot = 0.0
    this%fluxes_awen_act = 0.0
    this%fluxes_h = 0.0
    this%inhib = 1.0
    
  end subroutine setup

  ! Note storage order different from memory!!
  
  subroutine store_state(this, filename, fileunit)
    class(cascademodel_t), intent(in) :: this
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: fileunit
    
    real :: cstate(this%num_cls*num_c_pools), nstate(this%num_cls)
    integer :: unit

    if (present(fileunit)) then
       unit = fileunit
    else
       unit = fileunit_def
    end if
    
    open(unit, file=filename, action='write')
    call this%cstate(get=cstate)
    call this%nstate(get=nstate)
    write(unit, *) cstate, nstate
    close(unit) 

  end subroutine store_state

  subroutine load_state(this, filename, fileunit)
    class(cascademodel_t), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: fileunit

    integer :: unit
    real :: cstate(this%num_cls*num_c_pools), nstate(this%num_cls)
    
    if (present(fileunit)) then
       unit = fileunit
    else
       unit = fileunit_def
    end if
      
    open(unit, file=filename, action='read')
    read(unit, *) cstate, nstate
    close(unit)

    call this%cstate(set=cstate)
    call this%nstate(set=nstate)

    if (any(this%state < 0)) then
       print *, 'Negative values in state from file:', filename
       stop
    end if

    print *, this%state
  end subroutine load_state
  
  real function get_eff_awen(this, ind_cls) result(eff)
    class(cascademodel_t), intent(in) :: this
    integer, intent(in) :: ind_cls
    
    real :: cn_microb, ratio, n_hum

    if (sum(this%state(1:num_c_pools, ind_cls)) < 1e-15) then
       eff = 0.45
       return
    end if
    
    cn_microb = 1.0 / this%param(ind_nc_microb)

    n_hum = this%state(5,ind_cls) * this%param(ind_nc_humus)
    ratio = sum(this%state(1:4,ind_cls)) / (this%state(6,ind_cls) - n_hum)
    
    eff = 0.43 * (cn_microb / ratio) ** 0.6
    eff = min(eff, 0.6)
    eff = max(0.1, eff)
        
  end function get_eff_awen
  
  subroutine decomp_demand(this, climate, timestep_days)
    class(cascademodel_t), intent(inout) :: this
    real, intent(in) :: climate(:)
    real, intent(in) :: timestep_days

    real :: dt_years
    integer :: ind_cls
    
    dt_years = timestep_days / 365.0
    call get_matrix(this%param, climate, 0.0, 0.0, this%matr)
    do ind_cls = 1, this%num_cls
       call get_potfluxes(this%matr, this%state(:,ind_cls), dt_years, &
            this%fluxes_awen_pot(:,ind_cls), this%fluxes_h(:,ind_cls), &
            this%param(ind_nc_humus), this%param(ind_nc_microb), get_eff_awen(this,ind_cls), 0.45)
    end do

  end subroutine decomp_demand

  subroutine decomp_final(this, climate, timestep_days, n_alloc_soil)
    class(cascademodel_t), intent(inout) :: this
    real, intent(in) :: climate(:)
    real, intent(in) :: timestep_days
    real, intent(in) :: n_alloc_soil(:)

    real :: demand_awen
    integer :: ind_cls, ind1d
    
    do ind_cls = 1, this%num_cls
       demand_awen = this%fluxes_awen_pot(ind_dn,ind_cls)
       ind1d = (ind_cls-1)*2 + 1
       if (demand_awen > n_alloc_soil(1)) then
          this%inhib(ind_cls) = n_alloc_soil(ind1d) / demand_awen
          this%fluxes_awen_act(:,ind_cls) = this%inhib(ind_cls) * this%fluxes_awen_pot(:,ind_cls)
       else
          this%fluxes_awen_act(:,ind_cls) = this%fluxes_awen_pot(:,ind_cls)
          this%inhib(ind_cls) = 1.0
       end if
       !this%fluxes_h_act(:,ind_cls) = this%fluxes_h_pot(:,ind_cls) ! no inhibition for humus pool
    end do

    call this%eval_tend(climate, timestep_days)
  end subroutine decomp_final
  
  subroutine eval_tend(this, climate, timestep_days)
    class(cascademodel_t), intent(inout) :: this
    real, intent(in) :: climate(:)
    real, intent(in) :: timestep_days

    real, parameter :: littsize = 0.0, leach = 0.0
    
    real, dimension(num_c_pools, num_c_pools) :: matr_inhib
    real, dimension(this%num_cls) :: agetend_tmp
    
    real :: dt_years, netmin
    integer :: ind_cls, ind_c
    
    dt_years = timestep_days / 365.0

    ! Decomposition tendencies
    do ind_cls = 1, this%num_cls
       matr_inhib = this%matr * this%inhib(ind_cls)
       matr_inhib(5,5) = this%matr(5,5) ! no inhibition of humus decomp
       this%ctend(:,ind_cls) = matmul(matr_inhib, this%state(1:num_c_pools,ind_cls))*dt_years
       !print *, '***', ind_cls, this%ctend(:,ind_cls)
       netmin = -(this%fluxes_awen_act(ind_dn,ind_cls) + this%fluxes_h(ind_dn,ind_cls))
       this%ntend(ind_cls) = -netmin

    end do
    ! Aging

    do ind_c = 1, num_c_pools
       call get_age_tend(this%state(ind_c,:), this%ranges, timestep_days, agetend_tmp, this%num_cls)
       this%ctend(ind_c,:) = this%ctend(ind_c,:) + agetend_tmp
    end do
    call get_age_tend(this%state(ind_norg,:), this%ranges, timestep_days, agetend_tmp, this%num_cls)
    this%ntend = this%ntend + agetend_tmp
    
  end subroutine eval_tend

  subroutine update_state(this, cflux, nflux)
    class(cascademodel_t), intent(inout) :: this
    real, intent(in) :: cflux(:), nflux(:)

    integer :: ind_cls, ind_start, ind_end
    real :: cflux_cls(num_c_pools)

    do ind_cls = 1, this%num_cls
       ind_start = (ind_cls-1)*num_c_pools + 1
       ind_end = (ind_cls-1)*num_c_pools + num_c_pools
       this%state(1:num_c_pools,ind_cls) = this%state(1:num_c_pools,ind_cls) + this%ctend(:,ind_cls) + cflux(ind_start:ind_end)
       this%state(num_c_pools+1,ind_cls) = this%state(num_c_pools+1,ind_cls) + this%ntend(ind_cls) + nflux(ind_cls)
    end do
  end subroutine update_state

  subroutine force_n(this, cn)
    class(cascademodel_t), intent(inout) :: this
    real, intent(in) :: cn

    integer :: ind_cls
    real :: ccls
    
    do ind_cls = 1, this%num_cls
       ccls = sum(this%state(1:5,ind_cls))
       this%state(6,ind_cls) = ccls / cn
    end do
    
  end subroutine force_n

  subroutine casc_from_files(this, filename_param, filename_initcn, num_cls, incr_days)
    use nasso, only : num_parameters
    class(cascademodel_t), intent(out) :: this
    character(len=*), intent(in) :: filename_param
    character(len=*), intent(in) :: filename_initcn
    integer, intent(in) :: num_cls
    real, intent(in) :: incr_days

    real :: param(num_parameters), ranges(num_cls)
    integer :: ind_cls
    
    do ind_cls = 1, num_cls-1
       ranges(ind_cls) = incr_days
    end do
    ranges(num_cls) = 0.0
    
    call load_params(filename_param, param)
    call this%init(param, aux_int=(/num_cls/), aux_real=ranges)
    call this%load_state(filename_initcn)
  end subroutine casc_from_files

  subroutine cmatrix(this, clim, matr)
    class(cascademodel_t), intent(in) :: this
    real, intent(in) :: clim(:)
    real, intent(out) :: matr(:,:)

    integer :: indc, ind_cls, ind1d
    real :: ctmp(this%get_csize())
    class(cascademodel_t), allocatable :: this_copy
    real :: xmatr(num_c_pools, num_c_pools)
    
    allocate(this_copy)
    call this_copy%init(this%param, aux_int=(/this%num_cls/), aux_real=this%ranges)
    
    ! bruteforce
    do ind_cls = 1, this%num_cls
       do indc = 1, num_c_pools
          ind1d = (ind_cls-1)*num_c_pools + indc
          ctmp = 0.0
          ctmp(ind1d) = 1.0
          call this_copy%cstate(set=ctmp)
          call this_copy%decomp_demand(clim, timestep_days=1.0)
          call this_copy%decomp_final(clim, 1.0, this_copy%get_n_demand())
          matr(:,ind1d) = reshape(this_copy%ctend, (/num_c_pools*this%num_cls/))
       end do
    end do

    print *, 'matr ', matr(1,:)
    call get_matrix(this%param, clim, 0.0, 0.0, xmatr)
    print *, 'xmatr', xmatr(1,:) / 365.0
    !call get_matrix(this%param, clim, 0.0, 0.0, matr)
    
  end subroutine cmatrix
  
  subroutine cstate_cls(this, ind_cls, get, set, from)
    class(cascademodel_t) :: this
    integer, intent(in) :: ind_cls
    real, intent(out), optional :: get(:)
    real, intent(in), optional :: set(:)
    real, intent(inout), optional :: from(:,:)
    
    integer :: ind1d, ind_c

    if (present(from)) then
       if (present(get)) get(ind_c) = from(ind_c,ind_cls)
       if (present(set)) from(ind_c,ind_cls) = set(ind_c)
    else
       do ind_c = 1, num_c_pools
          if (present(get)) get(ind_c) = this%state(ind_c,ind_cls)
          if (present(set)) this%state(ind_c,ind_cls) = set(ind_c)
       end do
    end if
  end subroutine cstate_cls

  function get_ranges(this) result(ranges)
    class(cascademodel_t) :: this
    real :: ranges(this%num_cls)
    ranges = this%ranges
  end function get_ranges
  
  subroutine cstate(this, get, set)
    class(cascademodel_t) :: this
    real, intent(out), optional :: get(:)
    real, intent(in), optional :: set(:)
    
    integer :: ind1d, ind_c, ind_cls, ind_start, ind_end

    do ind_cls = 1, this%num_cls
       ind_start = (ind_cls-1)*num_c_pools + 1
       ind_end = (ind_cls-1)*num_c_pools + num_c_pools
       if (present(get)) get(ind_start:ind_end) = this%state(1:num_c_pools, ind_cls) !call this%cstate_cls(ind_cls, get=get(ind_start:ind_end))
       if (present(set)) this%state(1:num_c_pools, ind_cls) = set(ind_start:ind_end) !call this%cstate_cls(ind_cls, set=set(ind_start:ind_end))
       !do ind_c = 1, num_c_pools
       !   ind1d = (ind_cls-1)*num_c_pools + ind_c
       !   if (present(get)) get(ind1d) = this%state(ind_c,ind_cls)
       !   if (present(set)) this%state(ind_c,ind_cls) = set(ind1d)
       !end do
    end do
    
  end subroutine cstate

  
  subroutine nstate(this, get, set)
    class(cascademodel_t) :: this
    real, intent(out), optional :: get(:)
    real, intent(in), optional :: set(:)

    integer :: ind_cls

    do ind_cls = 1, this%num_cls
       if (present(get)) get(ind_cls) = this%state(num_c_pools+1,ind_cls)
       if (present(set)) this%state(num_c_pools+1,ind_cls) = set(ind_cls)
    end do

  end subroutine nstate
  
  subroutine init_1cls(this, params, aux_int, aux_real)
    class(nmodel1cls), intent(out) :: this
    real, intent(in) :: params(:)
    integer, intent(in), optional :: aux_int(:)
    real, intent(in), optional :: aux_real(:)

    if (present(aux_int)) then
       print *, 'aux_int must not be present for init_xcls'
       stop
    end if

    call this%cascademodel_t%init(params, aux_int=(/1/), aux_real=(/0.0/))
    
  end subroutine init_1cls

  subroutine init_2cls(this, params, aux_int, aux_real)
    class(nmodel2cls), intent(out) :: this
    real, intent(in) :: params(:)
    integer, intent(in), optional :: aux_int(:)
    real, intent(in), optional :: aux_real(:)

    if (present(aux_int)) then
       print *, 'aux_int must not be present for init_xcls'
       stop
    end if

    call this%cascademodel_t%init(params, aux_int=(/2/), aux_real=(/365.0, 0.0/))
    
  end subroutine init_2cls

  subroutine get_age_tend(values, ranges, timestep, tend, num_cls)
    real, intent(in) :: values(:)
    real, intent(in) :: ranges(:)
    real, intent(in):: timestep
    real, intent(out) :: tend(:)
    integer, intent(in) :: num_cls
    
    integer :: ind_cls

    do ind_cls = 1, num_cls - 1
       ! mass leaving
       tend(ind_cls) = -values(ind_cls) * timestep / ranges(ind_cls)
    end do
    tend(num_cls) = 0.0
    do ind_cls = 2, num_cls
       ! mass entering
       tend(ind_cls) = tend(ind_cls) + values(ind_cls-1) * timestep / ranges(ind_cls-1)
    end do
    
    ! no mass entering the first class, no mass leaving the last.
    
  end subroutine get_age_tend
  
 function litt_awenh_to_cflux(this, litt_awenh) result(cflux) 
    class(cascademodel_t) :: this
    real, intent(in) :: litt_awenh(:)
    real :: cflux(this%get_csize())

    cflux(1:num_c_pools) = litt_awenh
    if (this%num_cls > 1) cflux(num_c_pools+1:) = 0.0
    
  end function litt_awenh_to_cflux

  function litt_n_to_nflux(this, litt_n) result(nflux) 
    class(cascademodel_t) :: this
    real, intent(in) :: litt_n
    real :: nflux(this%get_nsize())

    nflux = 0.0
    nflux(1) = litt_n
    
  end function litt_n_to_nflux
  
  subroutine report(this)
    class(cascademodel_t), intent(in) :: this

    integer :: ind_cls
    real :: chum, nhum, cntot, cawen, nawen
    
    print *, '*** CASCADEMODEL ***'
    
    do ind_cls = 1, this%num_cls
       print *, 'CLASS', ind_cls
       print '(6A9)', 'A', 'W', 'E', 'N', 'H', 'NORG'
       print '(6E9.2)', this%state(:,ind_cls)
       print '(6A9)', 'CTOT', 'NTOT', 'CNTOT', 'CAWEN', 'NAWEN', 'CNAWEN'
       chum = this%state(5,ind_cls)
       nhum = chum*this%param(ind_nc_humus)
       cawen = sum(this%state(1:4,ind_cls))
       nawen = this%state(6,ind_cls) - nhum
       if (cawen > 0.0) then
          print '(2E9.2,F9.2,2E9.2,F9.2)', chum+cawen, this%state(6,ind_cls), (chum+cawen)/this%state(6,ind_cls), &
               cawen, nawen, cawen/nawen
       end if
    end do
    
  end subroutine report
  
end module cascademodel
