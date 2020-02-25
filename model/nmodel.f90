module nmodel
  use yassobase
  use nasso
  
  implicit none

  private

  public initcn2basgra_nmodel1
  
  public nmodel1
  type, extends(base_yasso_t) :: nmodel1
     
     real :: state(num_c_pools+1) ! the last index is a placeholder for nmin, but not used
     real :: ctend(num_c_pools)  ! per timestep
     real :: ntend               ! per timestep
     real :: matr(num_c_pools, num_c_pools)
     real :: clim(num_clim_par)
     real :: fluxes_awen_pot(3)
     real :: fluxes_h(3)
     real :: fluxes_awen_act(3)
     real :: inhib

   contains
     
     procedure, pass(this) :: get_n_demand
     procedure, pass(this) :: get_norg
     procedure, pass(this) :: get_netmin_act
     procedure, pass(this) :: get_totc
     procedure, pass(this) :: get_cn_awen
     procedure, pass(this) :: get_soilresp
     procedure, pass(this) :: map_to_basgra
     !procedure, pass(this) :: init_from_files
     procedure, pass(this) :: store_state
     procedure, pass(this) :: decomp_demand
     procedure, pass(this) :: decomp_final
     procedure, pass(this) :: update_state
     procedure, pass(this) :: force_n
     procedure, pass(this) :: cstate
     procedure, pass(this) :: nstate
     procedure, pass(this) :: cmatrix
     procedure, pass(this) :: load_state
     procedure, pass(this) :: init
     procedure, pass(this) :: litt_awenh_to_cflux
     procedure, pass(this) :: litt_n_to_nflux
     procedure, pass(this) :: report
     
     procedure, pass(this), private :: eval_tend
     
  end type nmodel1

  integer, public, save :: fileunit_def = 30
  integer, public, save :: days_per_year = 365
  
contains

   subroutine init(this, params, aux_int, aux_real)
     class(nmodel1), intent(out) :: this
     real, intent(in) :: params(:)
     integer, intent(in), optional :: aux_int(:)
     real, intent(in), optional :: aux_real(:)
     call this%set_params(params)
     call this%set_sizes(demand_size=2, csize=num_c_pools, nsize=1)
   end subroutine init
  
  function get_n_demand(this) result(demand)
    class(nmodel1), intent(in) :: this
    real :: demand(this%demand_size)
    demand(1) = this%fluxes_awen_pot(ind_dn)
    demand(2) = this%fluxes_h(ind_dn)
  end function get_n_demand
  
  real function get_norg(this) result(norg)
    class(nmodel1), intent(in) :: this
    norg = this%state(ind_norg)
  end function get_norg

  real function get_netmin_act(this) result(netmin)
    class(nmodel1), intent(in) :: this
    netmin = - (this%fluxes_awen_act(ind_dn) + this%fluxes_h(ind_dn))
  end function get_netmin_act

  real function get_totc(this) result(totc)
    class(nmodel1), intent(in) :: this
    totc = sum(this%state(1:num_c_pools))
  end function get_totc

  function get_cn_awen(this) result(ratio)
    class(nmodel1), intent(in) :: this
    real :: ratio
    
    real :: n_hum
    
    n_hum = this%state(5) * this%param(ind_nc_humus)
    ratio = sum(this%state(1:4)) / (get_norg(this) - n_hum)
  end function get_cn_awen

  function get_soilresp(this) result(resp)
    class(nmodel1), intent(in) :: this
    real :: resp
    resp = -sum(this%ctend)
  end function get_soilresp

   subroutine map_to_basgra(this, clitt, csomf, csoms, nlitt, nsomf, nsoms)
    class(nmodel1), intent(in) :: this
    real, intent(out) :: clitt, csomf, csoms, nlitt, nsomf, nsoms

    clitt = sum(this%state(1:3))
    csomf = this%state(4)
    csoms = this%state(5)

    nlitt = clitt / this%get_cn_awen()
    nsomf = csomf / this%get_cn_awen()
    nsoms = csoms * this%param(ind_nc_humus)
    
  end subroutine map_to_basgra

!!$  subroutine init_from_files(this, filename_param, filename_initcn, nmin)
!!$    class(nmodel1), intent(out) :: this
!!$    character(len=*), intent(in) :: filename_param
!!$    character(len=*), intent(in) :: filename_initcn
!!$    real, intent(out), optional :: nmin
!!$
!!$    real :: nmintmp, nhum
!!$    
!!$    this%demand_size = 2
!!$    this%cstate_size = num_c_pools
!!$    this%nstate_size = 1
!!$
!!$    open(fileunit, file=filename_param, action='read', status='old')
!!$    read(fileunit, *) this%param(1:num_parameters)
!!$    close(fileunit)
!!$
!!$    open(fileunit, file=filename_initcn, action='read', status='old')
!!$    read(fileunit, *) this%state(1:num_c_pools+1), nmintmp
!!$    if (present(nmin)) nmin = nmintmp
!!$    this%state(num_c_pools+2) = 0.0
!!$    
!!$    nhum = this%state(5) * this%param(ind_nc_humus)
!!$    if (nhum > this%state(6)) then
!!$       print *, 'Inconsistent N initialization'
!!$       stop
!!$    end if
!!$    
!!$  end subroutine init_from_files

  subroutine store_state(this, filename, fileunit)
    class(nmodel1), intent(in) :: this
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: fileunit

    integer :: unit
    
    if (present(fileunit)) then
       unit = fileunit
    else
       unit = fileunit_def
    end if
    open(unit, file=filename, action='write')
    write(unit, *) this%state
    close(unit)
    
  end subroutine store_state

  subroutine load_state(this, filename, fileunit)
    class(nmodel1), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: fileunit

    integer :: unit
    
    if (present(fileunit)) then
       unit = fileunit
    else
       unit = fileunit_def
    end if
    open(unit, file=filename, action='read')
    read(unit, *) this%state
    close(unit)
    
  end subroutine load_state

  real function get_eff_awen(yasso_inst) result(eff)
    type(nmodel1), intent(in) :: yasso_inst

    real :: cn_microb

    if (sum(yasso_inst%state(1:4)) < 1e-15) then
       eff = 0.45
       return
    end if
    
    cn_microb = 1.0 / yasso_inst%param(ind_nc_microb)
    
    eff = 0.43 * (cn_microb / get_cn_awen(yasso_inst)) ** 0.6
    eff = min(eff, 0.6)
    eff = max(0.1, eff)
        
  end function get_eff_awen

  
  subroutine decomp_demand(this, climate, timestep_days)
    class(nmodel1), intent(inout) :: this
    real, intent(in) :: climate(:)
    real, intent(in) :: timestep_days

    real :: dt_years

    dt_years = timestep_days / real(days_per_year)
    call get_matrix(this%param, climate, 0.0, 0.0, this%matr)
    call get_potfluxes(this%matr, this%state, dt_years, &
         this%fluxes_awen_pot, this%fluxes_h, &
         this%param(ind_nc_humus), this%param(ind_nc_microb), get_eff_awen(this), 0.45)
  end subroutine decomp_demand

  subroutine decomp_final(this, climate, timestep_days, n_alloc_soil)
    class(nmodel1), intent(inout) :: this
    real, intent(in) :: climate(:)
    real, intent(in) :: timestep_days
    real, intent(in) :: n_alloc_soil(:)

    real :: n_demand(2)

    n_demand = get_n_demand(this)
    if (n_demand(1) > n_alloc_soil(1)) then
       this%inhib = n_alloc_soil(1) / n_demand(1)
    else
       this%inhib = 1.0
    end if
    ! no inhibition for humus pool

    this%fluxes_awen_act = this%inhib * this%fluxes_awen_pot
    call this%eval_tend(climate, timestep_days)
    
  end subroutine decomp_final

  subroutine eval_tend(this, climate, timestep_days)
    real, intent(in) :: climate(:)
    real, intent(in) :: timestep_days
    class(nmodel1), intent(inout) :: this

    real, parameter :: littsize = 0.0, leach = 0.0
    
    real, dimension(num_c_pools, num_c_pools) :: matr_inhib
    real :: dt_years

    dt_years = timestep_days / 365.0

    matr_inhib = this%matr * this%inhib
    matr_inhib(5,5) = this%matr(5,5)
    this%ctend = matmul(matr_inhib, this%state(1:num_c_pools))*dt_years
    this%ntend = -this%get_netmin_act()
    
  end subroutine eval_tend
  
  subroutine update_state(this, cflux, nflux)
    class(nmodel1), intent(inout) :: this
    real, intent(in) :: cflux(:), nflux(:)

    this%state(1:num_c_pools) = this%state(1:num_c_pools) + this%ctend + cflux
    this%state(num_c_pools+1) = this%state(num_c_pools+1) + this%ntend + nflux(1)
    
  end subroutine update_state

  subroutine force_n(this, cn)
    class(nmodel1), intent(inout) :: this
    real, intent(in) :: cn
    this%state(ind_norg) = this%get_totc() / cn
  end subroutine force_n

   subroutine initcn2basgra_nmodel1(filename_param, filename_initcn, clitt, csomf, csoms, nlitt, nsomf, nsoms, nmin)
    character(len=*), intent(in) :: filename_initcn, filename_param
    real, intent(out) :: clitt, csomf, csoms, nlitt, nsomf, nsoms
    real, intent(out), optional :: nmin
    
    type(nmodel1) :: yasso_inst
    real :: nminloc
    
    call init_from_files(yasso_inst, filename_param, filename_initcn)
    call yasso_inst%map_to_basgra(clitt, csomf, csoms, nlitt, nsomf, nsoms)
    if (present(nmin)) nmin = nminloc
    
  end subroutine initcn2basgra_nmodel1

  subroutine cmatrix(this, clim, matr)
    class(nmodel1), intent(in) :: this
    real, intent(in) :: clim(:)
    real, intent(out) :: matr(:,:)

    call get_matrix(this%param, clim, 0.0, 0.0, matr)
    matr = matr/real(days_per_year)
    
  end subroutine cmatrix
    
  subroutine cstate(this, get, set)
    class(nmodel1) :: this
    real, intent(out), optional :: get(:)
    real, intent(in), optional :: set(:)

    if (present(get)) then
       get(1:num_c_pools) = this%state(1:num_c_pools)
    end if
    if (present(set)) then
       this%state(1:num_c_pools) = set(1:num_c_pools)
    end if
    
  end subroutine cstate

  subroutine nstate(this, get, set)
    class(nmodel1) :: this
    real, intent(out), optional :: get(:)
    real, intent(in), optional :: set(:)

    if (present(get)) then
       get(1) = this%state(num_c_pools+1)
    end if
    if (present(set)) then
       this%state(num_c_pools+1) = set(1)
    end if

  end subroutine nstate

  function litt_awenh_to_cflux(this, litt_awenh) result(cflux) 
    class(nmodel1) :: this
    real, intent(in) :: litt_awenh(:)
    real :: cflux(this%get_csize())

    cflux = litt_awenh
    
  end function litt_awenh_to_cflux

  function litt_n_to_nflux(this, litt_n) result(nflux) 
    class(nmodel1) :: this
    real, intent(in) :: litt_n
    real :: nflux(this%get_nsize())

    nflux(1) = litt_n
    
  end function litt_n_to_nflux
    
  subroutine report(this)
    class(nmodel1), intent(in) :: this
    print *, '*** NMODEL1 ***'
    print '(6A9)', 'A', 'W', 'E', 'N', 'H', 'NORG'
    print '(6E9.2)', this%state
    print *
  end subroutine report
  
  
end module nmodel
