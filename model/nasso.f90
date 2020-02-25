module nasso
  implicit none

  public

  integer, parameter, public :: state_size = 6
  integer, parameter, public :: ind_norg = 6, ind_nmin = 7
  integer, parameter, public :: num_c_pools = 5

  integer, parameter, public :: num_fluxes = 3
  integer, parameter, public :: ind_dn = 1, ind_resp = 2, ind_decomp = 3
  
  integer, parameter, public :: num_parameters = 38
  integer, parameter, public :: ind_nc_humus = 36, ind_nc_microb = 37, ind_eff = 38
  
  integer, parameter, public :: num_clim_par = 3
#ifndef F2PY
  integer, parameter, public :: rp = 8
#endif
  real(rp), parameter, public :: pi = 3.141592653589793_rp

contains

  integer function get_ind_norg()
    get_ind_norg = ind_norg
  end function get_ind_norg

  integer function get_ind_nmin()
    get_ind_nmin = ind_nmin
  end function get_ind_nmin

  integer function get_ind_nc_humus()
    get_ind_nc_humus = ind_nc_humus
  end function get_ind_nc_humus

  integer function get_ind_nc_microb()
    get_ind_nc_microb = ind_nc_microb
  end function get_ind_nc_microb

  integer function get_ind_eff()
    get_ind_eff = ind_eff
  end function get_ind_eff

  subroutine get_matrix(theta, climate, leach, diam, A)
    implicit none
    !f2py integer, intent(aux) :: num_clim_par
    REAL(rp),DIMENSION(35),INTENT(IN) :: theta ! parameters /J as listed above/
    ! climatic conditions, /J 1 temp, 2 prec, 3 temperature amplitude/
    REAL(rp),DIMENSION(num_clim_par),INTENT(IN) :: climate 
    real(rp), intent(in) :: leach, diam
    real(rp), dimension(5,5), intent(out) :: A

    real :: te(4), tem, temN, temH, size_dep
    integer :: i

    !#########################################################################
    ! Compute the coefficient matrix A for the differential equation

    ! temperature annual cycle approximation
    te(1) = climate(1)+4*climate(3)*(1/SQRT(2.0)-1)/pi
    te(2) = climate(1)-4*climate(3)/SQRT(2.0)/pi
    te(3) = climate(1)+4*climate(3)*(1-1/SQRT(2.0))/pi
    te(4) = climate(1)+4*climate(3)/SQRT(2.0)/pi

    tem = 0.0
    temN = 0.0
    temH = 0.0
    DO i = 1,4 ! Average temperature dependence
       tem = tem+EXP(theta(22)*te(i)+theta(23)*te(i)**2.0)/4.0 ! Gaussian
       temN = temN+EXP(theta(24)*te(i)+theta(25)*te(i)**2.0)/4.0
       temH = temH+EXP(theta(26)*te(i)+theta(27)*te(i)**2.0)/4.0
    END DO

    ! Precipitation dependence
    tem = tem*(1.0-EXP(theta(28)*climate(2)/1000.0))
    temN = temN*(1.0-EXP(theta(29)*climate(2)/1000.0))
    temH = temH*(1.0-EXP(theta(30)*climate(2)/1000.0))

    ! Size class dependence -- no effect if d == 0.0
    size_dep = MIN(1.0,(1.0+theta(33)*diam+theta(34)*diam**2.0)**(-ABS(theta(35))))

    ! Calculating matrix A (will work ok despite the sign of alphas)
    DO i = 1,3
       A(i,i) = -ABS(theta(i))*tem*size_dep
    END DO
    A(4,4) = -ABS(theta(4))*temN*size_dep

    A(1,2) = theta(5)*ABS(A(2,2))
    A(1,3) = theta(6)*ABS(A(3,3))
    A(1,4) = theta(7)*ABS(A(4,4))
    A(1,5) = 0.0 ! no mass flows from H -> AWEN
    A(2,1) = theta(8)*ABS(A(1,1))
    A(2,3) = theta(9)*ABS(A(3,3))
    A(2,4) = theta(10)*ABS(A(4,4))
    A(2,5) = 0.0
    A(3,1) = theta(11)*ABS(A(1,1))
    A(3,2) = theta(12)*ABS(A(2,2))
    A(3,4) = theta(13)*ABS(A(4,4))
    A(3,5) = 0.0
    A(4,1) = theta(14)*ABS(A(1,1))
    A(4,2) = theta(15)*ABS(A(2,2))
    A(4,3) = theta(16)*ABS(A(3,3))
    A(4,5) = 0.0
    A(5,5) = -ABS(theta(32))*temH ! no size effect in humus
    DO i = 1,4
       A(5,i) = theta(31)*ABS(A(i,i)) ! mass flows AWEN -> H (size effect is present here)
    END DO

    ! Leaching (no leaching for humus) /J this is where the correction term is considered/
    DO i = 1,4
       A(i,i) = A(i,i)+leach*climate(2) * 1e-3
    END DO

  end subroutine get_matrix

  subroutine get_potfluxes(matr, state, timestep, fluxes_awen, fluxes_h, nc_humus, nc_microb, eff_awen, eff_hum)
    implicit none
    !f2py integer, intent(aux) :: num_fluxes, num_c_pools, state_size
    real(rp), intent(in), dimension(num_c_pools,num_c_pools) :: matr
    real(rp), intent(in), dimension(state_size) :: state
    real(rp), intent(out), dimension(num_fluxes) :: fluxes_awen, fluxes_h
    real(rp), intent(in) :: timestep, nc_humus, eff_awen, eff_hum, nc_microb

    real(rp) :: awen, n_hum, nc_awen
    
    n_hum = state(5) * nc_humus
    if (sum(state(1:4)) <= maxval(state)*1e-5) then
       fluxes_awen = 0.0
    else       
       nc_awen = (state(ind_norg) - n_hum) / sum(state(1:4))
       fluxes_awen(ind_resp) = timestep * sum(matmul(matr(1:5,1:4), state(1:4)))
       fluxes_awen(ind_decomp) = fluxes_awen(ind_resp) /  (1.0 - eff_awen)
       fluxes_awen(ind_dn) = (nc_awen - nc_microb*eff_awen) * fluxes_awen(ind_decomp)
    end if
    fluxes_h(ind_resp) = timestep * matr(5,5) * state(5)
    fluxes_h(ind_decomp) = fluxes_h(ind_resp) / (1.0 - eff_hum)
    fluxes_h(ind_dn) = (nc_humus - nc_microb*eff_hum) * fluxes_h(ind_decomp)

    !print *, 'ncs', nc_awen, nc_microb, state(ind_norg),  n_hum
    
  end subroutine get_potfluxes

  subroutine get_actfluxes(state, fluxes_awen_pot, fluxes_h, matr, fluxes_awen_act, inhib, matr_inhib)
    implicit none
    !f2py integer, intent(aux) :: num_fluxes, num_c_pools, state_size
    real(rp), intent(in), dimension(state_size+1) :: state ! with nmin
    real(rp), intent(in), dimension(num_fluxes) :: fluxes_awen_pot, fluxes_h
    real(rp), intent(in), dimension(num_c_pools, num_c_pools) :: matr
    real(rp), intent(out), dimension(num_fluxes) :: fluxes_awen_act
    real(rp), intent(out) :: inhib
    real(rp), intent(out), dimension(num_c_pools, num_c_pools) :: matr_inhib

    if (fluxes_awen_pot(ind_dn) + fluxes_h(ind_dn) > state(ind_nmin)) then
       inhib = (state(ind_nmin) - fluxes_h(ind_dn)) / fluxes_awen_pot(ind_dn)
    else
       inhib = 1.0
    end if

    fluxes_awen_act = inhib*fluxes_awen_pot

    matr_inhib = matr * inhib
    matr_inhib(5,5) = matr(5,5)
    
  end subroutine get_actfluxes

  subroutine nmodel2(param, timestep, climate, littsize, leach, littc, littn, nuptake, state, state_new)
    implicit none
    !f2py integer, intent(aux) :: num_parameters, num_clim_par, num_c_pools, state_size
    !
    ! param: first 35 parameters as in mod5c + nc_humus, nc_microb, eff
    real(rp), intent(in), dimension(num_parameters) :: param
    real(rp), intent(in) :: timestep ! years
    real(rp), intent(in), dimension(num_clim_par) :: climate ! temp, prec, tempr amplitude
    real(rp), intent(in) :: littsize ! litter size
    real(rp), intent(in) :: leach ! leaching parameter
    real(rp), intent(in), dimension(num_c_pools) :: littc ! litter C input per timestep unit
    real(rp), intent(in) :: littn ! litter N input per timestep unit
    real(rp), intent(in) :: nuptake ! mineral N uptake per timestep unit
    real(rp), intent(in), dimension(state_size+1) :: state ! AWENH, NORG, NMIN
    real(rp), intent(out), dimension(size(state)) :: state_new 

    real(rp), dimension(num_c_pools, num_c_pools) :: matr, matr_inhib
    real(rp), dimension(num_fluxes) :: fluxes_awen_pot, fluxes_awen_act, fluxes_h
    !real(rp) :: state_new(state_size)
    real(rp) :: inhib, dn_tot

    call get_matrix(param(1:35), climate, leach, littsize, matr)
    call get_potfluxes(matr, state, timestep, fluxes_awen_pot, fluxes_h, &
         param(ind_nc_humus), param(ind_nc_microb), param(ind_eff), param(ind_eff))
    call get_actfluxes(state, fluxes_awen_pot, fluxes_h, matr, fluxes_awen_act, inhib, matr_inhib)
    dn_tot = fluxes_awen_act(ind_dn) + fluxes_h(ind_dn)
    
    state_new(1:num_c_pools) = state(1:num_c_pools) &
         + timestep * matmul(matr_inhib, state(1:num_c_pools))
    state_new(ind_norg) = state(ind_norg) + dn_tot + littn*timestep
    state_new(ind_nmin) = state(ind_nmin) - dn_tot - nuptake*timestep
    !print *, 'before', sum(state_new(1:num_c_pools)), timestep
    state_new(1:num_c_pools) = state_new(1:num_c_pools) + littc*timestep
    state_new(ind_nmin) = max(state_new(ind_nmin), 0.0_rp)
    !print *, 'after', sum(state_new(1:num_c_pools))

    !state_new(ind_nmin) = state_new(ind_nmin) - nuptake*timestep
    !state_new(ind_norg) = state_new(ind_norg) + littn*timestep
    !print *, state(ind_norg)+state(ind_nmin) - state_new(ind_norg) - state_new(ind_nmin)
    !print *, state_new
  end subroutine nmodel2


end module nasso
