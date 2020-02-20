program yassotest
  use yassobase
  use nmodel
  use cascademodel
  use nasso, only : num_clim_par, num_parameters

  implicit none

  abstract interface
     subroutine test(inst)
       import :: base_yasso_t
       class(base_yasso_t) :: inst
     end subroutine test
  end interface
  
  type(nmodel1) :: nmodel1_inst
  type(nmodel1cls) :: nmodel1cls_inst
  type(nmodel2cls) :: nmodel2cls_inst
  
  real :: clim(num_clim_par), param(num_parameters)
  character(len=*), parameter :: yasso_init_1cls = 'initialisation/yasso_init'
  print *, 'hoh hoh hoh'
    

  print *, '*** N RESOLVE ***'
  call test_resolve()

  call load_yasso_clim('weather/yasso_clim', clim)
  call load_params('parameters/yasso_parameters', param)

  print *, '*** NMODEL1 ***'
  call run_all_tests(nmodel1_inst, yasso_init_1cls)

  print *, '*** NMODEL1CLS ***'
  call run_all_tests(nmodel1cls_inst, yasso_init_1cls)

  call cascade_tests()
  
  !call nmodel1cls_inst%init_from_files('parameters/yasso_parameters', 'initialisation/yasso_init')
  !call test(nmodel1cls_inst)

  !call casc_from_files(nmodel2cls_inst, 'parameters/yasso_parameters', 'initialisation/yasso_init')

  
  !call init_yasso_from_files('parameters/yasso_parameters', 'initialisation/yasso_init', yasso_inst)

!### check N/C state ordering ####

  
contains
  
  subroutine run_all_tests(inst, filename_init)
    
    class(base_yasso_t) :: inst
    character(len=*), intent(in) :: filename_init

    real, allocatable :: cs(:), ns(:)
    real :: param(num_clim_par)
    
    call state_from_files(inst, filename_init, cs, ns)
    call run_test(test_zero, inst, cs, ns)
    call run_test(test_cn_state, inst, cs, ns)
    call run_test(test_cn_balance, inst, cs, ns)
    call run_test(test_map_to_basgra, inst, cs, ns)
    call run_test(test_matrix, inst, cs, ns)
    call run_test(test_self, inst, cs, ns)
    
  end subroutine run_all_tests

  subroutine run_test(testsub, inst, cs, ns)
    procedure(test) :: testsub
    class(base_yasso_t) :: inst
    real, intent(in) :: cs(:), ns(:)

    call setup(inst, param, cs, ns)
    call testsub(inst)
    call cleanup(inst)
    
  end subroutine run_test

  subroutine test_zero(inst)
    class(base_yasso_t) :: inst

    real :: cvec(inst%get_csize()), nvec(inst%get_nsize())

    cvec = 0.0
    nvec = 0.0
    call inst%cstate(set=cvec)
    call inst%nstate(set=nvec)
    
    call inst%decomp_demand(clim, 1.0)
    call inst%decomp_final(clim, 1.0, inst%get_n_demand())
    call inst%update_state(cvec, nvec)
    cvec = 1.0
    nvec = 1.0
    call inst%cstate(get=cvec)
    call check(all(cvec==0.0), 'Zeros C')
    call inst%nstate(get=nvec)
    call check(all(nvec==0.0), 'Zeros N')
  end subroutine test_zero
  
  subroutine test_self(inst)
    class(base_yasso_t) :: inst

    logical :: success

    if (inst%test()) then
       print *, 'PASS: self-test'
    else
       print *, 'FAIL: self-test'
    end if
    
  end subroutine test_self
  
  subroutine test_matrix(inst)
    class(base_yasso_t) :: inst

    real, allocatable :: matrix(:,:), ntmp(:), ctmp(:), cstate(:)
    real :: ndemand(inst%get_demand_size())
    integer :: csize
    
    csize = inst%get_csize()
    
    allocate(matrix(csize, csize), ctmp(csize), ntmp(inst%get_nsize()), cstate(csize))
    
    call inst%cmatrix(clim, matrix)
    call inst%cstate(get=cstate)
    call inst%decomp_demand(clim, 1.0)
    call inst%decomp_final(clim, 1.0, inst%get_n_demand())
    
    !tmp(1:csize) = matmul(matrix, cstate)
    !call inst%eval_tend(clim, 1e-9)
    ctmp = 0.0
    ntmp = 0.0
    call inst%update_state(ctmp, ntmp)
    call inst%cstate(get=ctmp)
    call check(all(almost(matmul(matrix,cstate), (ctmp(1:csize)-cstate))), 'Matrix')
    !print *, 'matrix', matrix

!!$    print *, 'tend from matrxi', matmul(matrix, cstate)
!!$    print *, 'tend from diff', (ctmp(1:csize)-cstate)
!!$    select type(inst)
!!$    type is(nmodel1)
!!$       print *, 'tend from state', inst%ctend
!!$    end select
!!$    !call check(all(almost(matmul(matrix, cstate)
    
  end subroutine test_matrix
  
  subroutine test_map_to_basgra(inst)
    class(base_yasso_t) :: inst

    real :: clitt, csomf, csoms, nlitt, nsomf, nsoms
    
    call inst%map_to_basgra(clitt, csomf, csoms, nlitt, nsomf, nsoms)
    call check(almost(inst%get_totc(), clitt+csomf+csoms), 'map_to_basgra: C')
    call check(almost(inst%get_norg(), nlitt+nsomf+nsoms), 'map_to_basgra: N')
    
  end subroutine test_map_to_basgra

  subroutine test_cn_state(inst)
    class(base_yasso_t) :: inst

    real :: cstate(inst%get_csize()), nstate(inst%get_nsize()), tmp(inst%get_csize())

    call inst%cstate(get=cstate)
    call check(almost(sum(cstate), inst%get_totc()), 'C-state total')

    call inst%nstate(get=nstate)
    call check(almost(sum(nstate), inst%get_norg()), 'N-state total')

    call random_number(cstate)
    call inst%cstate(set=cstate)
    call inst%cstate(get=tmp)
    call check(all(almost(cstate, tmp)), 'set C-state')
    
    call random_number(nstate)
    call inst%nstate(set=nstate)
    call inst%nstate(get=tmp)
    call check(all(almost(nstate, tmp(1:size(nstate)))), 'set N-state')
    
  end subroutine test_cn_state
  
  subroutine test_resolve()
    real :: demand_soil(4), demand_plant(2), alloc_soil(4), alloc_plant(2)
    real :: nmin
    
    demand_soil = 1.0
    nmin = 1e9
    demand_plant = 0.0
    call resolve_ndemand('demand_based', nmin, demand_soil, demand_plant=demand_plant, &
         alloc_soil=alloc_soil, alloc_plant=alloc_plant)
    call check(all(almost(alloc_soil, demand_soil)), 'Soil alloc no N-lim')
    call check(all(almost(alloc_plant, 0.0)), 'Plant alloc no demand')
    
    nmin = 0.0
    call resolve_ndemand('demand_based', nmin, demand_soil, demand_plant=demand_plant, &
         alloc_soil=alloc_soil, alloc_plant=alloc_plant)
    call check(all(almost(alloc_soil, 0.0)), 'Soil alloc no N avail')
    
    nmin = 0.01
    call resolve_ndemand('demand_based', nmin, demand_soil, demand_plant=demand_plant, &
         alloc_soil=alloc_soil, alloc_plant=alloc_plant)
    call check(all(almost(alloc_soil, nmin/size(alloc_soil))), 'Soil alloc N-lim')

    nmin = 1.0
    call resolve_ndemand('demand_based', nmin, demand_soil, demand_plant=(/sum(demand_soil)/), &
         alloc_soil=alloc_soil, alloc_plant=alloc_plant(1:1))
    ! 50/50 between plant and soil
    call check(almost(alloc_plant(1), sum(alloc_soil)), 'Soil alloc N-lim equal demand 1')
    ! 1/8 of total for each soil component
    call check(all(almost(alloc_soil, nmin/8.0)), 'Soil alloc N-lim equal demand 2')
    
  end subroutine test_resolve
  
  subroutine test_cn_balance(inst)
    class(base_yasso_t) :: inst

    real, dimension(:), allocatable :: cflux, nflux
    real :: dc, dn

    real, parameter :: dt = 1.0
    
    call inst%decomp_demand(clim, dt)
    call inst%decomp_final(clim, dt, inst%get_n_demand())

    allocate(cflux(inst%get_csize()), nflux(inst%get_nsize()))
    cflux = 1.0
    nflux = 1.0
    dc = inst%get_totc() !sum(inst%state(1:5))
    dn = inst%get_norg()
    call inst%update_state(cflux, nflux)
    dc = inst%get_totc() - dc !sum(inst%state(1:5)) - dc
    dn = inst%get_norg() - dn
    call check(almost(-inst%get_soilresp()*dt + sum(cflux), dc), 'C-balance')
    call check(almost(-inst%get_netmin_act()*dt + sum(nflux), dn), 'N-balance')
    
    call check(inst%get_soilresp() > 1e-3, 'Soilresp positive')
    
  end subroutine test_cn_balance

  subroutine cascade_tests()
    use nasso, only : num_c_pools
    
    type(cascademodel_t) :: inst
    integer, parameter :: num_cls = 4
    integer :: ind_cls, ii
    
    real :: tmp_cls(num_cls), tmp_cls_tend(num_cls), tmp_cstate(num_c_pools,num_cls), tmp_ctend(num_c_pools,num_cls), &
         ctmp(num_c_pools), net_immob, total_resp, cflux(num_c_pools*num_cls)
    
    call inst%init(param, (/num_cls/), (/(1.0, ind_cls = 1, num_cls)/))
    call random_number(tmp_cls)
    tmp_cls = tmp_cls + 1.0 ! avoid zeros
    
    call get_age_tend(tmp_cls, inst%get_ranges(), 0.1, tmp_cls_tend, num_cls)
    
    call check(almost(sum(tmp_cls_tend), 0.0), 'Agetend mass conservation')

    if (num_cls > 1) then
       call check(tmp_cls_tend(1) < 0, 'First cls tend < 0')
       call check(tmp_cls_tend(num_cls) > 0, 'Last cls tend > 0')
       tmp_cls = 0.0
       tmp_cls(1) = 1.0
       call get_age_tend(tmp_cls, inst%get_ranges(), 0.5, tmp_cls_tend, num_cls)
       call check(almost(tmp_cls_tend(1), -0.5), 'Unit pulse 1')
       call check(almost(tmp_cls_tend(2), 0.5), 'Unit pulse 2')
       call check(all(almost(tmp_cls_tend(3:), 0.0)), 'Unit pulse 3')
    end if

    cflux = 0.0
    call inst%cstate(set=cflux)
    call inst%cstate_cls(2, set=(/1.0, 0.0, 0.0, 0.0, 0.0/))
    tmp_cls = 0.0
    tmp_cls(2) = 0.5
    call inst%nstate(set=tmp_cls)


    ! A more heuristic test checking that the bins age in right direction and conserve
    ! mass.
    net_immob = 0.0
    total_resp = 0.0
    cflux = 0.0
    do ii = 1, 10*num_cls
       call inst%decomp_demand(clim, 0.5)
       net_immob = net_immob + sum(inst%get_n_demand())
       !print *, ii, net_immob, sum(inst%get_n_demand())
       !print *, '---', inst%get_n_demand()
       call inst%decomp_final(clim, 0.5, max(0.0, inst%get_n_demand()))
       total_resp = total_resp + inst%get_soilresp()
       call inst%update_state(cflux, tmp_cls*0)
    end do

    call inst%cstate_cls(num_cls, get=ctmp)
    !print *,  '*********'
    !print *, sum(ctmp), total_resp

    call check(sum(ctmp)+total_resp > 0.99, 'cascade timestep C transfer')

    call inst%nstate(get=tmp_cls)
    !print *, tmp_cls, net_immob
    call check(tmp_cls(num_cls) - net_immob > 0.99*0.5, 'cascade timestep N transfer')
    
    
  end subroutine cascade_tests
  
  subroutine state_from_files(yasso_inst, filename_init, cstate, nstate)
    use yassobase, only : init_from_files
    class(base_yasso_t) :: yasso_inst
    character(len=*), intent(in) :: filename_init
    real, dimension(:), allocatable, intent(out) :: cstate, nstate

    select type(yasso_inst)
    class default
       call init_from_files(yasso_inst, 'parameters/yasso_parameters', filename_init)
    end select
    allocate(cstate(yasso_inst%get_csize()), nstate(yasso_inst%get_nsize()))
    call yasso_inst%cstate(get=cstate)
    call yasso_inst%nstate(get=nstate)
    call yasso_inst%delete()
    
  end subroutine state_from_files

  subroutine setup(inst, params, cstate, nstate)
    class(base_yasso_t), intent(out) :: inst
    real, dimension(:), intent(in) :: params, cstate, nstate
    call inst%init(params)
    call inst%cstate(set=cstate)
    call inst%nstate(set=nstate)
  end subroutine setup

  subroutine cleanup(inst)
    class(base_yasso_t), intent(inout) :: inst
    call inst%delete()
  end subroutine cleanup

  
  subroutine check(cond, mesg)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: mesg
    if (.not. cond) then
       print *, 'FAIL:', mesg
    else
       print *, 'PASS:', mesg
    end if
    print *
  end subroutine check
  
  elemental logical function almost(val1, val2)
    real, intent(in) :: val1, val2
    real, parameter :: tol = 1e-5
    
    almost = abs(val1-val2) < tol
    
  end function almost
  
    
end program yassotest
