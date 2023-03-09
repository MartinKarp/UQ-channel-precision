module trunc
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: iso_c_binding
  implicit none
   
  type trunc_t 
     integer(kind=INT64) :: mask, rounding_bit
     integer(kind=INT64), allocatable :: work_array(:)
     integer :: n
     integer :: bits
   contains 
     procedure, pass(this) :: truncate => trunc_func
  end type trunc_t

  interface trunc_t
     module procedure trunc_init
  end interface trunc_t
  
contains

  function trunc_init(n, n_trunc_bits) result(this)
    type(trunc_t) :: this
    integer, intent(in) :: n, n_trunc_bits
    integer :: i

    call trunc_free(this)

    allocate(this%work_array(n))

    this%n = n
    this%bits = n_trunc_bits

    do i = 0, 63
       this%mask = ibset(this%mask,i)
       this%rounding_bit = ibclr(this%rounding_bit,i)
    end do
    do i = 0, n_trunc_bits-1
       this%mask = ibclr(this%mask,i)
    end do
    if(n_trunc_bits .gt. 0) this%rounding_bit = ibset(this%rounding_bit,n_trunc_bits-1)

  end function trunc_init

  subroutine trunc_free(this)
    type(trunc_t) :: this

    if(allocated(this%work_array)) deallocate(this%work_array)

  end subroutine trunc_free
 
  subroutine trunc_func(this,vector_to_trunc)
    class(trunc_t) :: this
    real(kind=REAL64), intent(inout) :: vector_to_trunc(this%n)
    integer(kind=INT64) :: bit_int, work_mask
    integer :: i

    this%work_array = transfer( vector_to_trunc, this%work_array, this%n)
    !Compute with rounding
    do i = 1, this%n
       bit_int = this%work_array(i) + this%rounding_bit
       work_mask = ieor(not(bit_int),this%mask)
       this%work_array(i) = iand(bit_int,work_mask)
    end do
    vector_to_trunc = transfer(this%work_array,vector_to_trunc, this%n)
  end subroutine trunc_func


end module trunc



module user
  use neko
  use trunc
  implicit none
  
  type(trunc_t) :: truncer
  integer, parameter :: n_pts = 9
  integer :: ele(n_pts) = (/1,1+1*18,1+2*18,1+3*18,1+4*18,1+5*18,1+6*18,1+7*18,1+8*18/)


contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%fluid_user_ic => user_ic
    u%user_mesh_setup => user_mesh_scale
    u%user_init_modules => user_initialize
    u%user_check => user_trunc
  end subroutine user_setup
  ! User defined initial condition

  ! User-defined initialization called just before time loop starts
  subroutine user_initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(param_t), intent(inout) :: params
    integer :: tstep, trunc_bits
    tstep = 0
    trunc_bits = 0 ! Number of bits to remove from the mantissa, bits left = 52-trunc_bits
    truncer = trunc_t(u%dof%size(),trunc_bits)
    call user_trunc(t, tstep, u, v, w, p, coef, params)
  end subroutine user_initialize

  subroutine user_trunc(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(param_t), intent(inout) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    integer :: ntot, i, n, el, e
    real(kind=rp) :: vv, sum_e1(1), e1, e2, sum_e2(1), oo

    n = u%dof%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(u%x, u%x_d, n, DEVICE_TO_HOST)
       call device_memcpy(v%x, v%x_d, n, DEVICE_TO_HOST)
       call device_memcpy(w%x, w%x_d, n, DEVICE_TO_HOST)
       call device_memcpy(p%x, p%x_d, n, DEVICE_TO_HOST)
    end if
    
    call truncer%truncate(u%x)
    call truncer%truncate(v%x)
    call truncer%truncate(w%x)
    call truncer%truncate(p%x)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(u%x, u%x_d, n, HOST_TO_DEVICE)
       call device_memcpy(v%x, v%x_d, n, HOST_TO_DEVICE)
       call device_memcpy(w%x, w%x_d, n, HOST_TO_DEVICE)
       call device_memcpy(p%x, p%x_d, n, HOST_TO_DEVICE)
    end if
    if (pe_rank .eq. 0) write(*,*) 'truncated ', truncer%bits,'bits'
    do e = 1, n_pts
       if (ele(e) .gt. coef%msh%offset_el .and. ele(e) .le. coef%msh%nelv+coef%msh%offset_el) then
          el = ele(e) - coef%msh%offset_el
          write(*,*) 'global-element-id:     ', 'x-coord   y-coord    z-coord    ',' u   v   w   p:',new_line('a'),&
                     ele(e), u%dof%x(4,4,4,el),u%dof%y(4,4,4,el),u%dof%z(4,4,4,el),&
                     u%x(4,4,4,el),v%x(4,4,4,el),w%x(4,4,4,el), p%x(4,4,4,el)
       end if
    end do
  end subroutine user_trunc


  ! Rescale mesh
  subroutine user_mesh_scale(msh)
    type(mesh_t), intent(inout) :: msh
    integer :: i, p, nvert
    real(kind=rp) :: d, y, viscous_layer, visc_el_h, el_h, center_el_h, dist_from_wall
    integer :: el_in_visc_lay, el_in_y
    real(kind=rp) :: pt(3)

    el_in_y = 18
    el_in_visc_lay = 2
    viscous_layer = 0.0888889
    el_h = 2.0_rp/el_in_y
    visc_el_h = viscous_layer/el_in_visc_lay
    center_el_h = (1.0_rp-viscous_layer)/(el_in_y/2-el_in_visc_lay)
    
    

    nvert = size(msh%points)
    do i = 1, nvert
       msh%points(i)%x(1) = pi*msh%points(i)%x(1) 
       y = msh%points(i)%x(2) 
       if ((1-abs(y)) .le. (el_in_visc_lay*el_h)) then
          dist_from_wall = (1-abs(y))/el_h*visc_el_h
       else
          dist_from_wall = viscous_layer + (1-abs(y)-el_in_visc_lay*el_h)/el_h*center_el_h
       end if
       if (y .gt. 0) msh%points(i)%x(2) = 1.0_rp - dist_from_wall
       if (y .lt. 0) msh%points(i)%x(2) = -1.0_rp + dist_from_wall
       msh%points(i)%x(3) = pi*msh%points(i)%x(3) 
    end do
    
  end subroutine user_mesh_scale


  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(param_t), intent(inout) :: params
    integer :: i
    real(kind=rp) :: uvw(3)

    do i = 1, u%dof%size()
       uvw = channel_ic(u%dof%x(i,1,1,1),u%dof%y(i,1,1,1),u%dof%z(i,1,1,1))
       u%x(i,1,1,1) = uvw(1)
       v%x(i,1,1,1) = uvw(2)
       w%x(i,1,1,1) = uvw(3)
    end do
  end subroutine user_ic
  
  function channel_ic(x, y, z) result(uvw)
    real(kind=rp) :: x, y, z
    real(kind=rp) :: uvw(3)
    real(kind=rp) :: ux, uy, uz, eps, Re_tau, yp, Re_b, alpha, beta
    real(kind=rp) :: C, k, kx, kz

      Re_tau = 300
      C      = 5.17
      k      = 0.41
      Re_b   = 5000

      yp = (1-y)*Re_tau
      if (y.lt.0) yp = (1+y)*Re_tau
      
      ! Reichardt function
      ux  = 1/k*log(1.0+k*yp) + (C - (1.0/k)*log(k)) * &
            (1.0 - exp(-yp/11.0) - yp/11*exp(-yp/3.0))
      ux  = ux * Re_tau/Re_b

      eps = 2e-2
      kx  = pi*4
      kz  = pi*4/3

      alpha = kx * 2*PI/8
      beta  = kz * 2*PI/4

      ! add perturbation to trigger turbulence 
      uvw(1)  = ux  + eps*beta  * sin(alpha*x)*cos(beta*z) 
      uvw(2)  =       eps       * sin(alpha*x)*sin(beta*z)
      uvw(3)  =      -eps*alpha * cos(alpha*x)*sin(beta*z)

  end function channel_ic

end module user
