program sys15f
    use, intrinsic :: iso_c_binding
    use types
    use fun

    implicit none

    namelist /param/ ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, &
        dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz, pitch

    integer(c_int) ne, nt, nz, i, j, freq_out
    real(c_double) tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz, pitch
    complex(c_double_complex), allocatable, target :: p(:, :), mean(:) !, oscill(:, :)
    real(c_double), allocatable, target :: tax(:), zax(:), u(:), eta(:, :), etag(:, :), w(:, :), f(:, :), w1(:, :)
    type(parametersf) paramf
    type(parametersp) paramp
    !complex(c_double_complex), pointer :: ff(:, :), pp(:, :)

    !interface
    !subroutine allocate_arrays(nz, nt, ne, f, p, u, t, z)!, oscill)
    !    use, intrinsic :: iso_c_binding
    !    implicit none
    !    integer, intent(in) :: nt, nz, ne
    !    complex(c_double_complex), allocatable, intent(inout) :: f(:, :), p(:, :)!, oscill(:, :)
    !    real(c_double), allocatable, intent(inout) :: t(:), z(:), u(:)
    !end subroutine allocate_arrays
    !
    !subroutine calc_u(u, zex, nz, zax)
    !    use, intrinsic :: iso_c_binding
    !    implicit none
    !    integer(c_int), intent(in) :: nz
    !    real(c_double), intent(in) :: zex, zax(nz)
    !    real(c_double), intent(out) :: u(:)
    !end subroutine
    !
    !subroutine read_param(ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz) bind(c,name='read_param')
    !    use, intrinsic :: iso_c_binding
    !    implicit none
    !    integer(c_int), intent(inout) :: ne
    !    real(c_double), intent(inout) :: tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz
    !end subroutine read_param

    !function dfdt(t, f, params) result(s)
    !    use, intrinsic :: iso_c_binding
    !    use types
    !    implicit none
    !    real(c_double) t
    !    complex(c_double_complex) f(:)
    !    complex(c_double_complex) s(size(f))
    !    type(parameters) params
    !end function dfdt

    !function dpdz(z, p, params) result(s)
    !    use, intrinsic :: iso_c_binding
    !    use types
    !    implicit none
    !    real(c_double) z
    !    complex(c_double_complex) p(:)
    !    complex(c_double_complex) s(size(p))
    !    type(parameters) params
    !
    !end function dpdz

    !subroutine ode4(dydt, y, neq, nt, t0, h, params)
    !    use, intrinsic :: iso_c_binding
    !    use types
    !    implicit none
    !    interface
    !        function dydt(t, y, par) result(s)
    !            use, intrinsic :: iso_c_binding
    !            use types
    !            implicit none
    !            real(c_double) t
    !            complex(c_double_complex) y(:)
    !            complex(c_double_complex) s(size(y))
    !            type(parameters) par
    !        end function dydt
    !    end interface
    !    type(parameters), intent(inout) :: params
    !    integer(c_int) nt, neq
    !    real(c_double) h, t0
    !    complex(c_double_complex), pointer :: y(:, :)
    !end subroutine ode4

    !subroutine ode4(dydt, y, neq, nt, t0, h, params)
    !use, intrinsic :: iso_c_binding
    !use types
    !implicit none
    !
        !!complex(c_double_complex), intent(out) :: y(:, :)
        !!complex(c_double_complex), pointer, intent(inout) :: y(:, :)
    !interface
    !    function dydt(t, y, par) result(s)
    !        use, intrinsic :: iso_c_binding
    !        use types
    !        implicit none
    !        real(c_double) t
    !        complex(c_double_complex) y(:)
    !        complex(c_double_complex) s(size(y))
    !        type(parameters) par
    !    end function dydt
    !end interface
    !
    !type(parameters), intent(inout) :: params
    !integer(c_int) nt, neq, i
    !real(c_double) h, t0, t
    !complex(c_double_complex), pointer :: y(:, :)
    !complex(c_double_complex) s1(size(y, 1)), s2(size(y, 1)), s3(size(y, 1)), s4(size(y, 1)), v(size(y, 1))

    !end subroutine ode4
    !end interface

    call read_param(ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, &
                    dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz, pitch)
    write (*, nml=param)

    !par%gamma = 0.1
    !par%omega = 1

    nt = tend/dt + 1
    nz = zex/dz + 1

    !call allocate_arrays(nz, nt, ne, f, p1, p2, u, tax, zax, oscill)
    call allocate_arrays(nz, nt, ne, f, p, u, tax, zax, mean, eta, etag, w, w1)

    f(1, 1) = f10
    f(3, 1) = f20
    f(5, 1) = f30

    paramp%ne = ne
    paramp%nz = nz
    paramf%nt = nt
    paramf%tend = tend
    paramp%zex = zex
    paramf%q(1) = q1
    paramf%q(2) = q2
    paramf%q(3) = q3
    paramf%i(1) = i1
    paramf%i(2) = i2
    paramf%th(1) = th1
    paramf%th(2) = th2
    paramf%a(1) = a1
    paramf%a(2) = a2
    paramp%dtr(1) = dtr1
    paramp%dtr(2) = dtr2
    paramf%dcir(1) = dcir1
    paramf%dcir(2) = dcir2
    paramf%r(1) = r1
    paramf%r(2) = r2
    paramf%dt = dt
    paramp%dz = dz
    paramf%f => f
    paramp%p => p
    paramp%mean => mean
    paramp%eta => eta
    paramp%etag => etag
    paramf%pitch = pitch
    paramf%w => w1

    !freq_out = 777
    !open (freq_out, file='W.dat')
    !paramf%freq_out = freq_out

    do i = 1, nt
        tax(i) = (i - 1)*dt
    end do

    do i = 1, nz
        zax(i) = (i - 1)*dz
    end do

    do i = 1, ne
        paramp%p(i, 1) = exp(ic*(i - 1)/dble(ne)*2*pi)
        paramp%p(ne + i, 1) = exp(ic*(i - 1)/dble(ne)*2*pi)
        !print *, exp(ic*(i - 1)/dble(ne)*2*pi)
    end do

    call calc_u(u, zex, nz, zax)

    paramp%u => u

    !oscill(1, 1) = 0.0d0 + ic*1.0d0

    write (*, '(/)')

    call ode4f(dfdt, paramf%f, 3, nt, 0.0d0, dt, paramf, paramp)

    !close (freq_out)

    do i = 2, nt
        do j = 1, 3
            w(j, i - 1) = imag(log(f(2*j - 1, i)*exp(ic*f(2*j, i))/(f(2*j - 1, i - 1)*exp(ic*f(2*j, i - 1)))))/dt
        end do
    end do

    write (*, '(/)')

    pause

    open (1, file='F.dat', err=101)
    do i = 1, nt
        write (1, '(4e17.8)') tax(i), abs(f(1, i)), abs(f(3, i)), abs(f(5, i))
    end do
    close (1)
    open (2, file='E.dat', err=101)
    do i = 1, nt
        write (2, '(5e17.8)') tax(i), eta(1, i), etag(1, i), eta(2, i), etag(2, i)
    end do
    close (2)
    open (3, file='W1.dat', err=101)
    do i = 1, nt - 1
        write (3, '(4e17.8)') tax(i + 1), w1(1, i), w1(2, i), w1(3, i)
    end do
    close (3)
    open (4, file='W.dat', err=101)
    do i = 1, nt - 1
        write (4, '(4e17.8)') tax(i + 1), w(1, i), w(2, i), w(3, i)
    end do
    close (4)
    stop
101 print *, 'error of file open.'
    pause
    stop
102 print *, 'error of file reading.'
    pause
    stop
103 print *, 'error of file writing.'
    pause
    stop
end program

!subroutine allocate_arrays(nz, nt, ne, f, p, u, t, z)!, oscill)
!    use, intrinsic :: iso_c_binding
!    implicit none
!
!    integer, intent(in) :: nz, nt, ne
!    complex(c_double_complex), allocatable, intent(inout) :: f(:, :), p(:, :)!, oscill(:, :)
!    real(c_double), allocatable, intent(inout) :: t(:), z(:), u(:)
!
!    integer(c_int) err_alloc
!
!    !allocate (f(nt, 3), p1(nz, ne), p2(nz, ne), u(nz), t(nt), z(nz), oscill(nt, 1), stat=err_alloc)
!    allocate (f(3, nt), p(2*ne, nz), u(nz), t(nt), z(nz), stat=err_alloc)
!
!    if (err_alloc /= 0) then
!        print *, "allocation error"
!        pause
!        stop
!    end if
!end subroutine allocate_arrays
!
!subroutine deallocate_arrays()
!    use, intrinsic :: iso_c_binding
!    implicit none
!
!    integer(c_int) err_dealloc
!    complex(c_double_complex), allocatable :: f(:, :), p(:, :), u(:)
!
!    deallocate (f, p, u, stat=err_dealloc)
!
!    if (err_dealloc /= 0) then
!        print *, "deallocation error"
!        pause
!        stop
!    end if
!end subroutine deallocate_arrays
!
!subroutine read_param(ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz) bind(c,name='read_param')
!    use, intrinsic :: iso_c_binding
!    implicit none
!
!    namelist /param/ ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz
!
!    integer(c_int), intent(inout) :: ne
!      real(c_double), intent(inout) :: tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz
!
!    open (unit=1, file='input_fortran.dat', status='old', err=101)
!    read (unit=1, nml=param, err=102)
!    close (unit=1)
!
!    return
!101 print *, "error of file open"; pause; stop
!102 print *, 'error of reading file "input_fortran.dat"'; pause; stop
!end subroutine read_param

!subroutine ode4(dydt, y, neq, nt, t0, h, params)
!    use, intrinsic :: iso_c_binding
!    use types
!    implicit none
!
!    !complex(c_double_complex), intent(out) :: y(:, :)
!    !complex(c_double_complex), pointer, intent(inout) :: y(:, :)
!    interface
!        function dydt(t, y, par) result(s)
!            use, intrinsic :: iso_c_binding
!            use types
!            implicit none
!            real(c_double) t
!            complex(c_double_complex) y(:)
!            complex(c_double_complex) s(size(y))
!            type(parameters) par
!        end function dydt
!    end interface
!
!    type(parameters), intent(inout) :: params
!    integer(c_int) nt, i, neq
!    real(c_double) h, t0, t
!    complex(c_double_complex), pointer :: y(:, :)
!    complex(c_double_complex) s1(size(y, 1)), s2(size(y, 1)), s3(size(y, 1)), s4(size(y, 1)), v(size(y, 1))
!
!    do i = 1, nt - 1
!        v = y(:, i)
!        t = t0 + (i - 1)*h
!        s1(:) = dydt(t, v, params)
!        s2(:) = dydt(t + h/2, v + h*s1(:)/2, params)
!        s3(:) = dydt(t + h/2, v + h*s2(:)/2, params)
!        s4(:) = dydt(t + h, v + h*s3(:), params)
!        y(:, i + 1) = v + h*(s1(:) + 2*s2(:) + 2*s3(:) + s4(:))/6
!    end do
!end subroutine ode4

!function dfdt(t, f, params) result(s)
!    use, intrinsic :: iso_c_binding
!    use types
!    implicit none
!
!    interface
!        function dpdz(z, p, par) result(s)
!            use, intrinsic :: iso_c_binding
!            use types
!            implicit none
!            real(c_double) z
!            complex(c_double_complex) p(:)
!            complex(c_double_complex) s(size(p))
!            type(parameters) par
!        end function dpdz
!        subroutine ode4(dydt, y, neq, nt, t0, h, params)
!            use, intrinsic :: iso_c_binding
!            use types
!            implicit none
!            interface
!                function dydt(t, y, par) result(s)
!                    use, intrinsic :: iso_c_binding
!                    use types
!                    implicit none
!                    real(c_double) t
!                    complex(c_double_complex) y(:)
!                    complex(c_double_complex) s(size(y))
!                    type(parameters) par
!                end function dydt
!            end interface
!            type(parameters), intent(inout) :: params
!            integer(c_int) nt, neq
!            real(c_double) h, t0
!            complex(c_double_complex), pointer :: y(:, :)
!        end subroutine ode4
!    end interface
!
!    integer(c_int) ne, nz
!    complex(c_double_complex), save :: ic = (0.0D0, 1.0D0)
!    real(c_double) t, dz
!    complex(c_double_complex) f(3), s(3)
!    type(parameters) params
!    !complex(c_double_complex), pointer :: p(:, :)
!
!    ne = params%ne
!    nz = params%nz
!    dz = params%dz
!
!    call ode4(dpdz, params%p, ne, nz, 0.0d0, dz, params)
!
!    s = 10
!
!    ![~, Pv] = osol(@ (z, Pv) MomentumODEv(z, Pv, F, SU, Par.D1, Par.D2, Idx), ZetaAxis, P0v, optionsP);
!    !P1(:, :) = Pv(:, Idx.Re1) + ic*Pv(:, Idx.Im1);
!    !P2(:, :) = Pv(:, Idx.Re2) + ic*Pv(:, Idx.Im2);
!    !Xi1 = Xi(SU(ZetaAxis), P1(:, :), ZetaAxis(:));
!    !Xi2 = Xi(SU(ZetaAxis), P2(:, :), ZetaAxis(:));
!    !s(1) = (ic*Par.I1*Xi1 - F(1))*(Par.Q3/Par.Q1) + (2*Par.R1*(Par.Q3/Par.Q1))*exp(-ic*Par.Th1)*F(3) - ic*Par.d1*2*Par.Q3*F(1);
!    !s(2) = (ic*Par.I2*Xi2 - F(2))*(Par.Q3/Par.Q2) + (2*Par.R2*(Par.Q3/Par.Q2))*exp(-ic*Par.Th2)*F(3) - ic*Par.d2*2*Par.Q3*F(2);
!    !s(3) = -F(3) + Par.A1*F(1) + Par.A2*F(2);
!end function dfdt

!function dpdz(z, p, params) result(s)
!    use, intrinsic :: iso_c_binding
!    use types
!    implicit none
!
!    type(parameters) params
!    integer(c_int) i
!    real(c_double) u, zex, z
!    complex(c_double_complex), save :: ic = (0.0D0, 1.0D0)
!    complex(c_double_complex), pointer :: p(:)
!    complex(c_double_complex) s(size(p,1)), f
!
!    zex = params%zex
!
!    do i = 1, 2
!        u = exp(-3*((z - zex/2)/(zex/2))**2)
!        s(:) = ic*(f*u - (params%dtr(i) + abs(p)**2 - 1)*p)
!    end do
!end function dpdz

!function dfdt(t, f, p) result(s)
!    use, intrinsic :: iso_c_binding
!    use types
!    implicit none
!
!    complex(c_double_complex) :: ic = (0.0D0, 1.0D0)
!    real(c_double) t
!    complex(c_double_complex) f(:), s(size(f))
!    type(parameters) p
!
!    s(:) = imag(f) - ic*(2*p%gamma*imag(f) + p%omega**2*real(f));
!end function dfdt

!subroutine calc_u(u, zex, nz, zax)
!    use, intrinsic :: iso_c_binding
!    implicit none
!
!    integer(c_int), intent(in) :: nz
!    real(c_double), intent(in) :: zex, zax(nz)
!    real(c_double), intent(out) :: u(:)
!
!    integer(c_int) i
!
!    do i = 1, nz
!        u(i) = exp(-3*((zax(i) - zex/2)/(zex/2))**2)
!    end do
!
!end subroutine

