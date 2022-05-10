module fun
    use types
    use ifcore

    complex(c_double_complex), parameter :: ic = (0.0d0, 1.0d0)
    real(c_double), parameter :: pi = 2.0D0*dacos(0.0D0)

contains
    subroutine allocate_arrays(nz, nt, ne, f, p, u, t, z, mean, eta, etag, w, w1)!, oscill)
        use, intrinsic :: iso_c_binding
        implicit none

        integer, intent(in) :: nz, nt, ne
        complex(c_double_complex), allocatable, intent(inout) :: p(:, :), mean(:)!, oscill(:, :)
        real(c_double), allocatable, intent(inout) :: t(:), z(:), u(:), eta(:, :), etag(:, :), w(:, :), f(:, :), w1(:, :)

        integer(c_int) err_alloc

        !allocate (f(nt, 3), p1(nz, ne), p2(nz, ne), u(nz), t(nt), z(nz), oscill(nt, 1), stat=err_alloc)
        allocate (f(6, nt), p(2*ne, nz), u(nz), t(nt), z(nz), mean(nz), eta(2, nt), etag(2, nt), w(3, nt - 1), w1(3, nt - 1), stat=err_alloc)

        if (err_alloc /= 0) then
            print *, "allocation error"
            pause
            stop
        end if
    end subroutine allocate_arrays

    subroutine deallocate_arrays()
        use, intrinsic :: iso_c_binding
        implicit none

        integer(c_int) err_dealloc
        complex(c_double_complex), allocatable :: f(:, :), p(:, :), u(:), t(:), z(:), mean(:)

        deallocate (f, p, u, t, z, mean, stat=err_dealloc)

        if (err_dealloc /= 0) then
            print *, "deallocation error"
            pause
            stop
        end if
    end subroutine deallocate_arrays

    subroutine read_param(ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz, pitch) bind(c,name='read_param')
        use, intrinsic :: iso_c_binding
        implicit none

namelist /param/ ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz, pitch

        integer(c_int), intent(inout) :: ne
        real(c_double), intent(inout) :: tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz, pitch

        open (unit=1, file='input_fortran.dat', status='old', err=101)
        read (unit=1, nml=param, err=102)
        close (unit=1)

        return
101     print *, "error of file open"; pause; stop
102     print *, 'error of reading file "input_fortran.dat"'; pause; stop
    end subroutine read_param

    subroutine ode4f(dydt, y, neq, nt, t0, h, paramf, paramp)
        import
        implicit none

        interface
            function dydt(t, y, parf, parp) result(s)
                use, intrinsic :: iso_c_binding
                use types
                implicit none
                real(c_double) t, y(:), s(size(y))
                type(parametersf) parf
                type(parametersp) parp
            end function dydt
        end interface

        type(parametersf), intent(inout) :: paramf
        type(parametersp), intent(inout) :: paramp
        integer(c_int) nt, neq, i, j
        real(c_double) h, t0, t
        real(c_double), pointer :: y(:, :)
        real(c_double) s1(size(y, 1)), s2(size(y, 1)), s3(size(y, 1)), s4(size(y, 1)), v(size(y, 1))
        logical(4) pressed
        character(1) key
        integer(c_int), parameter :: esc = 27

        !solve eq. at t=0
        paramp%f(1) = y(1, 1)*exp(ic*y(2, 1))
        paramp%f(2) = y(3, 1)*exp(ic*y(4, 1))
        call ode4p(dpdz, paramp%p, paramp%ne, paramp%nz, 0.0d0, paramp%dz, paramp)
        paramp%eta(:, 1) = eff(paramp%p(:, paramp%nz), paramp%ne)
        paramp%etag(:, 1) = paramf%pitch**2/(paramf%pitch**2 + 1)*paramp%eta(:, 1)

        do i = 1, nt - 1
            v = y(:, i)
            t = t0 + (i - 1)*h
            s1(:) = dydt(t, v, paramf, paramp)
            s2(:) = dydt(t + h/2, v + h*s1(:)/2, paramf, paramp)
            s3(:) = dydt(t + h/2, v + h*s2(:)/2, paramf, paramp)
            s4(:) = dydt(t + h, v + h*s3(:), paramf, paramp)
            y(:, i + 1) = v + h*(s1(:) + 2*s2(:) + 2*s3(:) + s4(:))/6

            paramp%eta(:, i + 1) = eff(paramp%p(:, paramp%nz), paramp%ne)
            paramp%etag(:, i + 1) = paramf%pitch**2/(paramf%pitch**2 + 1)*paramp%eta(:, i + 1)

            write (*, '(a,f12.7,a,f10.7,a,f10.7,a,f10.7,a,f10.7,a,f10.7,\,a)') 'Time = ', t + h, '   |F1| = ', abs(y(1, i + 1)), '   |F2| = ', abs(y(3, i + 1)), '   |F3| = ', abs(y(5, i + 1)), &
                '   Eff1 = ', paramp%eta(1, i + 1), '   Eff2 = ', paramp%eta(2, i + 1), char(13)

            pressed = peekcharqq()
            if (pressed) then
                key = getcharqq()
                if (ichar(key) .eq. esc) then
                    write (*, '(/,a)') 'Quit?'
                    key = getcharqq()
                    if (ichar(key) .eq. 121 .or. ichar(key) .eq. 89) then
                        open (1, file='F.dat')
                        do j = 1, i + 1
                            write (1, '(4e17.8)') (j - 1)*h, abs(y(1, j)), abs(y(2, j)), abs(y(3, j))
                        end do
                        close (2)
                        open (2, file='E.dat')
                        do j = 1, i + 1
                           write (2, '(5e17.8)') (j - 1)*h, paramp%eta(1, j), paramp%etag(1, j), paramp%eta(2, j), paramp%etag(2, j)
                        end do
                        close (2)
                        stop
                    end if
                end if
            end if
        end do
    end subroutine ode4f

    function eff(pex, ne) result(eta)
        import

        integer(c_int) ne
        real(c_double) eta(2)
        complex(c_double_complex), intent(in) :: pex(:)

        eta(1) = 1 - sum(abs(pex(1:ne))**2)/ne
        eta(2) = 1 - sum(abs(pex(ne + 1:2*ne))**2)/ne
    end function eff

    subroutine ode4p(dydt, y, neq, nz, z0, h, params)
        import
        implicit none

        interface
            function dydt(t, y, param) result(s)
                use, intrinsic :: iso_c_binding
                use types
                implicit none
                real(c_double) t
                complex(c_double_complex) y(:), s(size(y)), f
                type(parametersp) param
            end function dydt
        end interface

        type(parametersp), intent(inout) :: params
        integer(c_int) nz, neq, i
        real(c_double) h, z0, z
        complex(c_double_complex), pointer :: y(:, :)
        complex(c_double_complex) s1(size(y, 1)), s2(size(y, 1)), s3(size(y, 1)), s4(size(y, 1)), v(size(y, 1))

        do i = 1, nz - 1
            v = y(:, i)
            z = z0 + (i - 1)*h
            s1(:) = dydt(z, v, params)
            s2(:) = dydt(z + h/2, v + h*s1(:)/2, params)
            s3(:) = dydt(z + h/2, v + h*s2(:)/2, params)
            s4(:) = dydt(z + h, v + h*s3(:), params)
            y(:, i + 1) = v + h*(s1(:) + 2*s2(:) + 2*s3(:) + s4(:))/6
        end do
        !do i=1,size(y,1)
        !      print *, y(i,1)
        !enddo
        !pause
    end subroutine ode4p

    function dpdz(z, p, params) result(s)
        import
        implicit none

        type(parametersp) params
        integer(c_int) i, ne, idx(params%ne)
        real(c_double) u, zex, z
        complex(c_double_complex), parameter :: ic = (0.0D0, 1.0D0)
        complex(c_double_complex) p(:), s(size(p, 1)), f(2)

        zex = params%zex
        ne = params%ne
        f(1) = params%f(1)
        f(2) = params%f(2)

        do i = 1, 2
            u = exp(-3*((z - zex/2)/(zex/2))**2)
            idx = (/(i - 1)*ne + 1:ne + (i - 1)*ne/)
            s(idx) = ic*(f(i)*u - (params%dtr(i) + abs(p(idx))**2 - 1)*p(idx))
            !write(*, '(a,f17.8,/, a,f17.8,/, a,f17.8,/ a,f17.8,/ a,f17.8,/ a,f17.8,/)') 'z = ', z, 'p(3) = ', real(p(3)), 'p(16) = ' , real(p(16)), 'F = ' , real(f(1)), 'u(z) = ', u, 'D = ', params%dtr(1)
        end do
    end function dpdz

    complex(c_double_complex) function xi(u, p, dz, mean)
        integer(c_int) m, n, i
        real(c_double) dz
        complex(c_double_complex) u(:), p(:, :), mean(:)

        m = size(p, 1)
        n = size(mean)
        do i = 1, n
            mean(i) = sum(p(:, i), 1)/m
        end do

        mean = dconjg(u)*mean

        xi = (0.5d0*(mean(1) + mean(2)) + sum(mean(2:n - 1)))*dz
    end function

    function dfdt(t, f, paramf, paramp) result(s)
        implicit none

        !interface
        !    function dpdz(z, p, par) result(s)
        !        use, intrinsic :: iso_c_binding
        !        use types
        !        implicit none
        !        real(c_double) z
        !        complex(c_double_complex) p(:)
        !        complex(c_double_complex) s(size(p))
        !        type(parameters) par
        !    end function dpdz
        !    subroutine ode4(dydt, y, neq, nt, t0, h, params)
        !        use, intrinsic :: iso_c_binding
        !        use types
        !        implicit none
        !        interface
        !            function dydt(t, y, par) result(s)
        !                use, intrinsic :: iso_c_binding
        !                use types
        !                implicit none
        !                real(c_double) t
        !                complex(c_double_complex) y(:)
        !                complex(c_double_complex) s(size(y))
        !                type(parameters) par
        !            end function dydt
        !        end interface
        !        type(parameters), intent(inout) :: params
        !        integer(c_int) nt, neq
        !        real(c_double) h, t0
        !        complex(c_double_complex), pointer :: y(:, :)
        !    end subroutine ode4
        !end interface

        integer(c_int) :: ne, nz, i, j, iter_num = 1, time_num = 1
        real(c_double) t, dz, f(:), s(size(f)), &
            x1r, x1i, q31, i1, r1, th1, dcir1, cos1, sin1, &
            x2r, x2i, q32, i2, r2, th2, dcir2, cos2, sin2, q3, &
            f1, f2, f3, phi1, phi2, phi3, a1, a2
        complex(c_double_complex) x1, x2
        type(parametersf) paramf
        type(parametersp) paramp
        !logical(4) pressed
        !character(1) key
        !integer(c_int), parameter :: esc = 27

        ne = paramp%ne
        nz = paramp%nz
        dz = paramp%dz
        paramp%f(1) = f(1)*exp(ic*f(2))
        paramp%f(2) = f(3)*exp(ic*f(4))

        call ode4p(dpdz, paramp%p, ne, nz, 0.0d0, dz, paramp)

        ! test
        !paramp%f(1) = 0.033397679068901
        !paramp%f(2) = 0.033397679068901
        !call ode4p(dpdz, paramp%p, ne, nz, 0.0d0, dz, paramp)
        !paramp%eta(:, 1) = eff(paramp%p(:, paramp%nz), paramp%ne)
        !paramp%etag(:, 1) = paramf%pitch**2/(paramf%pitch**2 + 1)*paramp%eta(:, 1)
        !print *, paramp%eta(1, 1)
        !pause
        !end test

        !open (1, file='test.dat')
        !do i = 1, nz
        !    do j = 1, ne
        !        write (1, '(1p2e17.8\)') real(paramp%p(j, i)), imag(paramp%p(j, i))
        !    end do
        !    write (1, '(/\)')
        !end do
        !close (1)
        !pause

        x1 = xi(dcmplx(paramp%u), paramp%p(1:ne, :), dz, paramp%mean)
        x2 = xi(dcmplx(paramp%u), paramp%p(ne + 1:2*ne, :), dz, paramp%mean)

        x1r = real(x1)
        x1i = imag(x1)
        x2r = real(x2)
        x2i = imag(x2)

        f1 = f(1)
        phi1 = f(2)
        f2 = f(3)
        phi2 = f(4)
        f3 = f(5)
        phi3 = f(6)

        q31 = paramf%q(3)/paramf%q(1)
        i1 = paramf%i(1)
        r1 = paramf%r(1)
        th1 = paramf%th(1)
        dcir1 = paramf%dcir(1)
        cos1 = cos(f(2))
        sin1 = sin(f(2))

        q32 = paramf%q(3)/paramf%q(2)
        i2 = paramf%i(2)
        r2 = paramf%r(2)
        th2 = paramf%th(2)
        dcir2 = paramf%dcir(2)
        cos2 = cos(f(4))
        sin2 = sin(f(4))

        q3 = paramf%q(3)
        a1 = paramf%a(1)
        a2 = paramf%a(2)

        s(1) = (-f1 + i1*(-x1i*cos1 + x1r*sin1) + 2*r1*f3*cos(phi3 - phi1 - th1))*q31
        s(2) = -2*dcir1*q3 + (i1/f1*(x1r*cos1 + x1i*sin1) + 2*r1*(f3/f1)*sin(phi3 - phi1 - th1))*q31

        s(3) = (-f2 + i2*(-x2i*cos2 + x2r*sin2) + 2*r2*f3*cos(phi3 - phi2 - th2))*q32
        s(4) = -2*dcir2*q3 + (i2/f2*(x2r*cos2 + x2i*sin2) + 2*r2*(f3/f2)*sin(phi3 - phi2 - th2))*q32

        s(5) = -f3 + a1*f1*cos(phi1 - phi3) + a2*f2*cos(phi2 - phi3)
        s(6) = a1*f1/f3*sin(phi1 - phi3) + a2*f2/f3*sin(phi2 - phi3)

        if (mod(iter_num, 4) .eq. 0) then
            paramf%w(1, time_num) = s(2)
            paramf%w(2, time_num) = s(4)
            paramf%w(3, time_num) = s(6)
            time_num = time_num + 1
        end if
        iter_num = iter_num + 1

        !s(1) = (ic*paramf%i(1)*x1 - f(1))*(paramf%q(3)/paramf%q(1)) + (2*paramf%r(1)*(paramf%q(3)/paramf%q(1)))*exp(-ic*paramf%th(1))*f(3) - ic*paramf%dcir(1)*2*paramf%q(3)*f(1)
        !s(2) = (ic*paramf%i(2)*x2 - f(2))*(paramf%q(3)/paramf%q(2)) + (2*paramf%r(2)*(paramf%q(3)/paramf%q(1)))*exp(-ic*paramf%th(2))*f(3) - ic*paramf%dcir(2)*2*paramf%q(3)*f(2)
        !s(3) = -f(3) + paramf%a(1)*f(1) + paramf%a(2)*f(2)

        !do i=1,3
        !      write (*, '(2f17.8)') s(i)
        !enddo
        !pause

        !pressed = peekcharqq()
        !if (pressed) then
        !    key = getcharqq()
        !    if (ichar(key) .eq. esc) then
        !        print *, pressed
        !        pause
        !    end if
        !    print *, 'qq'
        !end if
    end function dfdt

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

    subroutine calc_u(u, zex, nz, zax)
        import
        implicit none

        integer(c_int), intent(in) :: nz
        real(c_double), intent(in) :: zex, zax(nz)
        real(c_double), intent(out) :: u(:)

        integer(c_int) i

        do i = 1, nz
            u(i) = exp(-3*((zax(i) - zex/2)/(zex/2))**2)
        end do

    end subroutine

end module fun
