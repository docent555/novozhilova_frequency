module types
    use, intrinsic :: iso_c_binding

    type parametersf
        integer(c_int) nt, freq_out
        real(c_double) tend, q(3), i(2), th(2), a(2), dcir(2), r(2), f10, f20, f30, dt, pitch
        real(c_double), pointer :: u(:), f(:, :), w(:, :)
        !real(c_double) gamma, omega
    end type parametersf

    type parametersp
        integer(c_int) ne, nz
        real(c_double) zex, dtr(2), dz
        real(c_double), pointer :: zax(:), u(:), eta(:, :), etag(:, :)
        complex(c_double_complex), pointer :: p(:, :), mean(:)
        complex(c_double_complex) f(2)
        !real(c_double) gamma, omega
    end type parametersp
end module types
