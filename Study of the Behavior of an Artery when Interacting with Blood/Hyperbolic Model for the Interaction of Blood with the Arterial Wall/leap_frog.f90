
subroutine leap_frog(u_out,dt,dx,gamma,Nt,Nx,f,u_0) bind(C, name="leap_frog")
    use iso_c_binding
    implicit none

    integer(c_int), intent(in), value :: Nt,Nx,f
    real(c_double), intent(in), value :: dt,dx,gamma
    real(c_double), intent(in) :: u_0(2,Nx)
    real(c_double), intent(inout) :: u_out(Nt,Nx)
    real(c_double) :: u(Nx),u_old(Nx),u_old2(Nx)

    real(c_double) :: lambda,t,x,v0(Nx),phi
    integer(c_int) :: i,j,n_t

    t = 0.0d0
    x = 0.0d0

    v0(:) = u_0(2,:)
    u_old2(:) = u_0(1,:)
    u_old = 0.0d0
    u = 0.0d0
    u_out = 0.0d0
    u_out(1,:) = u_old2(:)

    phi = (gamma*dt/dx)**2


    do i=2,Nx-1
        x = (i-1)*dx
        u_old(i) = u_old2(i) + dt*v0(i) +(phi/2.0d0)*(u_old2(i+1) - 2.0d0*u_old2(i) + u_old2(i-1)) + (dt**2/2)*func(t,x,f,gamma)
        u_out(2,i) = u_old(i)
    end do


    t = t + dt

    do n_t=3,Nt
        x = 0
        do j=2,Nx-1
            x = (j-1)*dx
            u(j) = 2.0d0*u_old(j) - u_old2(j) + phi*(u_old(j+1)-2.0d0*u_old(j)+u_old(j-1)) + (dt**2)*func(t,x,f,gamma)
            u_out(n_t,j) = u(j)

        end do

        u_old2(:) = u_old(:)
        u_old(:) = u(:)

        t = t + dt
    end do

    contains 



        real(c_double) function func(t,x,f,gamma)
            implicit none

            real(c_double), intent(in) :: x,t,gamma
            integer(c_int), intent(in) :: f
            real(c_double), parameter :: pi = 4.0d0*atan(1.0d0)

            real(c_double) :: dP,rho_w,H,omega_0,zeta,a,b,fx

            rho_w = 10**3
            H = 3.0d-4
            b = 133.32
            a = 10*b
            dP = 0.25*b
            omega_0 = 2.0d0*pi/0.8

            zeta = dP/(zeta)

            if (f==1) then
                fx = (1+(pi*gamma)**2)*exp(-t)*sin(pi*x)
            else if (f==2) then
                fx = x*dP*sin(omega_0*t)/(rho_w*H)
            end if 
            func = fx

        end function func

end subroutine leap_frog