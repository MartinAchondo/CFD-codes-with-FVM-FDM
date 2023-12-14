
subroutine newmark(u_out,dt,dx,gamma,Nt,Nx,f,beta,theta,u_0) bind(C, name="newmark")
    use iso_c_binding
    implicit none

    integer(c_int), intent(in), value :: Nt,Nx,f
    real(c_double), intent(in), value :: dt,dx,gamma,theta,beta
    real(c_double), intent(in) :: u_0(2,Nx)
    real(c_double), intent(inout) :: u_out(Nt,Nx)
    
    real(c_double) :: u(Nx),u_old(Nx),w(Nx),w_old(Nx),v(Nx),v_old(Nx)
    real(c_double) :: b(Nx-2), xx(Nx-2)
    real(c_double) :: t,x,phi,c1,c2,lambda
    integer(c_int) :: j,n_t

    t = 0.0d0
    x = 0.0d0

    xx = 0.0d0
    b = 0.0d0
    v_old(:) = u_0(2,:)
    u_old(:) = u_0(1,:)
 
    w_old = 0.0d0
    u = 0.0d0
    w = 0.0d0
    v = 0.0d0
    u_out = 0.0d0
    u_out(1,:) = u_old(:)

    lambda = dt/dx
    phi = (gamma*lambda)**2
    c1 = -phi*beta
    c2 = 1.0d0 + 2.0d0*phi*beta

    t = t 

    do j=2,Nx-1
        w_old(j) = u_old(j+1)-2.0d0*u_old(j)+u_old(j-1)
    end do

    do n_t=2,Nt

        x = 0
        do j=2,Nx-1
            x = (j-1)*dx

            b(j-1) = u_old(j) + dt*v_old(j)  + phi*(0.5d0-beta)*w_old(j) + &
            &dt**2*beta*func(t+dt,x,f,gamma) + dt**2*(0.5d0-beta)*func(t,x,f,gamma)
            
        end do

        call Thomas(xx,Nx-2,c1,c2,c1,b)

        do j=2,Nx-1
            u(j) = xx(j-1)
        end do

        do j=2,Nx-1
            x = (j-1)*dx
            w(j) = u(j+1)-2.0d0*u(j)+u(j-1)

            v(j) = v_old(j) + (phi/dt)*( (1-theta)*w_old(j) + theta*w(j)) +&
            & dt*(theta*func(t+dt,x,f,gamma) + (1-theta)*func(t,x,f,gamma))
        end do

        w_old(:) = w(:)
        v_old(:) = v(:)
        u_old(:) = u(:)

        u_out(n_t,:) = u(:)

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


        subroutine Thomas(x,N,e,a,c,b)
            implicit none
            
            integer(c_int), intent(in) :: N
            real(c_double), intent(in) :: a,c,e,b(N)
            real(c_double), intent(out) :: x(N)
            integer(c_int) :: k
            real(c_double) :: beta(N),alpha(N),y(N)

            alpha(1) = a

            do k=2,N
                beta(k) = e/alpha(k-1)
                alpha(k) = a - beta(k)*c
            end do

            y(1) = b(1)
            do k=2,N
                y(k) = b(k) - beta(k)*y(k-1)
            end do
            
            x(N) = y(N)/alpha(N)
            do k=N-1,1,-1
                x(k) = (y(k)-c*x(k+1))/alpha(k)
            end do

        end subroutine Thomas

end subroutine newmark