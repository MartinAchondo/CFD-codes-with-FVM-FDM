
subroutine theta_method(theta,y0,h,x,tf,t,y1,y2,N,alpha,beta,gamma) bind(C, name="theta_method")
    use iso_c_binding
    implicit none

    integer(c_int), intent(in), value :: N
    real(c_double), intent(in), value :: theta,h,x,alpha,beta,gamma,tf
    real(c_double), intent(in) :: y0(2)
    real(c_double), intent(out) :: y1(N),y2(N)
    real(c_double), intent(inout) :: t(N)

    real(c_double) :: y_old(2),y(2),v1,z,tx
    integer(c_int), parameter :: imax=5000
    integer(c_int) :: iter

    iter = 2
    tx = t(1)
    v1 = h*theta
    z = 1.0d0/(alpha*(v1)**2+beta*v1 +1)

    y_old(:) = y0(:)

    do while(tx.lt.tf)

        y(1) = z*((beta*v1+1)*y_old(1)+v1*y_old(2) + h*(1-theta)*(-alpha*v1*y_old(1)+y_old(2)) + h*gamma*v1*Pt_f(x,tx,theta,h))

        y(2) = z*(-alpha*v1*y_old(1)+y_old(2) + h*(1-theta)*(-alpha*y_old(1)-(alpha*v1+beta)*y_old(2))+ h*gamma*Pt_f(x,tx,theta,h))

        y1(iter) = y(1)
        y2(iter) = y(2)
        t(iter) = tx + h

        y_old(:) = y(:)

        tx = tx + h
        iter = iter + 1
        write(*,*) iter,N

        if (iter.gt.N) then
            exit
        end if

    end do 

    contains 

        real(c_double) function Pt_f(x,t,theta,h)
            implicit none
            real(c_double), intent(in) :: x,t,theta,h

            Pt_f = theta*Pt(x,t+h) + (1-theta)*Pt(x,t)

        end function Pt_f

        real(c_double) function Pt(x,t)
            implicit none

            real(c_double), intent(in) :: x,t
            real(c_double) :: b,a,dP,omega_0
            real(c_double), parameter :: pi = 4.0d0*atan(1.0d0)

            b = 133.32
            a = 10*b
            dP = 0.25*b
            omega_0 = 2.0d0*pi/0.8

            Pt = x*dP*(a+b*cos(omega_0*t))

        end function Pt

end subroutine theta_method