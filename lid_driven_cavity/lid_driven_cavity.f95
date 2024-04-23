! computer code for lid-driven cavity
! Vahid Mollania Malakshah
! vahid.m.malakshah1998@gmail.com
module variables
    integer :: m ,n, i, j, k,mstep,kk
    Real(8), Dimension(:,:,:), Allocatable :: f, feq
    Real(8), Dimension(:,:), Allocatable :: u, v, rho, strf
    Real(8), Dimension(0:8) :: w,cx,cy
    Real(8) :: uo=0.10d0, sumvelo=0d0, rhoo=5d0, alpha=0.02d0, sum=0d0, usum=0d0, vsum=0d0
    Real(8) :: Re, omega, t1, t2, rhon, rhoav, rhom
end module variables

program cavity

    use variables
    implicit none
    m= 100; n= 100
    allocate(f(0:8,0:n,0:m),feq(0:8,0:n,0:m)) 
    allocate(u(0:n,0:m),v(0:n,0:m),rho(0:n,0:m),strf(0:n,0:m))
    
    open(2,file='uvfield')
    open(3,file='uvely')
    open(4,file='vvelx')
    open(8,file='timeu')
    !

    Re=uo*m/alpha
    print *, "Re=", Re
    omega=1d0/(3.*alpha+0.5)
    mstep=40000
    
    call Initial
    ! main loop
    do kk=1,mstep
        call collesion
        call streaming
        ! ——————————–
        call sfbound
        call rhouv
        print *, kk,u(0,m/2),v(0,m/2)!,rho(0,m/2),u(n,m/2),v(n,m/2),rho(n,m/2)
        write(8,*) kk,u(n/2,m/2),v(n/2,m/2)
    END DO
    ! end of the main loop
    call result
    stop
    ! end of the main program
end program cavity

subroutine Initial
    use  variables
    implicit none

    w(0)=4./9.
    do i=1,4
        w(i)=1./9.
    end do
    do i=5,8
        w(i)=1./36.
    end do

    cx(0)=0 ;cx(1)=1 ;cx(2)=0 ;cx(3)=-1 ;cx(4)=0 ;cx(5)=1 ;cx(6)=-1 ;cx(7)=-1 ;cx(8)=1
    cy(0)=0 ;cy(1)=0 ;cy(2)=1 ;cy(3)=0 ;cy(4)=-1 ;cy(5)=1 ;cy(6)=1 ;cy(7)=-1 ;cy(8)=-1

    do j=0,m
        do i=0,n
            rho(i,j)=rhoo
            u(i,j)=0.0
            v(i,j)=0.0
        end do
    end do

    do i=1,n-1
        u(i,m)=uo
    end do
    return
end

subroutine collesion
    use  variables
    implicit none
    do i = 0,n
        do j = 0,m
            t1 = u(i,j)*u(i,j)+v(i,j)*v(i,j)
            do k = 0,8
                t2 = u(i,j)*cx(k)+ v(i,j)*cy(k)
                feq(k,i,j) = rho(i,j)*w(k)*(1.0+3.0*t2+4.50*t2*t2-1.50*t1)
                f(k,i,j) = f(k,i,j)*(1. - omega) + feq(k,i,j)*omega
            end do
        end do
    end do
    return
end

subroutine streaming
    use  variables
    implicit none
    ! streaming
    do j = 1,m
        do i = n,1,-1
            f(1,i,j) = f(1,i-1,j)!RIGHT TO LEFT
        end do
        do i = 0,n-1
            f(3,i,j) = f(3,i+1,j)!LEFT TO RIGHT
        end do
    end do

    do i = 0,n
        do j = m,1,-1
            f(2,i,j) = f(2,i,j-1)!TOP TO BOTTOM
        end do
        do j =  0, m-1
            f(4,i,j) = f(4,i, j+1)
        end do
    end do

    do i = n,1,-1
        do j = m,1,-1
            f(5,i,j) = f(5,i-1,j-1)
        end do
    end do

    do i =0,n-1
        do j=m,1,-1
            f(6,i,j) = f(6,i+1,j-1)
        end do
    end do

    do i = 0,n-1
        do j = 0,m-1
            f(7,i,j) = f(7,i+1,j+1)
        end do
    end do

    do j = 0,m-1
        do i = n,1,-1
            f(8,i,j) = f(8,i-1,j+1)
        end do
    end do
return
end


subroutine sfbound
    use  variables
    implicit none
    do j = 0,m
    ! bounce back on west boundary
        f(1,0,j) = f(3,0,j)
        f(5,0,j) = f(7,0,j)
        f(8,0,j) = f(6,0,j)
        ! bounce back on east boundary
        f(3,n,j) = f(1,n,j)
        f(7,n,j) = f(5,n,j)
        f(6,n,j) = f(8,n,j)
    end do
    ! bounce back on south boundary
    do i = 0,n
        f(2,i,0) = f(4,i,0)
        f(5,i,0)=f(7,i,0)
        f(6,i,0)=f(8,i,0)
    end do

    ! moving lid, north boundary
    do i = 1,n-1
        rhon = f(0,i,m) + f(1,i,m) + f(3,i,m) + 2.*(f(2,i,m) + f(6,i,m) + f(5,i,m))
        f(4,i,m ) = f(2,i,m)
        f(7,i,m) = f(5,i,m)  - rhon*uo/6.0
        f(8,i,m) = f(6,i,m)  + rhon*uo/6.0
    end do
    return
    end

subroutine rhouv
    use  variables
    implicit none
    do i = 0,n
        do j = 0,m
            sum = 0
            do k = 0,8
                sum = sum + f(k,i,j)
            end do
            rho(i,j) = sum
        end do
    end do

    do i = 1,n
        do j = 1,m-1
            usum = 0
            vsum = 0
            do k = 0,8
                usum = usum + f(k,i,j)*cx(k)
                vsum = vsum + f(k,i,j)*cy(k)
            end do
        u(i,j) = usum/rho(i,j)
        v(i,j) = vsum / rho(i,j)
        end do
    end do
return
end

subroutine result
    use  variables
    implicit none
    open(5, file='streamf')
    ! streamfunction calculations
    strf(0,0)=0.
    do i=0,n
        rhoav=0.5*(rho(i-1,0)+rho(i,0))
        if(i.ne.0) strf(i,0)=strf(i-1,0)-rhoav*0.5*(v(i-1,0)+v(i,0))
        do j=1,m
            rhom=0.5*(rho(i,j)+rho(i,j-1))
            strf(i,j)=strf(i,j-1)+rhom*0.5*(u(i,j-1)+u(i,j))
        end do
    end do
    ! ———————————–
    write(2,*)"VARIABLES =X, Y, U, V, S"
    write(2,*)"ZONE ","I=",n+1,"J=",m+1,",","F=BLOCK"
    do j=0,m
        write(2,*)(i,i=0,n)
    end do
    do j=0,m
        write(2,*)(j,i=0,n)
    end do
    do j=0,m
        write(2,*)(u(i,j),i=0,n)
    end do
    do j=0,m
        write(2,*)(v(i,j),i=0,n)
    end do
    do j=0,m
        write(2,*)(strf(i,j),i=0,n)
    end do
    do j=0,m
        write(3,*)j/float(m),u(n/2,j)/uo,u(n/4,j)/uo,u(3*n/4,j)/uo
    end do
    do i=0,n
        write(4,*) i/float(n),v(i,m/2)/uo
    end do
return
end
!============end of the program
