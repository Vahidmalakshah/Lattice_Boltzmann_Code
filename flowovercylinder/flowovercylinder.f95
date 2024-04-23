! computer code for flow over a cylinder
! Vahid Mollania Malakshah
! vahid.m.malakshah1998@gmail.com
module variables
    integer :: m ,n, i, j, k,mstep,kk
    Real(8), Dimension(:,:,:), Allocatable :: f, feq, DELTAK, DELTAKD
    Real(8), Dimension(:,:), Allocatable :: u, v, rho, strf, sf, xbound, dfydx, dfxdy
    Real(8), Dimension(0:8) :: w,cx,cy, wd
    Real(8) :: uo=0.10d0, sumvelo=0d0, rhoo=5d0, alpha=0.01d0, sum=0d0, usum=0d0, vsum=0d0 
    Real(8) :: radius = 10d0, xc = 100d0, yc = 50d0, Rc = 25d0
    Real(8) :: Re, omega, t1, t2, rhow, rhoav, rhom
end module variables

program flowovercylinder

    use variables
    implicit none
    m= 100; n= 500
    allocate(f(0:8,0:n,0:m),feq(0:8,0:n,0:m), DELTAK(0:8,0:n,0:m),DELTAKD(0:8,0:n,0:m)) 
    allocate(u(0:n,0:m),v(0:n,0:m),rho(0:n,0:m),strf(0:n,0:m), sf(0:n,0:m), xbound(0:n,0:m))
    open(2,file='uvfield')
    open(3,file='uvely')
    open(4,file='vvelx')
    open(8,file='timeu')
    !
    uo=0.20 ; sumvelo=0.0 ; rhoo=5.00 ; alpha=0.003
    Re = uo*m/alpha
    print *,"Re = ",Re
    omega= 0.56
    !1d0/(3.*alpha+0.5)
    print *,"omega = ",omega
    mstep= 30000
    
    call Initial

    ! main loop
    do kk=1,mstep
      call collesion
      call streaming
      ! ——————————–
      call sfbound
      call rhouv
      print *, kk, u(0,m/2),v(0,m/2),rho(0,m/2),u(n,m/2),v(n,m/2),rho(n,m/2)
      write(8,*) kk,u(n/2,m/2),v(n/2,m/2)
    END DO
    ! end of the main loop
    call result
stop
end program flowovercylinder
! end of the main program

subroutine Initial
    use  variables
    implicit none
    w(0) = 4./9.
    do i = 1,4
        w(i) = 1./9.
    end do
    do i = 5,8
        w(i) = 1./36.
    end do

    cx(0)=0 ;cx(1)=1 ;cx(2)=0 ;cx(3)=-1 ;cx(4)=0 ;cx(5)=1 ;cx(6)=-1 ;cx(7)=-1 ;cx(8)=1
    cy(0)=0 ;cy(1)=0 ;cy(2)=1 ;cy(3)=0 ;cy(4)=-1 ;cy(5)=1 ;cy(6)=1 ;cy(7)=-1 ;cy(8)=-1
    wd(0) = 0;wd(1) = 3;wd(2) = 4;wd(3) = 1;wd(4) = 2;wd(5) = 7;wd(6) = 8;wd(7) = 5;wd(8) = 6

    ! Initialize arrays
    sf(0:n,0:m) = 0.0
    xbound(0:n,0:m) = 0.0
    DELTAK(0:8,0:n,0:m) = 0.0
    DELTAKD(0:8,0:n,0:m) = 0.0

    ! Set up solid boundary
    do i = 0,n
        do j = 0,m
            if ((sqrt((i-xc)**2. +(j-yc)**2.)) < Rc)  then
                sf(i,j) = 1.0
            end if
        end do
    end do

    do i = 1,n-1
        do j = 1,m-1
            if (sf(i,j) == 0.0) then
                do k = 0,8
                    if (sf(i+int(cx(k)),j+int(cy(k))) == 1.0) then
                        xbound(i,j) = 1.
                    end if
                end do
            end if
        end do
    end do

    do i = 0,n
        do j = 0,m
            if (xbound(i,j) == 1.0) then
                do k = 0,8
                    if (sf(i+int(cx(k)),j+int(cy(k))) == 1.0) then
                        DELTAK(k,i,j) = k
                        DELTAKD(int(wd(k)),i,j) = wd(k)
                    end if
                end do
            end if
        end do
    end do

    ! Initialize flow field
    do i = 0,n
        do j = 0,m
            rho(i,j) = rhoo
            v(i,j) = 0.0
            u(i,j) = 0.0
        end do
    end do

    do j = 1,m-1
        u(0,j) = uo
    end do

    return
end

subroutine collesion
    use variables
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
end subroutine

subroutine streaming
    use variables
    implicit none
    ! Streaming
    do j = 0,m
        do i = n,1,-1
            f(1,i,j) = f(1,i-1,j) ! Right to left
        end do
        do i = 0,n-1
            f(3,i,j) = f(3,i+1,j) ! Left to right
        end do
    end do
    do i = 0,n
        do j = m,1,-1
            f(2,i,j) = f(2,i,j-1) ! Top to bottom
        end do
        do j =  0, m-1
            f(4,i,j) = f(4,i, j+1) ! Bottom to top
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
end subroutine

subroutine sfbound
    use variables
    implicit none
    do j = 0,m
        ! West boundary/ Inlet velocity
        rhow = (f(0,0,j) + f(2,0,j) + f(4,0,j) + 2.*(f(3,0,j)+f(6,0,j)+f(7,0,j) ))/(1.0-uo)
        f(1,0,j) = f(3,0,j) + 2.*rhow*uo/3.
        f(5,0,j) = f(7,0,j) + rhow*uo/6.
        f(8,0,j) = f(6,0,j) + rhow*uo/6.
    end do

    ! Bounce back on south boundary
    do i = 0,n
        f(2,i,0) = f(4,i,0)
        f(5,i,0) = f(7,i,0)
        f(6,i,0) = f(8,i,0)
    end do

    ! Bounce back on north boundary
    do i = 0,n
        f(4,i,m) = f(2,i,m)
        f(7,i,m) = f(5,i,m)  
        f(8,i,m) = f(6,i,m)  
    end do

    ! Bounce back on cylinder
    do i = 1,n-1
        do j = 1,m-1
            if (xbound(i,j) == 1.0) then
                do k = 0,8
                    if (sf(i+int(cx(k)),j+int(cy(k))) == 1.0) then
                        f(int(DELTAKD(k,i,j)),i,j) = f(int(DELTAK(k,i,j)),i,j)
                    end if
                end do
            end if
        end do
    end do

    ! Open boundary condition at the outlet
    do j = 1,m
        f(3,n,j) = f(3,n-1,j) 
        f(6,n,j) = f(6,n-1,j) 
        f(7,n,j) = f(7,n-1,j) 
    end do
    return
end subroutine

subroutine rhouv
    use variables
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
    do j = 1,m
        v(n,j) = 0.0
    end do
    do i = 0,n
        do j = 0,m
            if (sf(i,j) == 1.0)  then
                u(i,j) = 0.0
                v(i,j) = 0.0
            end if
        end do
    end do
    return
end subroutine

subroutine result
  use variables
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
    write(3,*)j/float(m),u(n/10,j)/uo,u(n/4,j)/uo,u(n/2,j)/uo,u(3*n/4,j)/uo
  end do
  do i=0,n
    write(4,*) i/float(n),v(i,m/2)/uo
  end do
return
end
!============end of the program