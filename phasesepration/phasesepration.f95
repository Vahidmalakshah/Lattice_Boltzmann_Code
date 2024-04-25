
module variables
    integer :: m ,n, i, j, k, mstep, kk
    Real(8) :: t1, t2
    Real(8) :: zed, beta, ka ,sigma, mobi, tphi , taoH, taoL, rhoH, rhoL,gsumx, gsumy, gsum
    Real(8) :: phiin=0d0, Uin=0.005d0
    Real(8), Dimension(:,:,:), Allocatable :: f, say, heq, h,geq,g,gteq,gama
    Real(8), Dimension(:,:), Allocatable :: u, v, rho, uyin, phi, phi1 , gphix, gphiy, nx, ny 
    Real(8), Dimension(:,:), Allocatable :: lapciphi,p,tao,sv,fsx,fsy, mu
    Real(8), Dimension(0:8) :: w,cx,cy,wd
end module variables

program Kanal
    use variables
    implicit none
    m= 100; n= 100
    allocate(f(0:8,0:n,0:m),heq(0:8,0:n,0:m) ,h(0:8,0:n,0:m) , say(0:8,0:n,0:m),geq(0:8,0:n,0:m),g(0:8,0:n,0:m),gteq(0:8,0:n,0:m))
    allocate(gama(0:8,0:n,0:m))
    allocate(u(0:n,0:m),v(0:n,0:m),rho(0:n,0:m),uyin(0:n,0:m))
    allocate(phi(-1:n+1,-1:m+1), phi1(-1:n+1,-1:m+1), gphix(0:n,0:m),nx(0:n,0:m),ny(0:n,0:m))
    allocate( gphiy(0:n,0:m), lapciphi(0:n,0:m),tao(0:n,0:m),sv(0:n,0:m),fsx(0:n,0:m),fsy(0:n,0:m), mu(0:n,0:m), p(0:n,0:m))
    !!!!!!!!!!!!!
    open(2,file="uvfield3.plt")
    mstep= 300000

    ! main loop
    Call Initial
    Call uP
    Call Results

    do kk=1,mStep
        CALL COLLESION
        CALL steaming
        Call boundaryconditions
        CALL uP
        call Results
    end do
    
end program Kanal

Subroutine Initial
    Use variables
    Implicit None
    Real(8) :: z, rrr
    w(0)=4.0/9.0
    do i=1,4
        w(i)=1.0/9.0
    end do
    do i=5,8
        w(i)=1.0/36.0
    end do

    !!!!!!!!!!!
    cx(0)=0.0; cx(1)=1.0; cx(2)=0.0; cx(3)=-1.0; cx(4)=0.0; cx(5)=1.0; cx(6)=-1.0; cx(7)=-1.0; cx(8)=1.0
    cy(0)=0.0; cy(1)=0.0; cy(2)=1.0; cy(3)=0.0; cy(4)=-1.0; cy(5)=1.0; cy(6)=1.0; cy(7)=-1.0; cy(8)=-1.0

    zed=3.0; sigma=0.005d0
    beta=12.0*sigma/zed
    ka=3.0*sigma*zed/2.0
    rhoL=800d0; rhoH=1000d0; mobi=0.1
	tphi=3d0*mobi
    taoH=0.1d0; taoL=0.1d0

    Do i=0,n
			Do j=0,m
				
                Call RANDOM_NUMBER(rrr)
				phi(i,j) = 0.7d0 + 0.02*rrr
                
            End do
    end do

    u=0d0; v=0d0; p=0d0
     Do i = 0,n
        Do j=0,m
            Do k=0,8
                h(k,i,j)=phi(i,j)*w(k)
                g(k,i,j)=0d0
            end do
        end do
    end do
End subroutine

Subroutine  COLLESION
    use  variables
    implicit none
   
    Do i=0,n
        Do j=0,m
            t1=u(i,j)*u(i,j)+v(i,j)*v(i,j)
            Do k=0,8

                t2=u(i,j)*cx(k)+v(i,j)*cy(k)

                say(k,i,j)=w(k)*(1.0+3.0*(t2)+4.50*(t2)**2-1.50*(t1))

                heq(k,i,j)=phi(i,j)*say(k,i,j)+3.0*w(k)*mobi*((4.0/zed)*phi(i,j)*(1-phi(i,j)))*(cx(k)*nx(i,j)+cy(k)*ny(i,j))

                f(k,i,j) = ((say(k,i,j)-w(k))*((rhoH-rhol)/3.0) + say(k,i,j)*mu(i,j)) * &
                           (((cx(k)-u(i,j))*gphix(i,j)) + ((cy(k)-v(i,j))*gphiy(i,j)))

                geq(k,i,j)=p(i,j)*w(k)+rho(i,j)*(say(k,i,j)-w(k))/3d0
                gteq(k,i,j)=geq(k,i,j)-0.50*f(k,i,j)
                gama(k,i,j)=-sv(i,j)*(g(k,i,j)-gteq(k,i,j))

                g(k,i,j)=f(k,i,j)+gama(k,i,j)+g(k,i,j)
                h(k,i,j)=h(k,i,j)-((h(k,i,j)-heq(k,i,j))/(tphi+0.50))

            End Do

        End Do
    End Do
   
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
Subroutine  steaming
    use  variables
    implicit none
   
    do j = 0,m
        do i = n,1,-1
            g(1,i,j) = g(1,i-1,j) ! Right to left
            h(1,i,j) = h(1,i-1,j)
        end do
        do i = 0,n-1
            g(3,i,j) = g(3,i+1,j) ! Left to right
            h(3,i,j) = h(3,i+1,j)
        end do
    end do
    do i = 0,n
        do j = m,1,-1
            g(2,i,j) = g(2,i,j-1) ! Top to bottom
            h(2,i,j) = h(2,i,j-1)
        end do
        do j =  0, m-1
            g(4,i,j) = g(4,i, j+1) ! Bottom to top
            h(4,i,j) = h(4,i, j+1)
        end do
    end do

    do i = n,1,-1
        do j = m,1,-1
            g(5,i,j) = g(5,i-1,j-1)
            h(5,i,j) = h(5,i-1,j-1)
        end do
    end do
    do i =0,n-1
        do j=m,1,-1
            g(6,i,j) = g(6,i+1,j-1)
            h(6,i,j) = h(6,i+1,j-1)
        end do
    end do

    do i = 0,n-1
        do j = 0,m-1
            g(7,i,j) = g(7,i+1,j+1)
            h(7,i,j) = h(7,i+1,j+1)
        end do
    end do

    do j = 0,m-1
        do i = n,1,-1
            g(8,i,j) = g(8,i-1,j+1)
            h(8,i,j) = h(8,i-1,j+1)
        end do
    end do
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine boundaryconditions
    use  variables
    implicit none

    do j=0,m
       !!!west boundry
       h(5,0,j)=h(5,n,j)
       g(5,0,j)=g(5,n,j)
       h(1,0,j)=h(1,n,j)
       g(1,0,j)=g(1,n,j)
       h(8,0,j)=h(8,n,j)
       g(8,0,j)=g(8,n,j)
       !!!!east boundry
       h(6,n,j)=h(6,0,j)
       g(6,n,j)=g(6,0,j)
       h(3,n,j)=h(3,0,j)
       g(3,n,j)=g(3,0,j)
       h(7,n,j)=h(7,0,j)
       g(7,n,j)=g(7,0,j)
    end do
    do i=0,n
        !north boundary
        h(7,i,m)=h(7,i,0)
        h(8,i,m)=h(8,i,0)
        h(4,i,m)=h(4,i,0)
        g(7,i,m)=g(7,i,0)
        g(8,i,m)=g(8,i,0)
        g(4,i,m)=g(4,i,0)

        !south boundary
        h(5,i,0)= h(5,i,m)
        h(6,i,0)=h(6,i,m)
        h(2,i,0)=h(2,i,m)
        g(5,i,0)= g(5,i,m)
        g(6,i,0)=g(6,i,m)
        g(2,i,0)=g(2,i,m)
    end do
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine  uP
    use  variables
    implicit none
    Real(8) :: temp1

    do j=0,m
        do i=0,n

            phi(i,j)=sum(h(:,i,j),1)
            rho(i,j)=rhoL+phi(i,j)*(rhoH-rhoL)
        end do
    end do
    call moshtagh
    do i=0,n
        do j=0,m
            temp1=sqrt(gphix(i,j)**2+gphiy(i,j)**2) + 1.d-32
            tao(i,j)=1d0/ ( (1d0/taoL) + phi(i,j)*( (1d0/taoH) - (1d0/taoL) ) )
            sv(i,j)=1.0/(tao(i,j)+0.50)   
            nx(i,j)=gphix(i,j)/temp1
            ny(i,j)=gphiy(i,j)/temp1
            mu(i,j)=4.0*beta*phi(i,j)*(phi(i,j)-1.0)*(phi(i,j)-0.50)-ka*lapciphi(i,j)
            fsx(i,j)=mu(i,j)*gphix(i,j)
            fsy(i,j)=mu(i,j)*gphiy(i,j)
            u(i,j)=(3.0/rho(i,j))*(g(1,i,j)+g(5,i,j)+g(8,i,j)-g(3,i,j)-g(7,i,j)-g(6,i,j))+(1.0/(2.0*rho(i,j)))*fsx(i,j)
            v(i,j)=(3.0/rho(i,j))*(g(2,i,j)+g(5,i,j)+g(6,i,j)-g(4,i,j)-g(7,i,j)-g(8,i,j))+(1.0/(2.0*rho(i,j)))*fsy(i,j)
            p(i,j)=sum(g(:,i,j),1)+(rhoH-rhoL)*(u(i,j)*gphix(i,j)+v(i,j)*gphiy(i,j))/6d0

        end do
    end do

End Subroutine

subroutine moshtagh
    use  variables
    implicit none

    phi1(0:n,0:m) = phi(0:n,0:m)
    phi1(-1,:) = phi(0,:)
    phi1(n+1,:) = phi(n,:)
    phi1(:,-1) = phi(:,0)
    phi1(:,m+1) = phi(:,m)

    do i=0,n
        do j=0,m
            gsumx=0.0; gsumy=0.0; gsum=0.0
            do k=0,8
				gsumx = gsumx + cx(k) * w(k) * phi1(INT(i + cx(k)), INT(j + cy(k)))
				gsumy = gsumy + cy(k) * w(k) * phi1(i + INT(cx(k)), j + INT(cy(k)))
				gsum  = gsum + w(k) * (phi1(i+INT(cx(k)),j+INT(cy(k)))-phi1(i,j))

            end do
            gphix(i,j)=3*gsumx
            gphiy(i,j)=3*gsumy
            lapciphi(i,j)=6.0*gsum
        end do
    end do

end subroutine

!!!!!!!!!!!!1

Subroutine Results
    Use variables
    Implicit None

     write(*,*) kk,sum(phi)!, u(0,50)/Uin

    if (MOD(kk,300)==0) then

        WRITE(2,*) 'VARIABLES = "X" , "Y" , "U" , "V" , "Phi" '
        WRITE(2,*) 'ZONE T="Rectangular zone"'
        WRITE(2,*) ' STRANDID=0, SOLUTIONTIME=0'
        WRITE(2,*) ' I=',n+1,',J=',m+1,',K=1, ZONETYPE=Ordered'
        WRITE(2,*) ' DATAPACKING=POINT'
        WRITE(2,*) ' DT=(SINGLE SINGLE SINGLE SINGLE SINGLE )'
        Do j=0,m
            DO i=0,n

               Write(2,*) i+1, j+1, u(i,j), v(i,j) ,phi(i,j)

            End do
        End do

    End if

End subroutine

!============end of the program
    
    

    
    
