    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IBM %%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine IBM(problem, imax, jmax, ppu, ppv, ppp, ibound, t, xg, yg, dx, dy, N, x, y, sc, nxm, nym, uflag, vflag, pflag, &
                    ua, va, ivp, jvp, ux, vx, uy, vy, xc, yc, invagu, invagv, invagp, invaguin, invagvin, gFIu, gFIv, gFIp, &
                    gFIuin, gFIvin, indu, indv, indp, induin, indvin, ghu, ghv, ghp, ghuin, ghvin, IMu, IMv, IMp, P, vos)
    use interfaces2
    use defs
    implicit none
    character (LEN=30), INTENT(IN) :: problem
    integer, intent(IN) :: imax, jmax, N
    real(RP), intent(IN):: t
    real(RP), dimension(0:), intent(IN) :: xg, yg, dx, dy    
    real(DP), dimension(:), intent(IN) :: x, y, sc
    real(DP), dimension(:), intent(OUT) :: nxm, nym
    integer, intent(OUT) :: ibound, ppu, ppv, ppp, ivp, jvp
    integer(I2B), dimension(0:,0:), intent(INOUT) :: uflag, vflag, pflag
    real(DP), dimension(:,:,:), allocatable, intent(OUT) :: invagu, invagv, invagp, invaguin, invagvin
    integer, dimension(:,:,:), allocatable, intent(OUT) :: gFIu, gFIv, gFIp, gFIuin, gFIvin
    integer, dimension(:,:), allocatable, intent(OUT) :: indu, indv, indp, induin, indvin
    real(DP), dimension(:,:), allocatable, intent(OUT) :: ghu, ghv, ghp, ghuin, ghvin
    integer, dimension(:), allocatable, intent(OUT) :: IMu, IMv, IMp
    integer, dimension(:,:), intent(OUT) :: ua, va
    real(DP), dimension(0:), intent(OUT) :: ux, vx, uy, vy, xc, yc
    real(DP), intent(IN) :: P
    real(DP), dimension(0:,0:), intent(INOUT) :: vos
    
    ! local Variables
    real(DP), dimension(N) :: ax, bx, cx, ay, by, cy
    real(DP) upacket(9*N,7), vpacket(9*N,7), ppacket(9*N,7)
    real(DP) a, b, s, xd, yd, dist, distmin, fs, fsd, sn, s0, res, xn, yn, nxx, nyy, lx, ly, delta, a0, an, &
             scip1, sci, scim1, xim1, xip1, yi, yip1, yim1, dist1, dist2, ratio, xnn, ynn
    real(DP) dxx, dyy, xi, yj, uv, wv, x1, x2, y1, y2, eps, xx, yy, xim, yim, h1, h2, xvpp, &
            yvpp, dh, teta, ss1, ss2, ss, root, ds, axx, ayy, bxx, byy, cxx, cyy, nxxp, nyyp
    integer i, j, k, ib, ic, jc, m, pl, iter, pp, ii, jj , ip1, im1, iii, jjj, i1, j1, i2, j2
    logical found
    integer, dimension(0:imax+1,0:jmax+1) :: flagu, flagv, flagp
    real(DP), dimension(3,3) :: r
    integer, dimension(imax,jmax) :: pa
    
    OPEN(11,FILE='grid.plt',STATUS='UNKNOWN')
    OPEN(13,FILE='inactiveu.plt',STATUS='UNKNOWN')
    OPEN(14,FILE='inactivev.plt',STATUS='UNKNOWN')
    OPEN(15,FILE='inactiveuvir.plt',STATUS='UNKNOWN')
    OPEN(16,FILE='inactivevvir.plt',STATUS='UNKNOWN')
    OPEN(17,FILE='flags.plt',STATUS='UNKNOWN')
    OPEN(18,FILE='ghostu.plt',STATUS='UNKNOWN')
    OPEN(19,FILE='ghostv.plt',STATUS='UNKNOWN')
    OPEN(20,FILE='ghostp.plt',STATUS='UNKNOWN')
    OPEN(21,FILE='ghostviru.plt',STATUS='UNKNOWN')
    OPEN(22,FILE='ghostvirum.plt',STATUS='UNKNOWN')
    OPEN(23,FILE='ghostvirv.plt',STATUS='UNKNOWN')
    OPEN(24,FILE='ghostvirvm.plt',STATUS='UNKNOWN')
    OPEN(25,FILE='ghostvirp.plt',STATUS='UNKNOWN')
    OPEN(26,FILE='ghostvirpm.plt',STATUS='UNKNOWN')
    OPEN(27,FILE='inactiveuvirm.plt',STATUS='UNKNOWN')
    OPEN(28,FILE='inactivevvirm.plt',STATUS='UNKNOWN')
    OPEN(29,FILE='activeu.plt',STATUS='UNKNOWN')
    OPEN(30,FILE='activev.plt',STATUS='UNKNOWN')
    OPEN(31,FILE='imghostu.plt',STATUS='UNKNOWN')
    OPEN(32,FILE='imghostv.plt',STATUS='UNKNOWN')
    OPEN(33,FILE='imghostp.plt',STATUS='UNKNOWN')
    OPEN(37,FILE='interfacep.plt',STATUS='UNKNOWN')
    OPEN(38,FILE='interfaceu.plt',STATUS='UNKNOWN')
    OPEN(39,FILE='interfacev.plt',STATUS='UNKNOWN')
    
    eps=1.0d-10
    if (.NOT.(immersed)) goto 150
    !
    !compute coefficients of a polynomial. Use forward points for the initial point
    !
    do i=1, N
        if (i==1) then      !use forward points
            im1=i
            ii=i+1
            ip1=i+2
        else if(i==N) then  !use backward points
            ip1=i
            ii=i-1
            im1=i-2
        else                
            ip1=i+1
            ii=i
            im1=i-1
        end if
        scip1=sc(ip1)
        sci=sc(ii)
        scim1=sc(im1)
        xi=x(ii)
        xim1=x(im1)
        xip1=x(ip1)
        yi=y(ii)
        yim1=y(im1)
        yip1=y(ip1)
        !coefficients for x polynomial        
        ax(i)=(scip1*(xi-xim1)+sci*(xim1-xip1)+scim1*(xip1-xi))/&
                (scim1-sci)/(scim1-scip1)/(sci-scip1)
        bx(i)=(scip1**2*(xim1-xi)+scim1**2*(xi-xip1)+&
                sci**2*(xip1-xim1))/(scim1-sci)/(scim1-scip1)/(sci-scip1)
        cx(i)=(scim1*scip1*(scip1-scim1)*xi+sci**2*(scip1*xim1-scim1*xip1)+&
                sci*(xip1*scim1**2-xim1*scip1**2))/(scim1-sci)/(scim1-scip1)/(sci-scip1)
        !coefficients for y polyno1ial
        ay(i)=(scip1*(yi-yim1)+sci*(yim1-yip1)+scim1*(yip1-yi))/&
                (scim1-sci)/(scim1-scip1)/(sci-scip1)
        by(i)=(scip1**2*(yim1-yi)+scim1**2*(yi-yip1)+sci**2*&
                (yip1-yim1))/(scim1-sci)/(scim1-scip1)/(sci-scip1)
        cy(i)=(scim1*scip1*(scip1-scim1)*yi+sci**2*(scip1*yim1-scim1*yip1)+&
                sci*(yip1*scim1**2-yim1*scip1**2))/(scim1-sci)/(scim1-scip1)/(sci-scip1)
    end do
    !compute nx and ny normals at marker points. These values will be used later to determine the VOS values at immersed boundary cells.
    do i=1, N
        sn=sc(i)
        xn=ax(i)*sn**2+bx(i)*sn+cx(i)     
        yn=ay(i)*sn**2+by(i)*sn+cy(i)     
        xd=2.0d0*ax(i)*sn+bx(i)           !d(xsn)/ds
        yd=2.0d0*ay(i)*sn+by(i)           !d(ysn)/ds
        !compute x and y normals at sn 
        nxm(i)=-yd/dsqrt(xd**2+yd**2)       !nx
        nym(i)=xd/dsqrt(xd**2+yd**2)        !ny
    end do
    
        
    !compute x and y coordinates of u, v and p points
    dxx=dble(dx(0))
    dyy=dble(dy(0))
    xc(0)=-dble(dx(0))/2.0d0        !x coordinate of cell center
    ux(0)=dble(xg(0))               !x coordinate of u-velocity point
    vx(0)=xc(0)                     !x coordinate of v-velocity point     
    do i=1, imax
        xc(i)=xc(i-1)+dble(dx(i-1)+dx(i))/2.0d0
        ux(i)=dble(xg(i))
        vx(i)=xc(i)
    end do
    xc(imax+1)=xc(imax)+dble(dx(imax)+dx(imax+1))/2.0d0
    ux(imax+1)=ux(imax)+dx(imax)
    vx(imax+1)=xc(imax+1)
    yc(0)=-dble(dy(0))/2.0d0        !y coordinate of cell center
    uy(0)=yc(0)                     !y coordinate of u-velocity point
    vy(0)=dble(yg(0))               !y coordinate of v-velocity point
    do j=1, jmax
        yc(j)=yc(j-1)+dble(dy(j-1)+dy(j))/2.0d0
        uy(j)=yc(j)
        vy(j)=dble(yg(j))
    end do
    yc(jmax+1)=yc(jmax)+dble(dy(jmax)+dy(jmax+1))/2.0d0
    uy(jmax+1)=yc(jmax+1)
    vy(jmax+1)=vy(jmax)+dy(jmax)
    !set all cell as fluid cells initially for each variable on a staggered grid
    do i=2, imax-1
        do j=2, jmax-1
            flagu(i,j)=0
            flagv(i,j)=0
            flagp(i,j)=0
        end do
    end do
    pp=0
    !search ghost points for u-velocity
    do i=1, N
        !
        !Find the closest grid point (ic,jc) of u-velocity to a marker point (ib)
        !
        distmin=1.0d100
        do j=1, imax
            do k=1, jmax
                dist=dsqrt((x(i)-ux(j))**2+(y(i)-uy(k))**2)
                if (dist<distmin) then 
                    distmin=dist
                    ic=j    !i index of the closest grid point 
                    jc=k    !j index o the closest grid point
                end if
            end do
        end do 
        !do while loop for the packet of points
        do j=ic-1, ic+1
            do k=jc-1, jc+1
                !Use Newton-Raphson algorithm to find the sn (normal from point to curve) 
                pp=pp+1
                s0=sc(i)
                iter=0
30              fs=(2.0d0*ay(i)**2+2.0d0*ax(i)**2)*s0**3+(3.0d0*ax(i)*bx(i)+3.0d0*ay(i)*by(i))*s0**2+(2.0d0*ay(i)*cy(i)+&
                    2.0d0*ax(i)*cx(i)+bx(i)**2+by(i)**2 -&
                    2.0d0*ay(i)*uy(k)-2.0d0*ax(i)*ux(j))*s0+(cx(i)*bx(i)+by(i)*cy(i)-ux(j)*bx(i)-by(i)*uy(k))
                fsd=3.0d0*(2.0d0*ay(i)**2+2.0d0*ax(i)**2)*s0**2+2.0d0*(3.0d0*ax(i)*bx(i)+3.0d0*ay(i)*by(i))*s0+&
                    (2.0d0*ay(i)*cy(i)+2.0d0*ax(i)*cx(i)+bx(i)**2+by(i)**2 -&
                    2.0d0*ay(i)*uy(k)-2.0d0*ax(i)*ux(j))
                sn=s0-fs/fsd
                res=(2.0d0*ay(i)**2+2.0d0*ax(i)**2)*sn**3+(3.0d0*ax(i)*bx(i)+3.0d0*ay(i)*by(i))*sn**2+(2.0d0*ay(i)*cy(i)+&
                    2.0d0*ax(i)*cx(i)+bx(i)**2+by(i)**2 -&
                    2.0d0*ay(i)*uy(k)-2.0d0*ax(i)*ux(j))*sn+(cx(i)*bx(i)+by(i)*cy(i)-ux(j)*bx(i)-by(i)*uy(k))
                iter=iter+1
                if (iter>1000) goto 40
                if (dabs(res)<eps) then 
                    goto 40
                else 
                    s0=sn
                    goto 30
                end if 
40              continue
                xn=ax(i)*sn**2+bx(i)*sn+cx(i)     
                yn=ay(i)*sn**2+by(i)*sn+cy(i)     
                xd=2.0d0*ax(i)*sn+bx(i)           !d(xsn)/ds
                yd=2.0d0*ay(i)*sn+by(i)           !d(ysn)/ds
                !compute x and y normals at sn 
                nxx=-yd/dsqrt(xd**2+yd**2)       !nx
                nyy=xd/dsqrt(xd**2+yd**2)        !ny
                !backup i ndex, j index, xn and yn values of packet points inorder to find virtual point later
                upacket(pp,1)=j         !i index
                upacket(pp,2)=k         !j index
                upacket(pp,3)=xn        
                upacket(pp,4)=yn
                upacket(pp,5)=nxx
                upacket(pp,6)=nyy
                upacket(pp,7)=sn
                !compute vector lambda
                lx=(xn-ux(j))/dsqrt((xn-ux(j))**2+(yn-uy(k))**2)     !Lambda_x
                ly=(yn-uy(k))/dsqrt((xn-ux(j))**2+(yn-uy(k))**2)     !Lambda_y
                delta=nxx*lx+nyy*ly
                if (delta>0.0d0) then     
                    flagu(j,k)=1             !(j,k) point lies inside the object (solid phase)
                else
                    flagu(j,k)=0             ! (j,k) point lies outside the object (fluid phase)
                end if
            end do           
        end do !end of packet of points loop
    end do
    !search ib points for v-cell
    pp=0
    do i=1, N
        distmin=1.0d100
        do j=1, imax
            do k=1, jmax
                dist=dsqrt((x(i)-vx(j))**2+(y(i)-vy(k))**2)
                if (dist<distmin) then 
                    distmin=dist
                    ic=j    !i index of the closest grid point 
                    jc=k    !j index o the closest grid point
                end if
            end do
        end do 
        !do while loop for the packet of points
        do j=ic-1, ic+1
            do k=jc-1, jc+1
                !Use Newton-Raphson algorithm to find the sn (normal from point to curve) 
                pp=pp+1
                s0=sc(i)
                iter=0
50              fs=(2.0d0*ay(i)**2+2.0d0*ax(i)**2)*s0**3+(3.0d0*ax(i)*bx(i)+3.0d0*ay(i)*by(i))*s0**2+&
                    (2.0d0*ay(i)*cy(i)+2.0d0*ax(i)*cx(i)+bx(i)**2+by(i)**2 -&
                    2.0d0*ay(i)*vy(k)-2.0d0*ax(i)*vx(j))*s0+(cx(i)*bx(i)+by(i)*cy(i)-vx(j)*bx(i)-by(i)*vy(k))
                fsd=3.0d0*(2.0d0*ay(i)**2+2.0d0*ax(i)**2)*s0**2+2.0d0*(3.0d0*ax(i)*bx(i)+3.0d0*ay(i)*by(i))*s0+&
                    (2.0d0*ay(i)*cy(i)+2.0d0*ax(i)*cx(i)+bx(i)**2+by(i)**2 -&
                    2.0d0*ay(i)*vy(k)-2.0d0*ax(i)*vx(j))
                sn=s0-fs/fsd
                res=(2.0d0*ay(i)**2+2.0d0*ax(i)**2)*sn**3+(3.0d0*ax(i)*bx(i)+3.0d0*ay(i)*by(i))*sn**2+&
                    (2.0d0*ay(i)*cy(i)+2.0d0*ax(i)*cx(i)+bx(i)**2+by(i)**2 -&
                    2.0d0*ay(i)*vy(k)-2.0d0*ax(i)*vx(j))*sn+(cx(i)*bx(i)+by(i)*cy(i)-vx(j)*bx(i)-by(i)*vy(k))
                iter=iter+1
                if (iter>1000) goto 60
                if (dabs(res)<eps) then 
                    goto 60
                else 
                    s0=sn
                    goto 50
                end if 
60              continue
                xn=ax(i)*sn**2+bx(i)*sn+cx(i)
                yn=ay(i)*sn**2+by(i)*sn+cy(i)
                xd=2.0d0*ax(i)*sn+bx(i)           !d(xsn)/ds
                yd=2.0d0*ay(i)*sn+by(i)           !d(ysn)/ds
                !compute x and y normals at sn 
                nxx=-yd/dsqrt(xd**2+yd**2)       !nx
                nyy=xd/dsqrt(xd**2+yd**2)        !ny
                !backup i ndex, j index, xn and yn values of packet points inorder to find virtual point later
                vpacket(pp,1)=j
                vpacket(pp,2)=k
                vpacket(pp,3)=xn
                vpacket(pp,4)=yn
                vpacket(pp,5)=nxx
                vpacket(pp,6)=nyy
                vpacket(pp,7)=sn
                !compute vector lambda
                lx=(xn-vx(j))/dsqrt((xn-vx(j))**2+(yn-vy(k))**2)     !Lambda_x
                ly=(yn-vy(k))/dsqrt((xn-vx(j))**2+(yn-vy(k))**2)     !Lambda_y
                delta=nxx*lx+nyy*ly
                if (delta>0.0d0) then     
                    flagv(j,k)=1             !(j,k) point lies inside the object (solid phase)
                else
                    flagv(j,k)=0            ! (j,k) point lies outside the object (fluid phase)
                end if
            end do           
        end do !end of packet of points loop
    end do     !end of ib point search
    !search ib points for p-cell
    pp=0
    do i=1, N
        distmin=1.0d100
        do j=1, imax
            do k=1, jmax
                dist=dsqrt((x(i)-xc(j))**2+(y(i)-yc(k))**2)
                if (dist<distmin) then 
                    distmin=dist
                    ic=j    !i index of the closest grid point 
                    jc=k    !j index o the closest grid point
                end if
            end do
        end do 
        !do while loop for the packet of points
        do j=ic-1, ic+1
            do k=jc-1, jc+1
                !Use Newton-Raphson algorithm to find the sn (normal from point to curve)
                pp=pp+1
                s0=sc(i)
                iter=0
70              fs=(2.0d0*ay(i)**2+2.d0*ax(i)**2)*s0**3+(3.0d0*ax(i)*bx(i)+3.0d0*ay(i)*by(i))*s0**2+&
                    (2.0d0*ay(i)*cy(i)+2.0d0*ax(i)*cx(i)+bx(i)**2+by(i)**2 -&
                    2.0d0*ay(i)*yc(k)-2.0d0*ax(i)*xc(j))*s0+(cx(i)*bx(i)+by(i)*cy(i)-xc(j)*bx(i)-by(i)*yc(k))
                fsd=3.0d0*(2.0d0*ay(i)**2+2.0d0*ax(i)**2)*s0**2+2.0d0*(3.0d0*ax(i)*bx(i)+3.0d0*ay(i)*by(i))*s0+&
                    (2.0d0*ay(i)*cy(i)+2.0d0*ax(i)*cx(i)+bx(i)**2+by(i)**2 -&
                    2.0d0*ay(i)*yc(k)-2.0d0*ax(i)*xc(j))
                sn=s0-fs/fsd
                res=(2.0d0*ay(i)**2+2.0d0*ax(i)**2)*sn**3+(3.0d0*ax(i)*bx(i)+3.0d0*ay(i)*by(i))*sn**2+&
                        (2.0d0*ay(i)*cy(i)+2.0d0*ax(i)*cx(i)+bx(i)**2+by(i)**2 -&
                        2.0d0*ay(i)*yc(k)-2.0d0*ax(i)*xc(j))*sn+(cx(i)*bx(i)+by(i)*cy(i)-xc(j)*bx(i)-by(i)*yc(k))
                iter=iter+1
                if (iter>1000) goto 80
                if (dabs(res)<eps) then 
                    goto 80
                else 
                    s0=sn
                    goto 70
                end if 
80              continue
                xn=ax(i)*sn**2+bx(i)*sn+cx(i)
                yn=ay(i)*sn**2+by(i)*sn+cy(i)
                xd=2.0d0*ax(i)*sn+bx(i)           !d(xsn)/ds
                yd=2.0d0*ay(i)*sn+by(i)           !d(ysn)/ds
                !compute x and y normals at sn 
                nxx=-yd/dsqrt(xd**2+yd**2)       !nx
                nyy=xd/dsqrt(xd**2+yd**2)        !ny
                !backup i ndex, j index, xn and yn values of packet points inorder to find virtual point later
                ppacket(pp,1)=j
                ppacket(pp,2)=k
                ppacket(pp,3)=xn
                ppacket(pp,4)=yn
                ppacket(pp,5)=nxx
                ppacket(pp,6)=nyy
                ppacket(pp,7)=sn
                !compute vector lambda
                lx=(xn-xc(j))/dsqrt((xn-xc(j))**2+(yn-yc(k))**2)     !Lambda_x
                ly=(yn-yc(k))/dsqrt((xn-xc(j))**2+(yn-yc(k))**2)     !Lambda_y
                delta=nxx*lx+nyy*ly
                if (delta>0.0d0) then     
                    flagp(j,k)=1             !(j,k) point lies inside the object (solid phase)
                else
                    flagp(j,k)=0            ! (j,k) point lies outside the object (fluid phase)
                    vos(j,k)=0.0d0
                end if
            end do           
        end do !end of packet of points loop
    end do     !end of ib point searchflagib
    ! identify internal boundaries to 1 for u, v and p    
    call ASSIGNINTERNALBOUND(problem, imax, jmax, flagu)  
    call ASSIGNINTERNALBOUND(problem, imax, jmax, flagv)  
    call ASSIGNINTERNALBOUND(problem, imax, jmax, flagp) 
    
    do i=1, imax
        do j=1, jmax
            if (flagp(i,j)==1) pflag(i,j)=C_B
            if (flagu(i,j)==1) uflag(i,j)=C_B
            if (flagv(i,j)==1) vflag(i,j)=C_B
        end do
    end do
150 continue   
    !
    ! flags for boundary cells
    !
    ibound=0
    do i = 1, imax
        do j = 1, jmax
            if ( IAND(pflag(i,j),C_F) /= C_F) then
                ibound = ibound + 1
            end if
            pflag(i,j) = pflag(i,j) + ( IAND(pflag(i-1,j),C_F) * B_W  &
            +   IAND(pflag(i+1,j),C_F) * B_E  &
            +   IAND(pflag(i,j-1),C_F) * B_S  &
            +   IAND(pflag(i,j+1),C_F) * B_N ) / C_F
        end do
    end do   
    do i = 1, imax
        do j = 1, jmax
            uflag(i,j) = uflag(i,j) + ( IAND(uflag(i-1,j),C_F) * B_W  &
            +   IAND(uflag(i+1,j),C_F) * B_E  &
            +   IAND(uflag(i,j-1),C_F) * B_S  &
            +   IAND(uflag(i,j+1),C_F) * B_N ) / C_F
        end do
    end do   
    do i = 1, imax
        do j = 1, jmax
            vflag(i,j) = vflag(i,j) + ( IAND(vflag(i-1,j),C_F) * B_W  &
            +   IAND(vflag(i+1,j),C_F) * B_E  &
            +   IAND(vflag(i,j-1),C_F) * B_S  &
            +   IAND(vflag(i,j+1),C_F) * B_N ) / C_F
        end do
    end do 
    !set active  points according to the algorithm defined by Faltinsen
    do i=1, imax
        do j=1, jmax
            if ((IAND(pflag(i,j),C_F)/=0).AND.(IAND(pflag(i+1,j), C_F)/=0)) then
                ua(i,j)=1
                write(29,*) ux(i), uy(j)
            else
                ua(i,j)=0
            end if
            
            if ((IAND(pflag(i,j),C_F)/=0).AND.(IAND(pflag(i,j+1), C_F)/=0)) then
                va(i,j)=1
                write(30,*) vx(i), vy(j)
            else
                va(i,j)=0
            end if
            if (IAND(pflag(i,j),C_F)/=0) then
                pa(i,j)=1
            else
                pa(i,j)=0
            end if
        end do
    end do
    if (t==0.0) then
        write (17,*) 'VARIABLES= ','"X",','"Z",','"flag"' 
        write (17,*) 'ZONE T="mesh " i=',imax+1,',j=',jmax+1,',   F=POINT'
        do j=0, jmax
            do i=0, imax
                write(17,*) xg(i), yg(j), 1.0
            end do
        end do
        write (17,*) 'VARIABLES= ','"X",','"Z",','"flag"' 
        write (17,*) 'ZONE T="u-cell " i=',imax,',j=',jmax,',   F=POINT'
        do j=1, jmax
            do i=1, imax
                write(17,*) ux(i), uy(j), flagu(i,j)
            end do
        end do
        write (17,*) 'VARIABLES= ','"X",','"Z",','"flag"' 
        write (17,*) 'ZONE T="v-cell"  i=',imax,',j=',jmax,',   F=POINT'
        do j=1, jmax
            do i=1, imax
                write(17,*) vx(i), vy(j), flagv(i,j)
            end do
        end do
        write (17,*) 'VARIABLES= ','"X",','"Z",','"flag"' 
        write (17,*) 'ZONE T="pcell" i=',imax,',j=',jmax,',   F=POINT'
        do j=1, jmax
            do i=1, imax
                write(17,*) xc(i), yc(j), flagp(i,j)
            end do
        end do
        write (17,*) 'VARIABLES= ','"X",','"Y",','"flagib"' 
        write (17,*) 'ZONE T="pflag"  i=',imax,',j=',jmax,',   F=POINT'
        do j=1, jmax
            do i=1, imax
                write(17,*) xc(i), yc(j), pflag(i,j)
            end do
        end do
        write (17,*) 'VARIABLES= ','"X",','"Y",','"flagib"' 
        write (17,*) 'ZONE T="uflag"  i=',imax,',j=',jmax,',   F=POINT'
        do j=1, jmax
            do i=1, imax
                write(17,*) ux(i), uy(j), uflag(i,j)
            end do
        end do
        write (17,*) 'VARIABLES= ','"X",','"Y",','"flagib"' 
        write (17,*) 'ZONE T="vflag"  i=',imax,',j=',jmax,',   F=POINT'
        do j=1, jmax
            do i=1, imax
                write(17,*) vx(i), vy(j), vflag(i,j)
            end do
        end do
    end if
    !determine number of ghost cells for u,v and p 
    ppu=0
    ppv=0
    ppp=0
    do i=1, imax
        do j=1, jmax
            if ( (uflag(i,j) >= B_N) .AND. (uflag(i,j) <= B_SWE) ) ppu=ppu+1
            if ( (vflag(i,j) >= B_N) .AND. (vflag(i,j) <= B_SWE) ) ppv=ppv+1
            if ( (pflag(i,j) >= B_N) .AND. (pflag(i,j) <= B_SWE) ) ppp=ppp+1
        end do
    end do
   
    allocate (invagu(ppu,3,3))
    allocate (invagv(ppv,3,3))
    allocate (invagp(ppp,3,3))
    allocate (gFIu(ppu,2,2))            !ideces of inerpolation points for u
    allocate (gFIv(ppv,2,2))            !ideces of inerpolation points for v
    allocate (gFIp(ppp,2,2))            !ideces of inerpolation points for p
    allocate (indu(ppu,2))              !indeces of ghost cell for u 
    allocate (indv(ppv,2))              !indeces of ghost cell for v 
    allocate (indp(ppp,2))              !indeces of ghost cell for p 
    allocate (ghu(ppu,7))               !x and y coordinates of the ghost cell for u
    allocate (ghv(ppv,7))               !x and y coordinates of the ghost cell for v
    allocate (ghp(ppp,8))               !x and y coordinates of the ghost cell for p
    allocate (IMu(ppu))                 !information of whether the imaginary of ghost pont is used for u
    allocate (IMv(ppv))                 !information of whether the imaginary of ghost pont is used for v
    allocate (IMp(ppp))                 !information of whether the imaginary of ghost pont is used for p
    do i=1, ppu
        IMu(i)=0
    end do
    do i=1, ppv
        IMv(i)=0
    end do
    do i=1, ppp
        IMp(i)=0
    end do
    !store interpolation coefficients and their addresses to arrays depending on the orientation of the normal vector for u
    ppu=0
    do i=1, imax
        do j=1, jmax
            if ( (uflag(i,j) >= B_N) .AND. (uflag(i,j) <= B_SWE) ) then
                ppu=ppu+1
                indu(ppu,1)=i            !i index of ghost cell
                indu(ppu,2)=j            !j index of ghost cell
                ghu(ppu,1)=ux(i)         !x coordinate of the ghost cell
                ghu(ppu,2)=uy(j)         !y coordinate of the ghost cell
                write (18,*) ghu(ppu,1), ghu(ppu,2)
                found=.FALSE.
                k=0
                !find the ghost point in the packet array to restore the detailed information about the point
                do while(.NOT.(found))
                    k=k+1
                    if (i==int(upacket(k,1)).AND.j==int(upacket(k,2))) then
                        found=.TRUE.
                        xn=upacket(k,3)         !x coordinate of the intersection point of normal and curve
                        yn=upacket(k,4)         !y coordinate of the intersection point of normal and curve
                        ghu(ppu,3)=xn
                        ghu(ppu,4)=yn
                        ghu(ppu,5)=upacket(k,7) !arclength coordinate of the intersection point of ghost cell with the immersed boundary
                        if (ghu(ppu,5).lt.0.0d0) then 
                            if (problem=="circle") then 
                                ghu(ppu,5)=P+ghu(ppu,5)
                            else
                                ghu(ppu,5)=dmax1(0.0d0,ghu(ppu,5))
                            end if
                        end if
                        r(1,1)=1.0d0
                        r(2,1)=1.0d0
                        r(3,1)=1.0d0
                        r(1,2)=xn
                        r(1,3)=yn
                        nxx=upacket(k,5)
                        nyy=upacket(k,6)
                        nxxp=DINT(nxx*1.d6)/1.d6
                        nyyp=DINT(nyy*1.d6)/1.d6
                        !if the immersed boundary contacts with the boundary (arch dam problem), do not seek the active-cell below the current cell
                        if (j==1.and.nyyp.lt.0.0d0) then 
                            nyyp=0.0d0
                        end if
                        
                        
                        !find the nearest fluid cell in x direction 
                        if (nxxp.gt.0.0d0) then
                            ii=FINDFLUIDX(imax, ua, i, j, 1)
                        else if (nxxp.lt.0.0d0) then
                            ii=FINDFLUIDX(imax, ua, i, j, -1)  
                        end if
                        if (nxxp.eq.0.0d0) then 
                            r(2,2)=1.0d20
                            r(2,3)=1.0d20
                        else
                            r(2,2)=ux(ii)
                            r(2,3)=uy(j)
                        end if
                        gFIu(ppu,1,1)=ii                !i index of first interpolation point
                        gFIu(ppu,1,2)=j                 !j index of first interpolation point
                         !find the nearest fluid cell in y direction
                        if (nyyp.gt.0.0d0) then
                            jj=FINDFLUIDY(jmax, ua, i, j, 1)
                        elseif (nyyp.lt.0.0d0) then
                            jj=FINDFLUIDY(jmax, ua, i, j, -1)
                        end if
                        if (nyyp.eq.0.0d0) then 
                            r(3,2)=1.0d20
                            r(3,3)=1.0d20
                        else
                            r(3,2)=ux(i)
                            r(3,3)=uy(jj)
                        end if
                        gFIu(ppu,2,1)=i                 !i index of second interpolation point
                        gFIu(ppu,2,2)=jj                !j index of second interpolation point
                    end if
                end do
                !compute the inversion of matrix A for every ghost cell
                call MATRIXINVERSA(ppu, r, invagu)
                !one of the fluid nodes used in the extrapolation
                !check whether the first interpolation fluid cell is near to solid boundary
                !
                found=.false.
                i1=gFIu(ppu,1,1)    !i index of first interpolation cell
                j1=gFIu(ppu,1,2)    !j index of first interpolation cell
                !find this point in the packet array
                k=0
                dist1=1.0d10
                do while(.NOT.(found).and.(k.lt.9*N))
                    k=k+1
                    if (i1==int(upacket(k,1)).AND.j1==int(upacket(k,2))) then
                        found=.true.
                        xn=upacket(k,3)         !x coordinate of the intersection point of normal and curve
                        yn=upacket(k,4)         !y coordinate of the intersection point of normal and curve
                        xx=ux(i1)
                        yy=uy(j1)
                        dist1=dsqrt((xx-xn)**2+(yy-yn)**2)
                    end if 
                end do
                !
                !check whether the second interpolation fluid cell is near to solid boundary
                !
                found=.false.
                i2=gFIu(ppu,2,1)    !i index of second interpolation cell
                j2=gFIu(ppu,2,2)    !j index of second interpolation cell
                !find this point in the packet 
                k=0
                dist2=1.0d10
                do while(.NOT.(found).and.(k.lt.9*N)) !if the point is not neighbour to the immersed points
                    k=k+1
                    if (i2==int(upacket(k,1)).AND.j2==int(upacket(k,2))) then
                        found=.true.
                        xn=upacket(k,3)         !x coordinate of the intersection point of normal and curve
                        yn=upacket(k,4)         !y coordinate of the intersection point of normal and curve
                        xx=ux(i2)
                        yy=uy(j2)
                        dist2=dsqrt((xx-xn)**2+(yy-yn)**2)
                    end if 
                end do
                !if one of the interpolation fluid cells is near to solid boundary then take the imagine of the ghost point 
                if ((dist1.le.0.1*min(dx(i1),dy(j1))).or.(dist2.le.0.1*min(dx(i2),dy(j2)))) then
                    x1=ux(i)
                    y1=uy(j)
                    x2=r(1,2)
                    y2=r(1,3)
                    xim=x2+(x2-x1)              !x coordinate of image of ghost point
                    yim=y2+(y2-y1)              !y coordinate of image of ghost point
                    IMu(ppu)=1                  !imagine of ghost point was used
                    ghu(ppu,1)=xim              !x coordinate of image of ghost cell
                    ghu(ppu,2)=yim              !y coordinate of image of ghost cell
                    write (31,*) xim, yim
                end if                
            end if
        end do
    end do
    !store interpolation coefficients and their addresses to arrays depending on the orientation of the normal vector for v
    ppv=0
    do i=1, imax
        do j=1, jmax
            if ( (vflag(i,j) >= B_N) .AND. (vflag(i,j) <= B_SWE) ) then
                ppv=ppv+1
                indv(ppv,1)=i
                indv(ppv,2)=j
                ghv(ppv,1)=vx(i)         !x coordinate of the ghost cell
                ghv(ppv,2)=vy(j)         !y coordinate of the ghost cell
                write (19,*) ghv(ppv,1), ghv(ppv,2)
                found=.FALSE.
                k=0
                !find the ghost cell point in the packet array
                do while(.NOT.(found))
                    k=k+1
                    if (i==int(vpacket(k,1)).AND.j==int(vpacket(k,2))) then
                        found=.TRUE.
                        xn=vpacket(k,3)         !x coordinate of the intersection point of normal and curve
                        yn=vpacket(k,4)         !y coordinate of the intersection point of normal and curve
                        ghv(ppv,3)=xn
                        ghv(ppv,4)=yn
                        ghv(ppv,5)=vpacket(k,7) !arclength coordinate of the intersection point of ghost cell with the immersed boundary
                        if (ghv(ppv,5).lt.0.0d0) then
                            if (problem=="circle") then 
                                ghu(ppu,5)=P+ghu(ppu,5)
                            else
                                ghu(ppu,5)=dmax1(0.0d0,ghu(ppu,5))
                            end if
                        end if
                        r(1,1)=1.0d0
                        r(2,1)=1.0d0
                        r(3,1)=1.0d0
                        r(1,2)=xn
                        r(1,3)=yn
                        nxx=vpacket(k,5)
                        nyy=vpacket(k,6)
                        nxxp=DINT(nxx*1.d6)/1.d6
                        nyyp=DINT(nyy*1.d6)/1.d6
                        !if the immersed boundary contacts with the boundary (arch dam problem), do not seek the active-cell below the current cell
                        if (j==1.and.nyyp.lt.0.0d0) then 
                            nyyp=0.0d0
                        end if
                        !find the active cell in y direction
                        if (nxxp.gt.0.0d0) then
                            ii=FINDFLUIDX(imax, va, i, j, 1)
                        elseif (nxxp.lt.0.0d0) then 
                            ii=FINDFLUIDX(imax, va, i, j, -1)  
                        end if
                        if (nxxp.eq.0.0d0) then 
                            r(2,2)=1.0d20
                            r(2,3)=1.0d20
                        else
                            r(2,2)=vx(ii)
                            r(2,3)=vy(j)
                        end if
                        gFIv(ppv,1,1)=ii                !i index of first interpolation point
                        gFIv(ppv,1,2)=j                 !j index of first interpolation point
                        if (nyyp.gt.0.0d0) then
                            jj=FINDFLUIDY(jmax, va, i, j, 1)
                        elseif (nyyp.lt.0.0d0) then
                            jj=FINDFLUIDY(jmax, va, i, j, -1)
                        end if
                        if (nyyp.eq.0.0d0) then 
                            r(3,2)=1.0d20
                            r(3,3)=1.0d20
                        else
                            r(3,2)=vx(i)
                            r(3,3)=vy(jj)
                        end if
                        gFIv(ppv,2,1)=i                 !i index of second interpolation point
                        gFIv(ppv,2,2)=jj                !j index of second interpolation point
                    end if
                end do
                !compute the inversion of matrix A at every forcing point
                call MATRIXINVERSA(ppv, r, invagv)
                !
                !check whether the first interpolation fluid cell is near to solid boundary
                !
                found=.false.
                i1=gFIv(ppv,1,1)    !i index of first interpolation cell
                j1=gFIv(ppv,1,2)    !j index of first interpolation cell
                !find this point in the packet array
                k=0
                dist1=1.0d10
                do while(.NOT.(found).and.(k.lt.9*N))
                    k=k+1
                    if (i1==int(vpacket(k,1)).AND.j1==int(vpacket(k,2))) then
                        found=.true.
                        xn=vpacket(k,3)         !x coordinate of the intersection point of normal and curve
                        yn=vpacket(k,4)         !y coordinate of the intersection point of normal and curve
                        xx=vx(i1)
                        yy=vy(j1)
                        dist1=dsqrt((xx-xn)**2+(yy-yn)**2)
                    end if 
                end do
                !
                !check whether the second interpolation fluid cell is near to solid boundary
                !
                found=.false.
                i2=gFIv(ppv,2,1)    !i index of second interpolation cell
                j2=gFIv(ppv,2,2)    !j index of second interpolation cell
                !find this point in the packet 
                k=0
                dist2=1.0d10
                do while(.NOT.(found).and.(k.lt.9*N))
                    k=k+1
                    if (i2==int(vpacket(k,1)).AND.j2==int(vpacket(k,2))) then
                        found=.true.
                        xn=vpacket(k,3)         !x coordinate of the intersection point of normal and curve
                        yn=vpacket(k,4)         !y coordinate of the intersection point of normal and curve
                        xx=vx(i2)
                        yy=vy(j2)
                        dist2=dsqrt((xx-xn)**2+(yy-yn)**2)
                    end if 
                end do
                !if one of the interpolation fluid cells is near to solid boundary then take the imagine of ghost point 
                if ((dist1.le.0.1*min(dx(i1),dy(j1))).or.(dist2.le.0.1*min(dx(i2),dy(j2)))) then
                    x1=vx(i)
                    y1=vy(j)
                    x2=r(1,2)
                    y2=r(1,3)
                    xim=x2+(x2-x1)              !x coordinate of image of ghost point
                    yim=y2+(y2-y1)              !y coordinate of image of ghost point
                    IMv(ppv)=1                  !imagine of ghost point is used
                    ghv(ppv,1)=xim              !x coordinate of imagine of ghost cell
                    ghv(ppv,2)=yim              !y coordinate of imagine of ghost cell
                    write (32,*) xim, yim
                end if
            end if
        end do
    end do
    !store interpolation coefficients and their addresses to arrays depending on the orientation of the normal vector for p
    ppp=0
    do i=1, imax
        do j=1, jmax
            if ( (pflag(i,j) >= B_N) .AND. (pflag(i,j) <= B_SWE) ) then
                ppp=ppp+1
                indp(ppp,1)=i
                indp(ppp,2)=j
                ghp(ppp,1)=xc(i)         !x coordinate of the ghost cell
                ghp(ppp,2)=yc(j)         !y coordinate of the ghost cell
                write (20,*) ghp(ppp,1), ghp(ppp,2)
                found=.FALSE.
                k=0
                !find the ghost cell point in the packet array
                do while(.NOT.(found))
                    k=k+1
                    if (i==int(ppacket(k,1)).AND.j==int(ppacket(k,2))) then
                        found=.TRUE.
                        xnn=ppacket(k,3)         !x coordinate of the intersection point of normal and curve
                        ynn=ppacket(k,4)         !y coordinate of the intersection point of normal and curve
                        ghp(ppp,3)=xnn
                        ghp(ppp,4)=ynn
                        ghp(ppp,5)=ppacket(k,7) !arclength coordinate of the intersection point of ghost cell with the immersed boundary
                        if (ghp(ppp,5).lt.0.0d0) then 
                            if (problem=="circle") then 
                                ghu(ppu,5)=P+ghu(ppu,5)
                            else
                                ghu(ppu,5)=dmax1(0.0d0,ghu(ppu,5))
                            end if
                        end if
                        r(1,1)=0.0d0
                        r(2,1)=1.0d0
                        r(3,1)=1.0d0
                        nxx=ppacket(k,5)
                        nyy=ppacket(k,6)  
                        ghp(ppp,7)=nxx
                        ghp(ppp,8)=nyy
                        !Tseng and Ferziger formulation (it gives the same with Balaras formulation)
                        !teta=nyy/nxx
                        !r(1,2)=-sind(teta)
                        !r(1,3)=cosd(teta)
                        r(1,2)=nxx
                        r(1,3)=nyy
                        nxxp=DINT(nxx*1.d6)/1.d6
                        nyyp=DINT(nyy*1.d6)/1.d6
                        !if the immersed boundary contacts with the boundary (arch dam problem), do not seek the active-cell below the current cell
                        if (j==1.and.nyyp.lt.0.0d0) then 
                            nyyp=0.0d0
                        end if
                        if (nxxp.gt.0.0d0) then
                            ii=FINDFLUIDX(imax, pa, i, j, 1)
                        elseif (nxxp.lt.0.0d0) then 
                            ii=FINDFLUIDX(imax, pa, i, j, -1)  
                        end if
                        if (nxxp.eq.0.0d0) then 
                            r(2,2)=1.0d20
                            r(2,3)=1.0d20
                        else
                            r(2,2)=xc(ii)
                            r(2,3)=yc(j)
                        end if
                        gFIp(ppp,1,1)=ii                !i index of first interpolation point
                        gFIp(ppp,1,2)=j                 !j index of first interpolation point
                        if (nyyp.gt.0.0d0) then
                            jj=FINDFLUIDY(jmax, pa, i, j, 1)
                        elseif (nyyp.lt.0.0d0) then
                            jj=FINDFLUIDY(jmax, pa, i, j, -1)
                        end if  
                        if (nyyp.eq.0.0d0) then 
                            r(3,2)=1.0d20
                            r(3,3)=1.0d20
                        else
                            r(3,2)=xc(i)
                            r(3,3)=yc(jj)
                        end if
                        gFIp(ppp,2,1)=i                 !i index of second interpolation point
                        gFIp(ppp,2,2)=jj                !j index of second interpolation point
                    end if
                end do
                !compute the inversion of matrix B at every forcing point
                call MATRIXINVERSB(ppp, r, invagp)
                !
                !check whether the first interpolation fluid cell is near to solid boundary
                !
                found=.false.
                i1=gFIp(ppp,1,1)    !i index of first interpolation cell
                j1=gFIp(ppp,1,2)    !j index of first interpolation cell
                !find this point in the packet array
                k=0
                dist1=1.0d10
                do while(.NOT.(found).and.(k.lt.9*N))
                    k=k+1
                    if (i1==int(ppacket(k,1)).AND.j1==int(ppacket(k,2))) then
                        found=.true.
                        xn=ppacket(k,3)         !x coordinate of the intersection point of normal and curve
                        yn=ppacket(k,4)         !y coordinate of the intersection point of normal and curve
                        xx=xc(i1)
                        yy=yc(j1)
                        dist1=dsqrt((xx-xn)**2+(yy-yn)**2)
                    end if 
                end do
                !
                !check whether the second interpolation fluid cell is near to solid boundary
                !
                found=.false.
                i2=gFIp(ppp,2,1)    !i index of second interpolation cell
                j2=gFIp(ppp,2,2)    !j index of second interpolation cell
                !find this point in the packet 
                k=0
                dist2=1.0d10
                do while(.NOT.(found).and.(k.lt.9*N))
                    k=k+1
                    if (i2==int(ppacket(k,1)).AND.j2==int(ppacket(k,2))) then
                        found=.true.
                        xn=ppacket(k,3)         !x coordinate of the intersection point of normal and curve
                        yn=ppacket(k,4)         !y coordinate of the intersection point of normal and curve
                        xx=xc(i2)
                        yy=yc(j2)
                        dist2=dsqrt((xx-xn)**2+(yy-yn)**2)
                    end if 
                end do
                !if one of the interpolation fluid cells is near to solid boundary then take the imagine of ghost point 
                if ((dist1.le.0.1*min(dx(i1),dy(j1))).or.(dist2.le.0.1*min(dx(i2),dy(j2)))) then
                    x1=xc(i)
                    y1=yc(j)
                    x2=xnn
                    y2=ynn
                    xim=x2+(x2-x1)              !x coordinate of image of ghost point
                    yim=y2+(y2-y1)              !y coordinate of image of ghost point
                    IMp(ppp)=1                  !imagine of ghost point is used
                    ghp(ppp,1)=xim              !x coordinate of imagine of ghost cell
                    ghp(ppp,2)=yim              !y coordinate of imagine of ghost cell
                    write (33,*) xim, yim
                end if
            end if
        end do
    end do
    !write intersection points of ghost points and solid boundary
    write (20,*) 'ZONE T="p xn yn p"'
    do i=1, ppp
        write (20,*) ghp(i,3), ghp(i,4)
    end do
    write (18,*) 'ZONE T="u xn yn u"'
    do i=1, ppu
        write (18,*) ghu(i,3), ghu(i,4)
    end do
    write (19,*) 'ZONE T="v xn yn v"'
    do i=1, ppv
        write (19,*) ghv(i,3), ghv(i,4)
    end do
    !determine the number of inactive cells which lie in the fluid for u and v 
    ivp=0
    jvp=0
    do i=2, imax-1
        do j=2, jmax-1
            if ((IAND(uflag(i,j),C_F)/=0).AND.(ua(i,j)==0)) then 
                ivp=ivp+1
            end if
            if ((IAND(vflag(i,j),C_F)/=0).AND.(va(i,j)==0)) then 
                jvp=jvp+1
            end if
        end do
    end do
    allocate (invaguin(ivp,3,3))
    allocate (invagvin(jvp,3,3))
    allocate (gFIuin(ivp,2,2))            !ideces of inerpolation points for u
    allocate (gFIvin(jvp,2,2))            !ideces of inerpolation points for v
    allocate (induin(ivp,2))              !indeces of inactive cell for u 
    allocate (indvin(jvp,2))              !indeces of inactive cell for v 
    allocate (ghuin(ivp,7))               !x and y coordinates of the inactive cell for u
    allocate (ghvin(jvp,7))               !x and y coordinates of the inactive cell for v
    
    !store interpolation coefficients and their addresses to arrays depending on the orientation of the normal vector for u (inactive cell)
    m=0
    do i=2, imax-1
        do j=2, jmax-1
            if ((IAND(uflag(i,j),C_F)/=0).AND.(ua(i,j)==0)) then
                m=m+1
                induin(m,1)=i            !i index of inactive cell
                induin(m,2)=j            !j index of inactive cell
                ghuin(m,1)=ux(i)         !x coordinate of the inactive cell
                ghuin(m,2)=uy(j)         !y coordinate of the inactive cell
                write (13,*) ux(i), uy(j)
                found=.FALSE.
                k=0
                !find the inactive point in the packet array to restore the detailed information about the point
                do while(.NOT.(found))
                    k=k+1
                    if (i==int(upacket(k,1)).AND.j==int(upacket(k,2))) then
                        found=.TRUE.
                        xn=upacket(k,3)         !x coordinate of the intersection point of normal and curve
                        yn=upacket(k,4)         !y coordinate of the intersection point of normal and curve
                        ghuin(m,3)=xn
                        ghuin(m,4)=yn
                        r(1,1)=1.0d0
                        r(2,1)=1.0d0
                        r(3,1)=1.0d0
                        r(1,2)=xn
                        r(1,3)=yn
                        nxx=upacket(k,5)
                        nyy=upacket(k,6)
                        nxxp=DINT(nxx*1.0d6)/1.0d6
                        nyyp=DINT(nyy*1.0d6)/1.0d6
                        !find the nearest active cell in x direction 
                        if (nxxp.gt.0.0d0) then
                            ii=FINDFLUIDX(imax, ua, i, j, 1)
                        elseif (nxxp.lt.0.0d0) then
                            ii=FINDFLUIDX(imax, ua, i, j, -1)  
                        end if
                        if (nxxp.eq.0.0d0) then 
                            r(2,2)=1.0d20
                            r(2,3)=1.0d20
                        else
                            r(2,2)=ux(ii)
                            r(2,3)=uy(j)
                        end if
                        gFIuin(m,1,1)=ii                !i index of first interpolation point
                        gFIuin(m,1,2)=j                 !j index of first interpolation point
                        !find the nearest active cell in y direction
                        if (nyyp.gt.0.0d0) then
                            jj=FINDFLUIDY(jmax, ua, i, j, 1)
                        elseif (nyyp.lt.0.0d0) then
                            jj=FINDFLUIDY(jmax, ua, i, j, -1)
                        end if  
                        if (nyyp.eq.0.0d0) then 
                            r(3,2)=1.0d20
                            r(3,3)=1.0d20
                        else
                            r(3,2)=ux(i)
                            r(3,3)=uy(jj)
                        end if
                        gFIuin(m,2,1)=i                 !i index of second interpolation point
                        gFIuin(m,2,2)=jj                !j index of second interpolation point
                    end if
                end do
                !compute the inversion of matrix A for every inactive cell
                call MATRIXINVERSA(m, r, invaguin)                
            end if
        end do
    end do
    !store interpolation coefficients and their addresses to arrays depending on the orientation of the normal vector for v
    m=0
    do i=2, imax-1
        do j=2, jmax-1
            if ((IAND(vflag(i,j),C_F)/=0).AND.(va(i,j)==0)) then 
                m=m+1
                indvin(m,1)=i
                indvin(m,2)=j
                ghvin(m,1)=vx(i)         !x coordinate of the inactive cell
                ghvin(m,2)=vy(j)         !y coordinate of the inactive cell
                write (14,*) vx(i), vy(j)
                found=.FALSE.
                k=0
                !find the inactive cell point in the packet array
                do while(.NOT.(found))
                    k=k+1
                    if (i==int(vpacket(k,1)).AND.j==int(vpacket(k,2))) then
                        found=.TRUE.
                        xn=vpacket(k,3)         !x coordinate of the intersection point of normal and curve
                        yn=vpacket(k,4)         !y coordinate of the intersection point of normal and curve
                        ghvin(m,3)=xn
                        ghvin(m,4)=yn
                        r(1,1)=1.0d0
                        r(2,1)=1.0d0
                        r(3,1)=1.0d0
                        r(1,2)=xn
                        r(1,3)=yn
                        nxx=vpacket(k,5)
                        nyy=vpacket(k,6)
                        nxxp=DINT(nxx*1.0d6)/1.0d6
                        nyyp=DINT(nyy*1.0d6)/1.0d6
                        !find the active cell in y direction
                        if (nxxp.gt.0.0d0) then
                            ii=FINDFLUIDX(imax, va, i, j, 1)
                        elseif (nxxp.lt.0.0d0) then  
                            ii=FINDFLUIDX(imax, va, i, j, -1)  
                        end if
                        if (nxxp.eq.0.0d0) then 
                            r(2,2)=1.0d20
                            r(2,3)=1.0d20
                        else
                            r(2,2)=vx(ii)
                            r(2,3)=vy(j)
                        end if
                        gFIvin(m,1,1)=ii                !i index of first interpolation point
                        gFIvin(m,1,2)=j                 !j index of first interpolation point
                        if (nyyp.gt.0.0d0) then
                            jj=FINDFLUIDY(jmax, va, i, j, 1)
                        elseif (nyyp.lt.0.0d0) then
                            jj=FINDFLUIDY(jmax, va, i, j, -1)
                        end if 
                        if (nyyp.eq.0.0d0) then 
                            r(3,2)=1.0d20
                            r(3,3)=1.0d20
                        else
                            r(3,2)=vx(i)
                            r(3,3)=vy(jj)
                        end if
                        gFIvin(m,2,1)=i                 !i index of second interpolation point
                        gFIvin(m,2,2)=jj                !j index of second interpolation point
                    end if
                end do
                !compute the inversion of matrix A at every forcing point
                call MATRIXINVERSA(m, r, invagvin)
            end if
        end do
    end do
    !
    !write arclength coordinates of intersection points for u, v and p
    !
    !write arclength coordinate of the intersection point for pressure
    do i=1, ppp
        xx=ghp(i,3)
        yy=ghp(i,4)
        ss=ghp(i,5)
        write (37,*) xx, yy, ss
    end do
    !
    !write arclength coordinate of the intersection point for u
    !
    do i=1, ppu
        xx=ghu(i,3)
        yy=ghu(i,4)
        ss=ghu(i,5)
        write (38,*) xx, yy, ss
    end do
    !
    !write arclength coordinate of the intersection point for v
    !
    do i=1, ppv
        xx=ghv(i,3)
        yy=ghv(i,4)
        ss=ghv(i,5)
        write (39,*) xx, yy, ss
    end do
      
    return
    end subroutine IBM

    subroutine FORCING(imax, jmax, ppu, ppv, U, V, F, G, dx, dy, delt, t, indu, indv, gFIu, gFIv, invagu, invagv, ghu, &
                        ghv, IMu, IMv, au, bu, cu, av, bv, cv, uxb, uyb)
    use nrtype
    use defs
    implicit none
    integer, intent(IN) :: imax, jmax, ppu, ppv     
    real(RP), dimension(0:,0:), intent(INOUT) :: U, V
    real(RP), dimension(0:,0:), intent(INOUT) :: F, G
    real(RP), dimension(0:), intent(IN) :: dx, dy
    real(RP), intent(IN):: delt, t, uxb, uyb
    integer, dimension(:,:), intent(IN) :: indu, indv
    integer, dimension(:,:,:), intent(IN) :: gFIu, gFIv
    real(DP), dimension(:,:,:), intent(IN) :: invagu, invagv
    real(DP), dimension(:,:), intent(IN) :: ghu, ghv
    integer, dimension(:), intent(IN) :: IMu, IMv
    real(DP), dimension(:), intent(OUT) :: au, bu, cu, av, bv, cv
    !local variables
    integer i, j, k, ii, jj, i1, j1, i2, j2
    real(DP) a0, a1, a2, fd, fl, cd, cl, fiu, fiv, um
    real(RP), dimension(imax,jmax) :: fx, fy
    real(DP), dimension(3):: fi
    !open(34,FILE='dragf.plt',STATUS='UNKNOWN')
    
    
    um=1.0d0
    do i=1, imax
        do j=1, jmax
            fx(i,j)=0.0
            fy(i,j)=0.0
        end do
    end do
    
    do i=1, ppu
        ii=indu(i,1)
        jj=indu(i,2)
        i1=gFIu(i,1,1)          !i index of first interpolation point
        j1=gFIu(i,1,2)          !j index of first interpolation point
        i2=gFIu(i,2,1)          !i index of second interpolation point
        j2=gFIu(i,2,2)          !j index of second interpolation point
        fi(1)=uxb               !horizontal solid velocity 
        fi(2)=U(i1,j1)          !velocity at first interpolation point
        fi(3)=U(i2,j2)          !velocity at second interpolation point
        !compute linear interpolation constants a0, a1 and a2
        a0=fi(1)+fi(2)*invagu(i,1,2)+fi(3)*invagu(i,1,3)
        a1=fi(1)+fi(2)*invagu(i,2,2)+fi(3)*invagu(i,2,3)
        a2=fi(1)+fi(2)*invagu(i,3,2)+fi(3)*invagu(i,3,3)
        !compute u variable at ghost point using linear interpolation
        fiu=a0+a1*ghu(i,1)+a2*ghu(i,2)
        !if imagine of ghost point is used then re-compute the ghost point value
        if (IMu(i)==1) fiu=2.0d0*fi(1)-fiu
        fx(ii,jj)=(fiu-F(ii,jj))/delt   
        au(i)=a0
        bu(i)=a1
        cu(i)=a2
    end do
    
    do i=1, ppv
        ii=indv(i,1)
        jj=indv(i,2)
        i1=gFIv(i,1,1)          !i index of first interpolation point
        j1=gFIv(i,1,2)          !j index of first interpolation point
        i2=gFIv(i,2,1)          !i index of second interpolation point
        j2=gFIv(i,2,2)          !j index of second interpolation point
        fi(1)=uyb               !vertical solid velocity 
        fi(2)=V(i1,j1)          !velocity at first interpolation point
        fi(3)=V(i2,j2)          !velocity at second interpolation point
        !compute linear interpolation constants a0, a1 and a2
        a0=fi(1)+fi(2)*invagv(i,1,2)+fi(3)*invagv(i,1,3)
        a1=fi(1)+fi(2)*invagv(i,2,2)+fi(3)*invagv(i,2,3)
        a2=fi(1)+fi(2)*invagv(i,3,2)+fi(3)*invagv(i,3,3)
        !compute u variable at ghost point using linear interpolation
        fiv=a0+a1*ghv(i,1)+a2*ghv(i,2)
        !if imagine of ghost point is used then re-compute the ghost point value
        if (IMv(i)==1) fiv=2.0d0*fi(1)-fiv
        fy(ii,jj)=(fiv-G(ii,jj))/delt
        av(i)=a0
        bv(i)=a1
        cv(i)=a2
    end do
    !step 3 
    !update the tentative velocity field (u*) satisfying boundary conditions
    do i=2, imax-1
        do j=2, jmax-1
            F(i,j)=F(i,j)+delt*fx(i,j)
            G(i,j)=G(i,j)+delt*fy(i,j)
        end do
    end do
    !fd=0.0d0
    !fl=0.0d0
    !do i=2, imax-1
    !    do j=2, jmax-1                        
    !        fd=fd+fx(i,j)*dx(i)*dy(j)
    !        fl=fl+fy(i,j)*dx(i)*dy(j)
    !    end do
    !end do
    !fd=(-1.0d0)*fd
    !fl=(-1.0d0)*fl
    !cd=fd/(0.5d0)/um**2/1.0d0
    !cl=fl/(0.5d0)/um**2/1.0d0
    !!write (34,*) t, cd, cl
   
    return    
    end subroutine FORCING


    subroutine ARCGEN(problem, Mx, My, dx, dy, x, y, sc, N, t, h0, xlen, zlen, P, xtoe, dam_height, sw)
    use nrtype
    use defs
    implicit none
    character (LEN=30), INTENT(IN) :: problem
    integer, intent(IN) :: Mx, My
    real(RP), dimension(0:), intent(IN) :: dx, dy 
    real(DP), dimension(:), allocatable, intent(OUT) :: x, y, sc
    integer, intent(OUT) :: N
    real(RP), intent(IN) :: t, h0, xlen, zlen, xtoe
    real(DP), intent(OUT) :: P
    real(RP), intent(OUT) :: dam_height, sw
    
    character (LEN=30) :: sol_problem
    integer i, j, N1, N2, N3, mm
    real(DP):: dxmin, dymin, r, a, b, ds, s, pi, a0, fs, fsd, an, res, s1, s2, s3, s4, s5, &
                d, hs, k, L, m, nn, alpha, beta, ds1, ds2, ds3, xx, zz, eps, dmin, ss, xtoei, &
                hr, cr, scale, ht, disp, period
    real(DP), dimension(:), allocatable :: xt, yt, st
    logical inclined, partiallyinclined, arch, try
    OPEN(10,FILE='circle.plt',STATUS='UNKNOWN')
    OPEN(12,FILE='arch.plt',STATUS='UNKNOWN')
    open(11,FILE='grid.plt',STATUS='UNKNOWN')

    pi=datan(1.d0)*4.d0
    sol_problem="shelf"
    dxmin=1.d10
    dymin=1.d10

    do i=1, Mx
        if (dx(i)<dxmin) dxmin=dx(i)
    end do
    do j=1, My
        if (dy(j)<dymin) dymin=dy(j)
    end do
    dxmin=dxmin
    dymin=dymin
    
    if (problem=="circle") then
        try=.false.
        if (try) then
            s1=1.0d0
            s2=0.5d0
            s3=1.0d0*dcos(pi/6.0d0)
            P=s1+s2+s3
            N=max(int(P/dxmin),int(P/dymin))    !number of markers
            allocate(x(N),y(N),sc(N))
            ds=P/N        !distance between two marker points
            x(1)=9.0d0    !x coordinate of the first marker      
            y(1)=14.0d0   !y coordinate of the first marker
            s=0.0d0       !arclength coordinate of the first marker
            sc(1)=s       !store the arclength coordinate of a point to an array
            if (t==0.0) write (10,*) x(1), y(1)
            do i=2, N
                s=s+ds
                sc(i)=s
                if (s<s1) then
                    x(i)=x(1)+s*dcos(pi/6.0d0)                    !x coordinate of the marker
                    y(i)=y(1)+s*dsin(pi/6.0d0)                    !y coordinate of the marker
                elseif (s<s1+s2) then
                    x(i)=x(1)+1.0d0*dcos(pi/6.0d0)
                    y(i)=y(i-1)-ds
                else
                    x(i)=x(i-1)-ds
                    y(i)=y(1)
                end if
                if (t==0.0) write (10,*) x(i), y(i)
            end do
        else
            !
            !parameters of the circle
            !
            period=1.0d0                            !period of the oscillating cylinder
            if (movbound) then 
                disp=-0.5d0/pi*dsin(2.0d0*pi*t/period)  !displacement of the cylinder
            else
                disp=0.d0
            end if
            r=0.5d0                                 !radius of circle
            !r=0.2d0                                !radius of circle
            P=2.d0*pi*r                             !perimeter of circle
            N=min(int(P/dxmin),int(P/dymin))        !number of markers
            allocate(x(N),y(N),sc(N))
            a=10.0d0                                !x coordinate of center
            b=15.0d0                                !y coordinate of center
            if (movbound) then 
                a=xlen/2.0d0
                b=zlen/2.0d0
            end if
            !a=0.5d0                                !x coordinate of center
            !b=0.5d0                                !y coordinate of center
            ds = P/N                                !distance between two marker points
            x(1)=a-r                                !x coordinate of the first marker
            y(1)=b                                  !y coordinate of the first marker
            s=0.0d0                                 !arclength coordinate of the first marker
            sc(1)=s                                 !store the arclength coordinate of a point to an array
            write (10,*) x(1), y(1), sc(1)
            do i=2, N
                s=s+ds
                sc(i)=s                             ! arclength coordinate
                x(i)=a+r*dcos(pi-s/r)               ! x coordinate of the marker
                y(i)=b+r*dsin(pi-s/r)               ! y coordinate of the marker
                write (10,*) x(i), y(i), sc(i)
            end do
        end if
        
    elseif (problem=="arch") then
        !Mx=100
        !My=100
        !xlen=10.d0    !length of domain
        !ylen=10.d0    !height of domain
        s1=4.64678d0   !length of part 1 
        s2=5.d0         !length of part 2
        s3=0.5d0
        s4=7.d0
        s5=4.5d0
        P=s1+s2+s3+s4+s5    !total length of parts
        N=max(int(P/dxmin),int(P/dymin))    !number of markers
        allocate(x(N),y(N),sc(N))
        ds=P/N        !distance between two marker points
        x(1)=1.d0     !x coordinate of the first marker      
        y(1)=5.5d0    !y coordinate of the first marker
        s=0.d0        !arclength coordinate of the first marker
        sc(1)=s       !store the arclength coordinate of a point to an array
        if (t==0.0) write (12,*) x(1), y(1)
        do i=2, N
            s=s+ds
            sc(i)=s
            if (s<s1) then
                a0=sc(i)
10              fs=1.0/2.d0*dsqrt(5.d0-8.d0*a0+4.d0*a0**2)*(a0-1.d0)+1.d0/4.d0*dasinh(2.d0*(a0-1.d0))-s
                fsd=dsqrt(5.0-8.0*a0+4.0*a0**2)
                an=a0-fs/fsd
                res=1.d0/2.d0*dsqrt(5.d0-8.d0*an+4.d0*an**2)*(an-1.d0)+1.d0/4.d0*dasinh(2.d0*(an-1.d0))-s
                if (dabs(res)<1.d-10) then 
                    goto 20
                else 
                    a0=an
                    goto 10
                end if 
20              continue
                x(i)=an                    !x coordinate of the marker
                y(i)=5.5d0-(x(i)-1.d0)**2   !y coordinate of the marker
            elseif (s<s1+s2) then
                x(i)=3.d0+(s-s1)
                y(i)=1.5d0
            elseif (s<s1+s2+s3) then 
                x(i)=8.d0
                y(i)=1.5d0-(s-(s1+s2))
            elseif(s<s1+s2+s3+s4) then 
                x(i)=8.d0-(s-(s1+s2+s3))
                y(i)=1.d0
            else
                x(i)=1.d0
                y(i)=1.d0+s-(s1+s2+s3+s4)
            end if
            if (t==0.0) write (12,*) x(i), y(i)
        end do
    elseif (problem=="dcav") then
        s1=0.01d0
        s2=s1
        s3=s1
        s4=s1
        P=s1+s2+s3+s4
        N=max(int(P/dxmin),int(P/dymin))    !number of markers
        allocate(x(N),y(N),sc(N))
        ds=P/N                              !distance between two marker points
        x(1)=0.005d0                        !x coordinate of the first marker      
        y(1)=0.015d0                        !y coordinate of the first marker
        s=0.d0                              !arclength coordinate of the first marker
        sc(1)=s                             !store the arclength coordinate of a point to an array
        if (t==0.0) write (10,*) x(1), y(1)
        do i=2, N
            s=s+ds
            sc(i)=s
            if (s<=s1) then
                x(i)=x(i-1)
                y(i)=y(i-1)-ds
            elseif (s<=s1+s2) then
                x(i)=x(i-1)+ds
                y(i)=y(i-1)
            elseif (s<=s1+s2+s3) then
                x(i)=x(i-1)
                y(i)=y(i-1)+ds
            else
                x(i)=x(i-1)-ds
                y(i)=y(i-1)
            end if
            write (10,*) x(i), y(i)
        end do
    elseif (problem=="solitary") then
        if (sol_problem=="shelf") then
            m=1.0/20.0
            alpha=atan(m)
            d = 0.0381
            s1 = d/dsin(alpha)
            s2 = xlen-2.0-d/m
            P = s1+s2
            N=max(int(P/dxmin),int(P/dymin))    !number of markers
            allocate(x(N),y(N),sc(N))
            ds=P/N                              !distance between two marker points
            x(1)=2.0d0                       !x coordinate of the first marker      
            y(1)=0.d0
            s=0.d0                              !arclength coordinate of the first marker
            sc(1)=s                             !store the arclength coordinate of a point to an array
            if (t==0.0) write (10,*) x(1), y(1)
            do i=2, N
                s=s+ds
                sc(i)=s
                if (s.lt.s1) then
                    x(i)=x(i-1)+ds*dcos(alpha)
                    y(i)=y(i-1)+ds*dsin(alpha)
                else
                    x(i)=x(i-1)+ds
                    y(i)=d
                end if
                write (10,*) x(i), y(i)
            end do

        else
            m=1.0/19.85
            alpha=atan(m)
            d = 1.0			!initial water depth
            hs = 0.3        !wave height A 
            k = sqrt(3.0/4.0*hs/d**3)
            L = 2.0/k*acosh(sqrt(1.0/0.05))
            s1=1.5d0/dsin(alpha)    ! length of the plane beach
            s2=1.5d0
            s3=1.5d0/m
            P=s1+s2+s3
            N=max(int(P/dxmin),int(P/dymin))    !number of markers
            allocate(x(N),y(N),sc(N))
            ds=P/N                              !distance between two marker points
            x(1)=2.0d0*L                        !x coordinate of the first marker      
            !y(1)=3.0d0*dy(1)                    !y coordinate of the first marker
            y(1)=0.d0
            s=0.d0                              !arclength coordinate of the first marker
            sc(1)=s                             !store the arclength coordinate of a point to an array
            if (t==0.0) write (10,*) x(1)-38.21, y(1)
            do i=2, N
                s=s+ds
                sc(i)=s
                if (s.lt.s1) then
                    x(i)=x(i-1)+ds*dcos(alpha)
                    y(i)=y(i-1)+ds*dsin(alpha)
                elseif (s<s1+s2) then 
                    x(i)=x(1)+s3
                    y(i)=y(i-1)-ds
                else
                    x(i)=x(i-1)-ds
                    y(i)=y(1)
                end if
                write (10,*) x(i)-38.21, y(i)
            end do
        end if
        
   elseif (problem=="reservoir") then
       partiallyinclined=.true.
       inclined=.false.
       arch=.false.
       !if reservoir is partially inclined 
       if (partiallyinclined) then 
           xtoei=19.0d0
           cr=0.5d0
           hr=cr*h0
           alpha=pi/4.0d0
           sw=hr/dcos(alpha)+(h0-hr)                                    !wetted perimeter
           dam_height=1.1d0*dble(h0)
           s1=dble(hr)/dcos(alpha)                                      ! length of the plane beach
           s2=dam_height-hr                                             ! length of the vertical part
           dmin=dmin1(dxmin,dymin)
           r=5.0d0*dmin
           a=dble(xtoei)+hr*dtan(alpha)+r
           b=dam_height
           s3=2.0d0*pi*r/4.0d0                                          !length of fam crest
           s4=dble(xlen)-dble(xtoei)-hr*dtan(alpha)-r                   ! length of fam crest
           P=s1+s2+s3+s4
           N=max(int(P/(dxmin)),int(P/(dymin)))                 !number of markers
           allocate(x(N),y(N),sc(N))
           ds=P/N                                                       !distance between two marker points
           x(1)=xtoei                                                   !x coordinate of the first marker      
           y(1)=0.0d0                                                   !y coordinate of the first marker
           s=0.d0                                                       !arclength coordinate of the first marker
           sc(1)=s                                                      !store the arclength coordinate of a point to an array
           write (12,*) x(1), y(1)
           do i=2, N
               s=s+ds
               sc(i)=s
               if (s.le.s1) then
                   x(i)=x(i-1)+ds*dsin(alpha)
                   y(i)=y(i-1)+ds*dcos(alpha)
               elseif (s.le.s1+s2) then
                   x(i)=x(i-1)
                   y(i)=y(i-1)+ds
               elseif (s.le.s1+s2+s3) then 
                   
                   x(i)=a+r*dcos(pi-ss/r)
                   y(i)=b+r*dsin(pi-ss/r)
                   ss=ss+ds
               else
                   x(i)=x(i-1)+ds
                   y(i)=y(i-1)
               end if
               write (12,*) x(i), y(i)
           end do
       elseif (inclined) then
           xtoei=8.0d0
           alpha=pi/4.0d0
           sw=h0/dsin(alpha)
           dam_height=1.1d0*dble(h0)
           s1=dble(dam_height)/dsin(alpha)              ! length of the plane beach
           s2=dble(xlen)-dble(xtoei)-s1*dcos(alpha)     ! length of fam crest
           P=s1+s2
           N=max(int(P/(0.5*dxmin)),int(P/(0.5*dymin))) !number of markers
           allocate(x(N),y(N),sc(N))
           ds=P/N                                       !distance between two marker points
           x(1)=xtoei                                   !x coordinate of the first marker      
           y(1)=0.0d0                                   !y coordinate of the first marker
           s=0.d0                                       !arclength coordinate of the first marker
           sc(1)=s                                      !store the arclength coordinate of a point to an array
           write (12,*) x(1), y(1)
           do i=2, N
               s=s+ds
               sc(i)=s
               if (s.lt.s1) then
                   x(i)=x(i-1)+ds*dcos(alpha)
                   y(i)=y(i-1)+ds*dsin(alpha)
               else  
                   x(i)=x(i-1)+ds
                   y(i)=y(i-1)
               end if
               write (12,*) x(i), y(i)
           end do
       elseif (arch) then
           dxmin=dxmin/h0
           dymin=dymin/h0
           scale=1.0d0/6.0d0
           s1=6.153170d0
           sw=(s1)*scale*h0
           s1=6.66d0
           !s2=dble(xlen)/scale-(xtoe/scale+0.90d0)
           s2=dble(xlen/h0)/scale-(xtoe/scale/h0)
           s3=7.0d0*dy(1)/h0
           dmin=dmin1(dxmin,dymin)
           r=5.0d0*dmin/scale
           s4=2.0d0*pi*r/4.0d0
           !P=s1+s2+s3+s4
           P=s1+s2+s4
           N=max(int(P/(0.5*dxmin/scale)),int(P/(0.5*dymin/scale))) !number of markers
           N1=dint(s1*N/P) !number of markers at the toe of the arch dam
           N2=0            
           N=N+N2
           allocate(x(N),y(N),sc(N))
           allocate(xt(N),yt(N),st(N))
           ds=P/N                                       !distance between two marker points
           !xt(1)=0.2675d0                                   !x coordinate of the first marker      
           !yt(1)=0.0d0                                   !y coordinate of the first marker
           !s=0.d0                                       !arclength coordinate of the first marker
           !st(1)=s                                      !store the arclength coordinate of a point to an array
           s=-ds
           do i=1, N1-1
               s=s+ds
               st(i)=s
               a0=sc(i)
100            fs=4.64757d0-2.25d0*dsqrt(1.2025d0-0.09d0*a0+0.01d0*a0**2)+&
                    0.5d0*a0*dsqrt(1.2025d0-0.09d0*a0+0.01d0*a0**2)-5.0d0*dasinh(0.45d0-0.1d0*a0)-s
               fsd=dsqrt(1.2025d0-0.09d0*a0+0.01d0*a0**2)
               an=a0-fs/fsd
               res=4.64757d0-2.25d0*dsqrt(1.2025d0-0.09d0*an+0.01d0*an**2)+&
                    0.5d0*an*dsqrt(1.2025d0-0.09d0*an+0.01d0*an**2)-5.0d0*dasinh(0.45d0-0.1d0*an)-s
               if (dabs(res)<1.d-10) then 
                   goto 110
               else 
                   a0=an
                   goto 100
               end if 
110            continue
               yt(i)=an+s3                                          !y coordinate of the marker
               xt(i)=1.28d0-0.05d0*(an-4.5d0)**2                    !x coordinate of the marker
           end do
           !do i=N1,N1+N2
           !    s=s+ds
           !    yt(i)=yt(i-1)+ds+s3
           !    xt(i)=xt(i-1)
           !    st(i)=s
           !end do
           N1=N1+N2
           do i=1, N1
               xt(i)=xtoe/h0-xt(i)*scale
               yt(i)=(yt(N1-1)-yt(i))*scale
               st(i)=st(i)*scale
           end do
           do i=1, N1-1
               sc(i)=st(i)
               x(i)=xt(N1-i)
               y(i)=yt(N1-i)
           end do
           deallocate(xt,yt,st)
           ds=ds*scale
           r=r*scale
           a=x(N1-1)+r
           b=y(N1-1)
           do i=N1, N
               sc(i)=sc(i-1)+ds
               if (sc(i).lt.(s1+s3+s4)*scale) then 
                   x(i)=a+r*dcos(pi-ss/r)
                   y(i)=b+r*dsin(pi-ss/r)
                   ss=ss+ds
               else
                   x(i)=x(i-1)+ds
                   y(i)=y(i-1)
               end if 
           end do
           do i=1, N
               x(i)=x(i)*h0
               y(i)=y(i)*h0
               sc(i)=sc(i)*h0
               write (12,*) x(i), y(i)
           end do
       else
           xtoei=9.0d0
           sw=h0
           dam_height=1.05d0*dble(h0)
           s1=dam_height
           dmin=dmin1(dxmin,dymin)
           r=5.0d0*dmin
           a=dble(xtoe)+r
           b=s1
           s2=2.0d0*pi*r/4.0d0     !length of fam crest
           s3=dble(xlen)-dble(xtoe)-r
           P=s1+s2+s3
           N=max(int(P/dxmin),int(P/dymin))  !number of markers
           allocate(x(N),y(N),sc(N))            
           ds=P/N                              !distance between two marker points
           x(1)=dble(xtoe)                     !x coordinate of the first marker      
           y(1)=0.0d0                          !y coordinate of the first marker
           s=0.d0                              !arclength coordinate of the first marker
           sc(1)=s                             !store the arclength coordinate of a point to an array
           ss=0.0d0
           write (12,*) x(1), y(1)
           do i=2, N
               s=s+ds
               sc(i)=s
               
               if (s.le.s1) then
                   x(i)=x(i-1)
                   y(i)=y(i-1)+ds
               elseif (s.le.s1+s2) then 
                   ss=ss+(sc(i)-sc(i-1))
                   x(i)=a+r*dcos(pi-ss/r)
                   y(i)=b+r*dsin(pi-ss/r)
                   
               else
                   x(i)=x(i-1)+ds
                   y(i)=y(i-1)
               end if
               

               write (12,*) x(i), y(i)
           end do
       end if
       
   elseif (problem=="wavef") then 
       m=1.0/6.0          !seaward slope 
       nn=1.0/3.0          !landward slope
       alpha=atan(m)
       beta=atan(nn)
       d=0.8d0  		    !crest height
       s1=d/dsin(alpha)     !length of the plane beach
       s2=0.3d0             !crest width
       s3=0.2d0/dcos(beta)           
       s4=d-0.2d0*nn
       s5=4.8d0+0.5d0
       P=s1+s2+s3+s4+s5
       N=max(int(P/dxmin),int(P/dymin))    !number of markers
       allocate(x(N),y(N),sc(N))
       ds=P/N                              !distance between two marker points
       x(1)=1.0d0                          !x coordinate of the first marker      
       y(1)=0.0                             !y coordinate of the first marker
       s=0.d0                              !arclength coordinate of the first marker
       sc(1)=s                             !store the arclength coordinate of a point to an array
       if (t==0.0) write (10,*) x(1), y(1)
       do i=2, N
           s=s+ds
           sc(i)=s
           if (s.lt.s1) then
               x(i)=x(i-1)+ds*dcos(alpha)
               y(i)=y(i-1)+ds*dsin(alpha)
           elseif (s.lt.s1+s2) then 
               x(i)=x(i-1)+ds
               y(i)=y(i-1)
           elseif (s.lt.s1+s2+s3) then 
               x(i)=x(i-1)+ds*dcos(beta)
               y(i)=y(i-1)-ds*dsin(beta)
           elseif (s.lt.s1+s2+s3+s4) then 
               x(i)=x(i-1)
               y(i)=y(i-1)-ds
           else
               x(i)=x(i-1)-ds
               y(i)=y(i-1)
           end if
           write (10,*) x(i), y(i)
       end do

   end if
   
   write (11,*) 'VARIABLES= ','"X",','"Z"'
   write (11,*) 'ZONE i=',Mx+1,',j=',My+1
   zz=-dy(0)
   do j=1, My+1
       zz = zz+ dy(j-1)
       xx=-dx(0)
       do i=1, Mx+1
           xx=xx+dx(i-1)
           write(11,*) xx, zz
           
       end do
   end do
    

    return 

    end subroutine ARCGEN


    subroutine ASSIGNINTERNALBOUND(problem, imax, jmax, flag)
    use nrtype
    use defs
    implicit none
    character (LEN=30), INTENT(IN) :: problem
    integer, intent(IN) :: imax, jmax
    integer, dimension(0:,0:), intent(INOUT) :: flag

    integer i, j, k, t, jib
    logical ib, zero


    !identify all u-velocity points inside objects with 1
    do i=1, imax
        do j=1, jmax
            if (flag(i,j)==1) then
                !check whether there is a solid cell at the north of the cell
                if (problem=="reservoir".or.problem=="solitary") then 
                    do k=i, imax
                        flag(k,j)=1
                    end do
                else
                    ib=.false.
                    zero=.false.
                    do k=j+1, jmax
                        if (flag(i,k)==0) zero=.true.
                        if (flag(i,k)==1) then 
                            ib=.true.
                            exit
                        end if
                    end do
                    k=j+1
                    do while (zero.and.ib) 
                        flag(i,k)=1
                        if (flag(i,k+1)==0) exit
                        k=k+1
                    end do
                    !check whether there is a solid cell at the south of the cell
                    ib=.false.
                    zero=.false.
                    do k=j-1, 1, -1
                        if (flag(i,k)==0) zero=.true.
                        if (flag(i,k)==1) then 
                            ib=.true.
                            exit
                        end if
                    end do
                    k=j-1
                    do while (zero.and.ib) 
                        flag(i,k)=1
                        if (flag(i,k-1)==0) exit
                        k=k-1
                    end do
                end if

            end if
        end do
    end do

    return 
    end subroutine ASSIGNINTERNALBOUND
    
    
    
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FINDPOINT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    integer function FINDPOINT(x, vp)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use defs
    implicit none
    real(DP), INTENT(IN):: vp
    real(DP), DIMENSION(0:), INTENT(IN):: x
    !
    !-----------------------------------------------------------------------
    !
    integer:: k
    logical found

    found=.FALSE.
    k=0
    do while(.NOT.(found))
        k=k+1
        if (x(k)>vp.AND.x(k-1)<=vp) then
            FINDPOINT=k
            found=.TRUE.
            return
        end if 
    end do

    return
    end function FINDPOINT
    
    !!
    !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FINDFLUIDX %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !!
    !integer function FINDFLUIDX(flag, i, j, const)
    !!
    !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !!
    !! Find a fluid cell in x direction and return the i index of the fluid cell 
    !!
    !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !!
    !use nrtype
    !use defs
    !implicit none
    !integer, intent(IN) :: i, j, const
    !integer(I2B), dimension(0:,0:), intent(IN) :: flag
    !!
    !!-----------------------------------------------------------------------
    !!
    !integer:: k
    !logical found
    !
    !found=.FALSE.
    !k=i
    !do while(.NOT.(found))
    !    k=k+const
    !    if (IAND(flag(k,j),  C_F)/=0) then
    !        FINDFLUIDX=k
    !        found=.TRUE.
    !        return
    !    end if 
    !end do
    !
    !return
    !end function FINDFLUIDX
    !
    !!
    !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FINDFLUIDX %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !!
    !integer function FINDFLUIDY(flag, i, j, const)
    !!
    !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !!
    !! Find a fluid cell in y direction and return the j index of the fluid cell 
    !!
    !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !!
    !use nrtype
    !use defs
    !implicit none
    !integer, intent(IN) :: i, j, const
    !integer(I2B), dimension(0:,0:), intent(IN) :: flag
    !!
    !!-----------------------------------------------------------------------
    !!
    !integer:: k
    !logical found
    !
    !found=.FALSE.
    !k=j
    !do while(.NOT.(found))
    !    k=k+const
    !    !if the cell is active cell
    !    if (IAND(flag(i,k),  C_F)/=0) then
    !        FINDFLUIDY=k
    !        found=.TRUE.
    !        return
    !    end if 
    !end do
    !
    !return
    !end function FINDFLUIDY
    
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FINDFLUIDX %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    integer function FINDFLUIDX(imax, ac, i, j, const)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    ! Find a fluid cell in x direction and return the i index of the fluid cell 
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use defs
    implicit none
    integer, intent(IN) :: imax, i, j, const
    integer, dimension(:,:), intent(IN) :: ac
    !
    !-----------------------------------------------------------------------
    !
    integer:: k
    logical found
    
    found=.FALSE.
    k=i
    do while(.NOT.(found))
        k=k+const        
        if (ac(k,j)==1) then
            FINDFLUIDX=k
            found=.TRUE.
            return
        elseif (k==1.or.k==imax) then 
            FINDFLUIDX=i
            return
        end if 
    end do
    
    return
    end function FINDFLUIDX
    
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FINDFLUIDX %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    integer function FINDFLUIDY(jmax, ac, i, j, const)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    ! Find a fluid cell in y direction and return the j index of the fluid cell 
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use defs
    implicit none
    integer, intent(IN) :: jmax, i, j, const
    integer, dimension(:,:), intent(IN) :: ac
    !
    !-----------------------------------------------------------------------
    !
    integer:: k
    logical found
    
    found=.FALSE.
    k=j
    do while(.NOT.(found))
        k=k+const
        !if the cell is active cell
        if (ac(i,k)==1) then
            FINDFLUIDY=k
            found=.TRUE.
            return
        elseif (k==1.or.k==jmax) then 
            FINDFLUIDY=j
            return
        end if 
    end do
    
    return
    end function FINDFLUIDY



    subroutine MATRIXINVERSA(pp, r, a)
    use nrtype
    implicit none
    integer, intent(IN) :: pp
    real(DP), dimension(:,:), intent(IN) :: r
    real(DP), dimension(:,:,:), intent(OUT) :: a
    a(pp,1,1)=(-r(2,3)*r(3,2)+r(2,2)*r(3,3))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(2,3)*r(3,2)-r(1,2)*r(3,3)+r(2,2)*r(3,3))
    a(pp,1,2)=(r(1,3)*r(3,2)-r(1,2)*r(3,3))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(2,3)*r(3,2)-r(1,2)*r(3,3)+r(2,2)*r(3,3))
    a(pp,1,3)=(-r(1,3)*r(2,2)+r(1,2)*r(2,3))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(2,3)*r(3,2)-r(1,2)*r(3,3)+r(2,2)*r(3,3))
    a(pp,2,1)=(r(2,3)-r(3,3))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(2,3)*r(3,2)-r(1,2)*r(3,3)+r(2,2)*r(3,3))
    a(pp,2,2)=(-r(1,3)+r(3,3))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(2,3)*r(3,2)-r(1,2)*r(3,3)+r(2,2)*r(3,3))
    a(pp,2,3)=(r(1,3)-r(2,3))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(2,3)*r(3,2)-r(1,2)*r(3,3)+r(2,2)*r(3,3))
    a(pp,3,1)=(-r(2,2)+r(3,2))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(2,3)*r(3,2)-r(1,2)*r(3,3)+r(2,2)*r(3,3))
    a(pp,3,2)=(r(1,2)-r(3,2))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(2,3)*r(3,2)-r(1,2)*r(3,3)+r(2,2)*r(3,3))
    a(pp,3,3)=(-r(1,2)+r(2,2))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(2,3)*r(3,2)-r(1,2)*r(3,3)+r(2,2)*r(3,3))
    
    return
    end subroutine MATRIXINVERSA
    
    subroutine MATRIXINVERSB(pp, r, a)
    use nrtype
    implicit none
    integer, intent(IN) :: pp
    real(DP), dimension(:,:), intent(IN) :: r
    real(DP), dimension(:,:,:), intent(OUT) :: a
    a(pp,1,1)=(-r(2,3)*r(3,2)+r(2,2)*r(3,3))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(1,2)*r(3,3))
    a(pp,1,2)=(r(1,3)*r(3,2)-r(1,2)*r(3,3))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(1,2)*r(3,3))
    a(pp,1,3)=(-r(1,3)*r(2,2)+r(1,2)*r(2,3))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(1,2)*r(3,3))
    a(pp,2,1)=(r(2,3)-r(3,3))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(1,2)*r(3,3))
    a(pp,2,2)=(-r(1,3))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(1,2)*r(3,3))
    a(pp,2,3)=(r(1,3))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(1,2)*r(3,3))
    a(pp,3,1)=(-r(2,2)+r(3,2))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(1,2)*r(3,3))
    a(pp,3,2)=(r(1,2))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(1,2)*r(3,3))
    a(pp,3,3)=(-r(1,2))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(1,2)*r(3,3))
    
    return
    end subroutine MATRIXINVERSB
    
    subroutine COMP_SURFACE_FORCE(problem, imax, jmax, nm, ppp, ppu, ppv, invagp, invagu, invagv, gFIp, gFIu, gFIv, &
                                    indp, indu, indv, ghp, ghu, ghv, P, U, V, s, Re, xm, ym, t, au, bu, cu, av, bv, cv, &
                                    ap, bp, cp, uxb, uyb, ro, GZ, h0, alphag, sw, pflag, dx, dy, xc, yc)
    use defs
    implicit none
    character (LEN=30), intent(IN) :: problem
    integer, INTENT(IN) :: imax, jmax, nm, ppp, ppu, ppv
    real(DP), dimension(:,:,:), intent(IN) :: invagp, invagu, invagv
    integer, dimension(:,:,:), intent(IN) :: gFIp, gFIu, gFIv
    integer, dimension(:,:), intent(IN) :: indp, indu, indv
    real(DP), dimension(:,:), intent(INOUT) :: ghp, ghu, ghv
    real(RP), dimension(0:,0:), intent(IN) :: P, U, V
    real(DP), dimension(:), intent(IN) :: s, xm, ym
    real(RP), intent(IN):: Re, t, uxb, uyb, ro, GZ, h0, alphag, sw
    real(DP), dimension(:), intent(IN) :: au, bu, cu, av, bv, cv, ap, bp, cp
    integer(I2B), dimension(0:,0:), intent(IN) :: pflag
    real(RP), dimension(0:), intent(IN) :: dx, dy 
    real(DP), dimension(0:), intent(IN) :: xc, yc
    
    
    
    integer:: i, j, i1, j1, i2, j2, im1, ip1, ic, ju, jd, imin, mm, ii, ip, iu, id, N, K, jc, i0, j0
    real(DP):: a0, a1, a2, pi, D, distu, distd, distmin, smin, scim1, &
                sci, scip1, xim1, xi, xip1, dist, fd, fl, cd, cl, su, sd, ss, distminu, &
                distmind, ds, X, DXX, ALPHA, BETA, GAMMA, ETA, dsx, dsy, &
                pm, dudxm, dudzm, dwdxm, dwdzm, nx, ny, nxp, nyp, ptoe, ps, pinf, xx, zz, &
                nxxp, nyyp, nxx, nyy, cdd, Acy, KC, Vcy, period
    real(DP), dimension(3):: fi
    real(DP), dimension(3,3) :: r,invr
    real(DP), dimension(:), allocatable:: st, spp, spu, spv, a, b, c, ao, bo, co, XII, FII, P22
    real(DP), dimension(nm):: fx, fy
    real(DP), dimension(ppp):: pp
    real(DP), dimension(ppu,2):: uu
    real(DP), dimension(ppv,2):: vv
    real(DP), dimension(nm,5):: force
    logical found
    OPEN(41,FILE='ps.plt',STATUS='UNKNOWN')
    OPEN(42,FILE='arranged.plt',STATUS='UNKNOWN')
    open(34,FILE='dragf.plt',STATUS='UNKNOWN')

    pi=4.0d0*datan(1.0d0)
    D=1.0d0     !circle diameter
    do i=1, ppp
        distmin=1.0d10
        do j=1, imax
            do k=1, jmax
                dist=dsqrt((ghp(i,3)-xc(j))**2+(ghp(i,4)-yc(k))**2)
                if ((dist.lt.distmin).and.(IAND(pflag(j,k),C_F)/=0)) then 
                    distmin=dist
                    ic=j
                    jc=k
                end if
            end do
        end do
        r(1,1)=1.0d0
        r(2,1)=1.0d0
        r(3,1)=1.0d0
        r(1,2)=xc(ic)
        r(1,3)=yc(jc)
        nxx=ghp(i,7)
        nyy=ghp(i,8)
        nxxp=DINT(nxx*1.0d6)/1.0d6
        nyyp=DINT(nyy*1.0d6)/1.0d6
        if (nxxp.ge.0.0d0) then
            i1=ic+1
        elseif (nxxp.lt.0.0d0) then 
            i1=ic-1
        end if
        j1=jc
        if (nxxp.eq.0.0d0) then 
            r(2,2)=1.0d20
            r(2,3)=1.0d20
        else
            r(2,2)=xc(i1)
            r(2,3)=yc(j1)
        end if
        if (nyyp.ge.0.0d0) then
            j2=jc+1
        elseif (nyyp.lt.0.0d0) then
            j2=jc-1
        end if 
        i2=ic
        if (nyyp.eq.0.0d0) then 
            r(3,2)=1.0d20
            r(3,3)=1.0d20
        else
            r(3,2)=xc(i2)
            r(3,3)=yc(j2)
        end if
        fi(1)=P(ic,jc)
        fi(2)=P(i1,j1)
        fi(3)=P(i2,j2)
        a0=(fi(3)*(-r(1,3)*r(2,2)+r(1,2)*r(2,3)))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(2,3)*r(3,2)-& 
            r(1,2)*r(3,3)+r(2,2)*r(3,3))+(fi(2)*(r(1,3)*r(3,2)-r(1,2)*r(3,3)))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+&
            r(1,3)*r(3,2)-r(2,3)*r(3,2)- & 
            r(1,2)*r(3,3)+r(2,2)*r(3,3))+(fi(1)*(-r(2,3)*r(3,2)+r(2,2)*r(3,3)))/(-r(1,3)*r(2,2)+ r(1,2)*r(2,3)+&
            r(1,3)*r(3,2)-r(2,3)*r(3,2) -  &
            r(1,2)*r(3,3)+r(2,2)*r(3,3))
        a1=(fi(3)*(r(1,3)-r(2,3)))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(2,3)*r(3,2)-r(1,2)*r(3,3)+ & 
            r(2,2)*r(3,3))+(fi(1)*(r(2,3)-r(3,3)))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(2,3)*r(3,2)-&
            r(1,2)*r(3,3) + & 
            r(2,2)*r(3,3))+(fi(2)*(-r(1,3)+r(3,3)))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(2,3)*r(3,2)-&
            r(1,2)*r(3,3) + r(2,2)*r(3,3))
        a2=(fi(3)*(-r(1,2)+r(2,2)))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(2,3)*r(3,2)-r(1,2)*r(3,3)+ &
            r(2,2)*r(3,3))+(fi(2)*(r(1,2)-r(3,2)))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-r(2,3)*r(3,2)-&
            r(1,2)*r(3,3) + &
            r(2,2)*r(3,3))+(fi(1)*(-r(2,2)+r(3,2)))/(-r(1,3)*r(2,2)+r(1,2)*r(2,3)+r(1,3)*r(3,2)-&
            r(2,3)*r(3,2)-r(1,2)*r(3,3)+r(2,2)*r(3,3))
        !compute u variable at ghost point using linear interpolation
        ghp(i,6)=a0+a1*ghp(i,3)+a2*ghp(i,4)
    end do
    
    do i=1, ppu
        i1=gFIu(i,1,1)          !i index of first interpolation point
        j1=gFIu(i,1,2)          !j index of first interpolation point
        i2=gFIu(i,2,1)          !i index of second interpolation point
        j2=gFIu(i,2,2)          !j index of second interpolation point
        fi(1)=uxb               !horizontal solid velocity 
        fi(2)=U(i1,j1)          !velocity at first interpolation point
        fi(3)=U(i2,j2)          !velocity at second interpolation point
        !compute linear interpolation constants a1 and a2
        a1=fi(1)+fi(2)*invagu(i,2,2)+fi(3)*invagu(i,2,3)
        a2=fi(1)+fi(2)*invagu(i,3,2)+fi(3)*invagu(i,3,3)
        !compute du/dx variable at intersection point 
        ghu(i,6)=a1
        !compute du/dz variable at intersection point 
        ghu(i,7)=a2        
    end do  
    do i=1, ppv
        i1=gFIv(i,1,1)          !i index of first interpolation point
        j1=gFIv(i,1,2)          !j index of first interpolation point
        i2=gFIv(i,2,1)          !i index of second interpolation point
        j2=gFIv(i,2,2)          !j index of second interpolation point
        fi(1)=uyb               !vertical solid velocity 
        fi(2)=V(i1,j1)          !velocity at first interpolation point
        fi(3)=V(i2,j2)          !velocity at second interpolation point
        !compute linear interpolation constants a1 and a2
        a1=fi(1)+fi(2)*invagv(i,2,2)+fi(3)*invagv(i,2,3)
        a2=fi(1)+fi(2)*invagv(i,3,2)+fi(3)*invagv(i,3,3)
        !compute dw/dx variable at intersection point 
        ghv(i,6)=a1
        !compute dw/dz variable at intersection point 
        ghv(i,7)=a2        
    end do
    mm=max(ppp, ppu, ppv)
    allocate (st(mm))
    allocate (spp(mm))
    allocate (spu(mm))
    allocate (spv(mm))
    allocate (a(mm))
    allocate (b(mm))
    allocate (c(mm))
    allocate (ao(mm))
    allocate (bo(mm))
    allocate (co(mm))
    !!arrange the pressure points on the solid surface according to the arclength coordinate
    !write (42,*) 'ZONE T="p"'
    do i=1, ppp
        st(i)=ghp(i,5)
    end do
    j=0
    do j=1, ppp
        !find minimum arclength coordinate
        smin=1.0d10
        do i=1, ppp
            if (st(i)<smin) then 
                smin=st(i)
                imin=i
            end if
        end do
        spp(j)=st(imin)
        pp(j)=ghp(imin,6)
        st(imin)=1.0d10
        !write (42,*) j, spp(j), pp(j)
    end do
    !write (42,*) 'ZONE T="u"'
    !arrange the u points on the solid surface according to the arclength coordinate
    do i=1, ppu
        st(i)=ghu(i,5)
    end do
    j=0
    do j=1, ppu
        !find minimum arclength coordinate
        smin=1.0d10
        do i=1, ppu
            if (st(i)<smin) then 
                smin=st(i)
                imin=i
            end if
        end do
        spu(j)=st(imin)
        uu(j,1)=ghu(imin,6)
        uu(j,2)=ghu(imin,7)
        st(imin)=1.0d10
        !write (42,*) j, spu(j)
    end do
    !write (42,*) 'ZONE T="w"'
    !arrange the v points on the solid surface according to the arclength coordinate
    do i=1, ppv
        st(i)=ghv(i,5)
    end do
    j=0
    do j=1, ppv
        !find minimum arclength coordinate
        smin=1.0d10
        do i=1, ppu
            if (st(i)<smin) then 
                smin=st(i)
                imin=i
            end if
        end do
        spv(j)=st(imin)
        vv(j,1)=ghv(imin,6)
        vv(j,2)=ghv(imin,7)
        st(imin)=1.0d10
        !write (42,*) j, spv(j)
    end do  
    !compute a,b and c polynomial coefficients for the interpolation at marker point for pressure
    do i=1, ppp
        if (i==1) then      !use forward points
            im1=i
            ii=i+1
            ip1=i+2
        else if(i==ppp) then  !use backward points
            ip1=i
            ii=i-1
            im1=i-2
        else                
            ip1=i+1
            ii=i
            im1=i-1
        end if
        scip1=spp(ip1)
        sci=spp(ii)
        scim1=spp(im1)
        xi=pp(ii)
        xim1=pp(im1)
        xip1=pp(ip1)
        !coefficients for x polynomial        
        a(i)=(scip1*(xi-xim1)+sci*(xim1-xip1)+scim1*(xip1-xi))/&
              (scim1-sci)/(scim1-scip1)/(sci-scip1)
        b(i)=(scip1**2*(xim1-xi)+scim1**2*(xi-xip1)+&
              sci**2*(xip1-xim1))/(scim1-sci)/(scim1-scip1)/(sci-scip1)
        c(i)=(scim1*scip1*(scip1-scim1)*xi+sci**2*(scip1*xim1-scim1*xip1)+&
              sci*(xip1*scim1**2-xim1*scip1**2))/(scim1-sci)/(scim1-scip1)/(sci-scip1)
    end do
    if (problem=="circle") then 
        pinf=P(imax,jmax/2) !free-stream pressure
        if (movbound) pinf=0.0
        ptoe=1.0d0
    end if
    do i=1, nm
        distmin=1.0d10
        do j=1, ppp
            dist=dabs(s(i)-spp(j))
            if (dist.lt.distmin) then 
                distmin=dist
                ic=j
            end if
        end do
        !find static pressure at the marker point for reservoir problem
        if (problem=="reservoir") then 
            ps=ro*(-GZ)*(h0-ym(i))
            pinf=ps
        end if
        !hydrodynamic pressure at the marker point
        force(i,1)=a(ic)*s(i)**2+b(ic)*s(i)+c(ic)-pinf
    end do
    
    !compute a,b and c polynomial coefficients for the interpolation at marker point for u velocity
    do i=1, ppu
        if (i==1) then      !use forward points
            im1=i
            ii=i+1
            ip1=i+2
        else if(i==ppu) then  !use backward points
            ip1=i
            ii=i-1
            im1=i-2
        else                
            ip1=i+1
            ii=i
            im1=i-1
        end if
        scip1=spu(ip1)
        sci=spu(ii)
        scim1=spu(im1)
        xi=uu(ii,1)
        xim1=uu(im1,1)
        xip1=uu(ip1,1)       
        a(i)=(scip1*(xi-xim1)+sci*(xim1-xip1)+scim1*(xip1-xi))/&
              (scim1-sci)/(scim1-scip1)/(sci-scip1)
        b(i)=(scip1**2*(xim1-xi)+scim1**2*(xi-xip1)+&
              sci**2*(xip1-xim1))/(scim1-sci)/(scim1-scip1)/(sci-scip1)
        c(i)=(scim1*scip1*(scip1-scim1)*xi+sci**2*(scip1*xim1-scim1*xip1)+&
              sci*(xip1*scim1**2-xim1*scip1**2))/(scim1-sci)/(scim1-scip1)/(sci-scip1)
        xi=uu(ii,2)
        xim1=uu(im1,2)
        xip1=uu(ip1,2)      
        ao(i)=(scip1*(xi-xim1)+sci*(xim1-xip1)+scim1*(xip1-xi))/&
              (scim1-sci)/(scim1-scip1)/(sci-scip1)
        bo(i)=(scip1**2*(xim1-xi)+scim1**2*(xi-xip1)+&
              sci**2*(xip1-xim1))/(scim1-sci)/(scim1-scip1)/(sci-scip1)
        co(i)=(scim1*scip1*(scip1-scim1)*xi+sci**2*(scip1*xim1-scim1*xip1)+&
              sci*(xip1*scim1**2-xim1*scip1**2))/(scim1-sci)/(scim1-scip1)/(sci-scip1)
    end do
    do i=1, nm
        distmin=1.0d10
        do j=1, ppu
            dist=dabs(s(i)-spu(j))
            if (dist.lt.distmin) then 
                distmin=dist
                ic=j
            end if
        end do
        !du/dx at the marker point
        force(i,2)=a(ic)*s(i)**2+b(ic)*s(i)+c(ic)
        !du/dz at the marker point
        force(i,3)=ao(ic)*s(i)**2+bo(ic)*s(i)+co(ic)
    end do
    !compute a,b and c polynomial coefficients for the interpolation at marker point for v velocity
    do i=1, ppv
        if (i==1) then      !use forward points
            im1=i
            ii=i+1
            ip1=i+2
        else if(i==ppv) then  !use backward points
            ip1=i
            ii=i-1
            im1=i-2
        else                
            ip1=i+1
            ii=i
            im1=i-1
        end if
        scip1=spv(ip1)
        sci=spv(ii)
        scim1=spv(im1)
        xi=vv(ii,1)
        xim1=vv(im1,1)
        xip1=vv(ip1,1)       
        a(i)=(scip1*(xi-xim1)+sci*(xim1-xip1)+scim1*(xip1-xi))/&
              (scim1-sci)/(scim1-scip1)/(sci-scip1)
        b(i)=(scip1**2*(xim1-xi)+scim1**2*(xi-xip1)+&
              sci**2*(xip1-xim1))/(scim1-sci)/(scim1-scip1)/(sci-scip1)
        c(i)=(scim1*scip1*(scip1-scim1)*xi+sci**2*(scip1*xim1-scim1*xip1)+&
              sci*(xip1*scim1**2-xim1*scip1**2))/(scim1-sci)/(scim1-scip1)/(sci-scip1)
        xi=vv(ii,2)
        xim1=vv(im1,2)
        xip1=vv(ip1,2)      
        ao(i)=(scip1*(xi-xim1)+sci*(xim1-xip1)+scim1*(xip1-xi))/&
              (scim1-sci)/(scim1-scip1)/(sci-scip1)
        bo(i)=(scip1**2*(xim1-xi)+scim1**2*(xi-xip1)+&
              sci**2*(xip1-xim1))/(scim1-sci)/(scim1-scip1)/(sci-scip1)
        co(i)=(scim1*scip1*(scip1-scim1)*xi+sci**2*(scip1*xim1-scim1*xip1)+&
              sci*(xip1*scim1**2-xim1*scip1**2))/(scim1-sci)/(scim1-scip1)/(sci-scip1)
    end do
    do i=1, nm
        distmin=1.0d10
        do j=1, ppu
            dist=dabs(s(i)-spv(j))
            if (dist.lt.distmin) then 
                distmin=dist
                ic=j
            end if
        end do
        !dw/dx at the marker point
        force(i,4)=a(ic)*s(i)**2+b(ic)*s(i)+c(ic)
        !dw/dz at the marker point
        force(i,5)=ao(ic)*s(i)**2+bo(ic)*s(i)+co(ic)
    end do
    cd=0.0d0
    cl=0.0d0
    do i=1, nm
        ip=i+1
        if (i==nm) ip=1 !for the case of circular cylinder (closed solid boundary)
        !find the nearest intersection point for the pressure point to determine the nx and ny pointing out from the 
        !boundary to the fluid
        distmin=1.0d10
        do j=1, ppp
            dist=dabs(ghp(j,5)-s(i))
            if (dist.lt.distmin) then 
                distmin=dist
                imin=j
            end if
        end do
        nxp=ghp(imin,7)
        nyp=ghp(imin,8)
        ds=dsqrt((ym(ip)-ym(i))**2+(xm(ip)-xm(i))**2)
        ds=dmax1(ds,1.0d-20)
        nx=dsign(dabs(ym(ip)-ym(i))/ds,nxp)
        ny=dsign(dabs(xm(ip)-xm(i))/ds,nyp)
        dsx=ds*ny
        dsy=ds*nx
        pm=(force(i,1)+force(ip,1))/2.0d0
        dudxm=(force(i,2)+force(ip,2))/2.0d0
        dudzm=(force(i,3)+force(ip,3))/2.0d0
        dwdxm=(force(i,4)+force(ip,4))/2.0d0
        dwdzm=(force(i,5)+force(ip,5))/2.0d0 
        
        if (problem=="reservoir".and.s(i).gt.sw) then 
            continue
        else
            cd=cd+(-pm+2.0d0/Re*dudxm)*dsy+1.0d0/Re*(dudzm+dwdxm)*dsx
            cl=cl+(-pm+2.0d0/Re*dwdzm)*dsx+1.0d0/Re*(dudzm+dwdxm)*dsy
        end if
    end do
    if (problem=="circle") then
        if (movbound) then 
            !KC=5.0d0
            !Acy=KC*D/2.0d0/pi
            !period=1.0d0                            !period of the oscillating cylinder
            !Vcy=Acy*2.0d0*pi/period
            !cd=cd/(0.5d0*1.0d0*Vcy**2)
            !cl=cl/(0.5d0*1.0d0*Vcy**2)
            !cd=2.0*cd
            !cl=2.0*cl
        else
            cd=2.0d0*cd
            cl=2.0d0*cl
        end if
    elseif (problem=="reservoir") then 
        cd=cd/(ro*alphag*(-GZ)*h0**2)
        cl=cl/(ro*alphag*(-GZ)*h0**2)
    end if
    cdd=dsqrt(cd**2+cl**2)
    write (34,*) t, cd, cl
    !
    !write (41,*) 'ZONE T="pressure"'
    !write (41,*) 'ZONE T="',t,'"'
    !if (problem=="reservoir") then 
    !    ptoe=ro*(-GZ)*h0*alphag
    !    do i=1, ppp
    !        ps=ro*(-GZ)*(h0-ghp(i,4))
    !        write (41,*) (ghp(i,6)-ps)/ptoe, ghp(i,4)
    !    end do
    !elseif (problem=="circle") then 
    !    alpha=(pi*D-s(i))*360.0d0/pi/D+180.0d0
    !    write (41,*) alpha, 2.0d0*force(i,1)
    !end if 
    !write (41,*) 'ZONE T="pressure"'
    write (41,*) 'ZONE T=', '"',t,'"'
    if (problem=="reservoir") then 
        ptoe=ro*(-GZ)*h0*alphag
        do i=1, nm
            if (problem=="reservoir".and.s(i).gt.sw) then 
                continue
            else
                write (41,*) force(i,1)/ptoe, ym(i)/h0
            end if
            
        end do
    elseif (problem=="circle") then
        do i=1, nm
            alpha=(pi*D-s(i))*360.0d0/pi/D+180.0d0
            write (41,*) alpha, 2.0d0*force(i,1)
        end do
    end if 
    !write (41,*) 'ZONE T="vorticity"'
    do i=1, nm
        if (problem=="circle") then 
            alpha=(pi*D-s(i))*360.0d0/pi/D+180.0d0
            write (41,*) alpha, force(i,4)-force(i,3)
        elseif (problem=="reservoir") then
            !if (problem=="reservoir".and.s(i).gt.sw) then 
            !    continue
            !else
            !    write (41,*) s(i), force(i,4)-force(i,3)
            !end if
            
        end if
        
    end  do
    !
    !write (41,*) 'ZONE T="dw/dx"'
    ! do i=1, nm
    !    alpha=(pi*D-s(i))*360.0d0/pi/D+180.0d0
    !    write (41,*) alpha, force(i,4)
    ! end  do
    ! 
    ! write (41,*) 'ZONE T="du/dz"'
    ! do i=1, nm
    !    alpha=(pi*D-s(i))*360.0d0/pi/D+180.0d0
    !    write (41,*) alpha, force(i,3)
    !end do
    !
  
    
    return     
    end subroutine COMP_SURFACE_FORCE
                                    
    subroutine CUBIC_SPLINE (N, XI, FI, P2)
    use defs
    implicit none
!
! Function to carry out the cubic-spline approximation
! with the second-order derivatives returned.
!
  integer :: I
  integer, intent (IN) :: N
  real(DP), intent (IN), dimension (N+1):: XI, FI
  real(DP), intent (OUT), dimension (N+1):: P2
  real(DP), dimension (N):: G, H
  real(DP), dimension (N-1):: D, B, C
!
! Assign the intervals and function differences
!
  do I = 1, N
      H(I) = XI(I+1) - XI(I)
      G(I) = FI(I+1) - FI(I)
  end do
!
! Evaluate the coefficient matrix elements
  do I = 1, N-1
      D(I) = 2*(H(I+1)+H(I))
      B(I) = 6*(G(I+1)/H(I+1)-G(I)/H(I))
      C(I) = H(I+1)
  end do
!
! Obtain the second-order derivatives
!
  call TRIDIAGONAL_LINEAR_EQ (N-1, D, C, C, B, G)
  P2(1) = 0
  P2(N+1) = 0
  do I = 2, N 
      P2(I) = G(I-1)
  end do
  
  return
  end subroutine CUBIC_SPLINE
!
subroutine TRIDIAGONAL_LINEAR_EQ (L, D, E, C, B, Z)
    use defs
    implicit none
!
! Functione to solve the tridiagonal linear equation set.
!
  integer, intent (IN) :: L
  integer :: I
  real(DP), intent (IN), dimension (L):: D, E, C, B
  real(DP), intent (OUT), dimension (L):: Z
  real(DP), dimension (L):: Y, W
  real(DP), dimension (L-1):: V, T
!
! Evaluate the elements in the LU decomposition
!
  W(1) = D(1)
  V(1)  = C(1)
  T(1)  = E(1)/W(1)
  do I = 2, L - 1
      W(I) = D(I)-V(I-1)*T(I-1)
      V(I) = C(I)
      T(I) = E(I)/W(I)
  end do
  W(L) = D(L)-V(L-1)*T(L-1)
!
! Forward substitution to obtain y
!
  Y(1) = B(1)/W(1)
  do I = 2, L
      Y(I) = (B(I)-V(I-1)*Y(I-1))/W(I)
  end do
!
! Backward substitution to obtain z
  Z(L) = Y(L)
  do I = L-1, 1, -1
      Z(I) = Y(I) - T(I)*Z(I+1)
  end do
  
  return  
end subroutine TRIDIAGONAL_LINEAR_EQ

    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%% WR_CIRCLE_POS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine WR_CIRCLE_POS(t, N, xlen, zlen, step)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  Computes stream function and writes to a file
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    USE nrtype
    use defs
    implicit none
    integer, intent(IN)::  N
    integer, intent(INOUT):: step
    real(RP), intent(IN):: t, xlen, zlen
    !
    !-----------------------------------------------------------------------
    !
    integer:: i, j
    real (DP):: A, pi, D, KC, period, disp, r, s, P, ds, xc, yc
    real(DP), dimension(N) :: xm, ym, sc
    
    KC=5.d0
    D=1.0d0
    pi=4.0d0*datan(1.0d0)
    A=KC*D/2.0d0/pi
    period=1.0d0                            !period of the oscillating cylinder
    disp=-A*dsin(2.0d0*pi*t/period)         !displacement of the cylinder
    xc=xlen/2.0+disp                           !x coordinate of center
    yc=zlen/2.0                                !y coordinate of center
    !xc=10.0+disp                           !x coordinate of center
    !yc=15.0                                !y coordinate of center
    r=0.5d0                                 !radius of circle
    P=2.d0*pi*r                             !perimeter of circle
    ds=P/N                                  !distance between two marker points
    xm(1)=xc-r                               !x coordinate of the first marker
    ym(1)=yc                                 !y coordinate of the first marker
    s=0.0d0                                 !arclength coordinate of the first marker
    sc(1)=s                                 !store the arclength coordinate of a point to an array
    do i=2, N
        s=s+ds
        sc(i)=s                             ! arclength coordinate
        xm(i)=xc+r*dcos(pi-s/r)              ! x coordinate of the marker
        ym(i)=yc+r*dsin(pi-s/r)              ! y coordinate of the marker
    end do

    
    !write (4,*)	'ZONE T="',t,'"'
    !write (4,*) 'GEOMETRY X=',xc, 'Y=',yc, 'T=CIRCLE, C=RED, FC=BLUE'
    !write (4,*) '10'
    !do i=1, N
    !    write (4,*) xm(i), ym(i), sc(i)
    !end do
    
    if (t.ge.200.0) then
        step=step+1
        open(4,FILE='stream.plt')
        write (4,*) 'GEOMETRY X=', xc, 'Y=',yc, 'T=CIRCLE, C=BLACK, ZN=', step,'FC=BLACK,CS=GRID 0.5'
        open(9,FILE='uwp.plt')
        write (9,*) 'GEOMETRY X=', xc, 'Y=',yc, 'T=CIRCLE, C=BLACK, ZN=', step,'FC=BLACK,CS=GRID 0.5'
    end if
    
    
    return
    end subroutine WR_CIRCLE_POS
                                    
                                    