!  Regula.f90 
!
!  FUNCTIONS:
!  Regula - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Regula
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    !program Regula
    !
    !implicit none
    !integer:: iFind, iFound, i
    !double precision :: a, b  
    ! double precision, dimension(:), allocatable:: root
    !a = -5.d0
    !b = 5.d0
    !iFind=10
    !
    !call RegulaFalsi(a, b, root, iFind, iFound)
    !do i=1, iFound 
    !    print *, i, root(i)
    !end do
    !
    !
    !
    !! Variables
    !
    !! Body of Regula
    !print *, 'Hello World'
    !
    !end program Regula
    subroutine RegulaFalsi(fnf, a, b, roots, iFind, iFound, sec)
    use nrtype
    implicit none
    real(DP), intent(IN):: a, b, sec
    integer, intent(IN):: iFind
    integer, intent(OUT):: iFound
    real(DP), dimension(:), allocatable, intent(OUT) :: roots
    real(DP), external:: fnf

    ! local Variables
    integer:: iSubint, Maxhalve, iHata, k
    real(DP):: eps, D, x1, x2, f1, f2, xNew, ak, bk, aL, r, xk, fxk, fk1, epsf, iF1, iF2

    allocate (roots(iFind))
    iSubint = 200
    Maxhalve = 50
    epsf = 1.d-8
    iFound = 0
    iHata=0
    if (a>b) return
    if (a==b) iSubint = 1
    D = (b-a)/iSubint
    x2 = a
    do while (x2<b) 
        if (iFound>=iFind) return 
        x1 = x2
        x2 = x1 + D
        if (x2>b) x2 = b
        ! any root in [x1,x2] ? 
        f1 = fnf(x1, sec)
        if (dabs(f1)<epsf) then 
            xNew = x1
            iFound = iFound + 1
            roots(iFound) = xNew
            goto 6
        end if
        f2 = fnf(x2, sec)
        if (dabs(f2)<epsf) then 
            xNew = x2
            iFound = iFound + 1 
            roots(iFound)=xNew
            goto 6
        end if
        iF1=dsign(1.d0,f1)
        iF2=dsign(1.d0,f2)
        if (iF1*iF2==1.d0) then 
            goto 6
        end if

        ! iteration in the interval [ak, bk]
        ak = x1
        bk = x2
        aL = f1 
        r = f2
        xk = x1
        fxk = aL
        do k=1, Maxhalve
            if (dabs(r-aL)<epsf) then 
                xNew = ak 
            else
                xNew = (ak*r-bk*aL)/(r-aL)
            end if
            fk1 = fnf(xNew, sec)
            if (dabs(fk1)<epsf) then 
                iFound = iFound + 1
                roots(iFound) = xNew
                goto 6 
            end if
            if (fnf(ak, sec)*fk1<0.d0) then 
                bk = xNew
                r = fk1 
                if (fxk*r>0.d0) aL = 0.5d0*aL
            else
                ak = xNew 
                aL = fk1 
                if (fxk*aL>0.d0) r = 0.5d0*r
            end if
            xk = xNew 
            fxk=fk1
        end do
6       continue 
    end do
    end subroutine RegulaFalsi
    
    
    !real(DP) function fxcir(x, y)
    !use nrtype
    !implicit none
    !real(DP), intent(IN):: x, y
    !
    !fxcir=(x-10.0d0)**2+(y-15.0d0)**2-0.5d0**2
    !end function
    !
    !real(DP) function fycir(y, x)
    !use nrtype
    !implicit none
    !real(DP), intent(IN):: x, y
    !
    !fycir=(x-10.0d0)**2+(y-15.0d0)**2-0.5d0**2
    !end function
    
    real(DP) function fxcir(x,y)
    use nrtype
    implicit none
    real(DP), intent(IN):: x, y
    real(DP):: fycir
    
    fxcir=(x-10.0d0)**2+(y-15.0d0)**2-0.5d0**2
    return
    
    entry fycir(y,x)
    fycir=(x-10.0d0)**2+(y-15.0d0)**2-0.5d0**2
    return
    
    end function fxcir
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIND_INTERS_OBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine FIND_INTERS_OBS (problem, imax, jmax, nm, x, z, xm, ym, sm, nxm, nym, vos)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !   This subroutine finds intersection points between the solid object and grid.
    !   x and y coordinates of intersection points are stored to ints(:,:,:) array
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use interfaces2
    use defs
    implicit none
    character (LEN=30), INTENT(IN) :: problem
    integer, intent(IN):: imax, jmax, nm
    real(RP), dimension(0:), intent(IN):: x, z
    real(DP), dimension(:), intent(IN) :: xm, ym, sm
    real(DP), dimension(:), intent(IN) :: nxm, nym
    real(DP), dimension(0:,0:), intent(INOUT) :: vos
    !
    !local variables
    !
    integer:: iFind, iFound, iiFound, i, j, imin, im1, ii, ip1, p, reg, k, jj, jm
    real(DP):: a, b, sec, distmin, dist, px, py, scip1, sci, scim1, xi, xim1, xip1, yi, yim1, yip1, &
                s1, s2, si, minr, ds, dist1, dist2, nx, ny, xx, zz, north, south, east, west, eps, dx, dz, &
                north2, south2, west2, east2, dmin, smid, alpha, alpha2, hx, hy, c1, c2, area, smm, sc, s, nz, &
                Fc, dxc, dzc, delx, delz, hz
    !real(RP):: nx, ny
    real(DP), external:: fxcir, fycir
    real(DP), dimension(:), allocatable:: root
    real(DP), dimension(:,:), allocatable:: roots, roots2
    integer, dimension(:,:), allocatable:: cell
    !real(DP), dimension(:,:), allocatable:: ispt, isp !intersection points
    logical internal, xmid, ymid
    OPEN(43,FILE='obs.plt',STATUS='UNKNOWN')
    OPEN(44,FILE='obs2.plt',STATUS='UNKNOWN')
    OPEN(45,FILE='vos.plt',STATUS='UNKNOWN')
    !
    !find the x coordinates of intersection points between the solid and horizontal lines
    !
    allocate(roots(200,3))
    iFind=10
    a=x(1)
    b=x(imax)
    iiFound=0
    do j=2, jmax-1
        sec = z(j)
        call RegulaFalsi(fxcir, a, b, root, iFind, iFound, sec)
        do i=1, iFound
            iiFound=iiFound+1
            roots(iiFound,1)=root(i)
            roots(iiFound,2)=sec
        end do
    end do
    a=z(1)
    b=z(jmax)
    do i=2, imax-1
        sec = x(i)
        call RegulaFalsi(fycir, a, b, root, iFind, iFound, sec)
        do j=1, iFound
            iiFound=iiFound+1
            roots(iiFound,1)=sec
            roots(iiFound,2)=root(j)
        end do
    end do
    deallocate(root)
    do i=1, iiFound
        px=roots(i,1)
        py=roots(i,2)
        !find in which two marker points the intersection points is located
        do j=1, nm
            if (j.eq.nm) then 
                ds=dsqrt((xm(1)-xm(j))**2+(ym(1)-ym(j))**2)
                dist1=dsqrt((xm(1)-px)**2+(ym(1)-py)**2)
                dist2=dsqrt((px-xm(j))**2+(py-ym(j))**2)
            else
                ds=dsqrt((xm(j+1)-xm(j))**2+(ym(j+1)-ym(j))**2)
                dist1=dsqrt((xm(j+1)-px)**2+(ym(j+1)-py)**2)
                dist2=dsqrt((px-xm(j))**2+(py-ym(j))**2)
            end if
            minr=ds/20.d0   !maximum error 
            if (dabs(dist1+dist2-ds).lt.minr) then 
                ii=j
                exit
            end if
        end do
        !Find the arclength coordinate of the intersection point.
        if(ii==nm) then 
            ds=sm(ii)-sm(ii-1)
            dist1=dsqrt((xm(1)-xm(ii))**2+(ym(1)-ym(ii))**2)
            dist2=dsqrt((px-xm(ii))**2+(py-ym(ii))**2)
            si=sm(ii-1)+dist2/dist1*ds
        else
            ds=sm(ii+1)-sm(ii)
            dist1=dsqrt((xm(ii+1)-xm(ii))**2+(ym(ii+1)-ym(ii))**2)
            dist2=dsqrt((xm(ii)-px)**2+(ym(ii)-py)**2)
            si=sm(ii)+dist2/dist1*ds
        end if
        roots(i,3)=si
    end do
    
    do i=1, iiFound
        write (43,*) roots(i,1), roots(i,2), roots(i,3)
    end do
    allocate(roots2(iiFound,3))
    !sorth the intersection points with respect to arclength coordinate
    do j=1, iiFound
        minr=1.0d10
        do i=1, iiFound
            if (roots(i,3).lt.minr) then 
                ii=i
                minr=roots(i,3)
            end if
        end do
        roots2(j,1)=roots(ii,1)
        roots2(j,2)=roots(ii,2)
        roots2(j,3)=roots(ii,3)
        roots(ii,3)=1.0d10
    end do
    deallocate(roots)
    allocate(cell(iiFound,2))
    !find inwhich cell the intersection point is located
    do k=1, iiFound
        if (k==iiFound) then 
            xi=roots2(k,1)
            yi=roots2(k,2)
            xip1=roots2(1,1)
            yip1=roots2(1,2)
        else
            xi=roots2(k,1)
            yi=roots2(k,2)
            xip1=roots2(k+1,1)
            yip1=roots2(k+1,2)
        end if
        eps=0.01d0
        do i=1, imax
            dx=x(i)-x(i-1)
            do j=1, jmax
                dz=z(j)-z(j-1)
                dmin=min(abs(dx),abs(dz))
                xmid=.false.
                ymid=.false.
                if ((xi.ge.x(i-1)).and.(xi.le.x(i))) xmid=.true.
                if ((yi.ge.z(j-1)).and.(yi.le.z(j))) ymid=.true.
                ii=0
                north=dabs(z(j)-yi)/dmin
                south=dabs(z(j-1)-yi)/dmin
                east=dabs(x(i)-xi)/dmin
                west=dabs(x(i-1)-xi)/dmin
                if((north.le.eps).and.xmid) ii=ii+1
                if((south.le.eps).and.xmid) ii=ii+1
                if((east.le.eps).and.ymid) ii=ii+1
                if((west.le.eps).and.ymid) ii=ii+1
                
                xmid=.false.
                ymid=.false.
                if ((xip1.ge.x(i-1)).and.(xip1.le.x(i))) xmid=.true.
                if ((yip1.ge.z(j-1)).and.(yip1.le.z(j))) ymid=.true.
                jj=0
                north2=dabs(z(j)-yip1)/dmin
                south2=dabs(z(j-1)-yip1)/dmin
                east2=dabs(x(i)-xip1)/dmin
                west2=dabs(x(i-1)-xip1)/dmin
                if((north2.le.eps).and.xmid) jj=jj+1
                if((south2.le.eps).and.xmid) jj=jj+1
                if((east2.le.eps).and.ymid) jj=jj+1
                if((west2.le.eps).and.ymid) jj=jj+1
                if (ii==1.and.jj==1) then 
                    cell(k,1)=i
                    cell(k,2)=j
                    goto 10
                end if
            end do
        end do
10      continue
    end do
    
    do i=1, iiFound 
        !write (44,*) i, roots2(i,1), roots2(i,2), roots2(i,3), cell(i,1), cell(i,2)
        write (44,*) roots2(i,1), roots2(i,2), roots2(i,3)
        !write (44,*) i, cell(i,1), cell(i,2)
    end do
    !initialize vos values 
    do i=1, iiFound
        if (i==88) then 
            continue
        end if
        ii=cell(i,1)
        jj=cell(i,2)
        
        delx=x(ii)-x(ii-1)
        delz=z(jj)-z(jj-1)
        if (i==iiFound) then 
            ds=roots2(i,3)-roots2(i-1,3)
            smid=roots2(i,3)+ds/2.0d0
        else
            smid=(roots2(i,3)+roots2(i+1,3))/2.0d0
        end if
        !find the nearest marker point to smid
        distmin=1.0d10
        do j=1, nm
            dist=dabs(smid-sm(j))
            if (dist.lt.distmin) then 
                distmin=dist
                jm=j
            end if
        end do
        xi=dabs(roots2(i,1)-x(ii-1))
        yi=dabs(roots2(i,2)-z(jj-1))
        nx=nxm(jm)
        nz=nym(jm)
        if (nx.lt.0.0d0) xi=x(ii)-roots2(i,1)
        if (nz.lt.0.0d0) yi=z(jj)-roots2(i,2)
        nx=dabs(nx)+1.0d-50
        nz=dabs(nz)+1.0d-50
        alpha=xi*nx+yi*nz
        hx=1.0d0
        hz=1.0d0
        if (alpha-nx*delx.lt.0.0d0) hx=0.0d0
        if (alpha-nz*delz.lt.0.0d0) hz=0.0d0
        area=alpha**2/2.0d0/nx/nz*(1.0d0-hx*((alpha-nx*delx)/alpha)**2-&
            hz*((alpha-nz*delz)/alpha)**2)/(delx*delz)
        vos(ii,jj)=area
        continue
    end do
    
    call ASSIGNINTERNALVOS(problem, imax, jmax, vos)
    
    write (45,*) 'VARIABLES= ','"X",','"Z",','"vos"'
    write (45,*) 'ZONE i=',imax,',j=',jmax
    do i=1, imax
        dx=x(i)-x(i-1)
        do j=1, jmax
            dz=z(j)-z(j-1)
            write (45,*) x(i)-dx/2.0, z(j)-dz/2.0, vos(i,j)
        end do
    end do
    
    
   
    return
    end subroutine
    
    
subroutine ASSIGNINTERNALVOS(problem, imax, jmax, vos)
    use nrtype
    use defs
    implicit none
    integer, intent(IN) :: imax, jmax
    real(DP), dimension(0:,0:), intent(INOUT) :: vos
    character (LEN=30), INTENT(IN) :: problem

    
    integer i, j, k, t, jib
    logical ib, zero


    !identify all u-velocity points inside objects with 1
    do i=1, imax
        do j=1, jmax
            if (vos(i,j)/=0.0d0) then
                !check whether there is a solid cell at the north of the cell
                if (problem=="reservoir".or.problem=="solitary") then 
                    do k=i, imax
                        vos(k,j)=1.0d0
                    end do
                else
                    ib=.false.
                    zero=.false.
                    do k=j+1, jmax
                        if (vos(i,k)==0.0d0) zero=.true.
                        if (vos(i,k)/=0.0d0) then 
                            ib=.true.
                            exit
                        end if
                    end do
                    k=j+1
                    do while (zero.and.ib) 
                        vos(i,k)=1.0d0
                        if (vos(i,k+1)/=0.0d0) exit
                        k=k+1
                    end do
                    !check whether there is a solid cell at the south of the cell
                    ib=.false.
                    zero=.false.
                    do k=j-1, 1, -1
                        if (vos(i,k)==0.0d0) zero=.true.
                        if (vos(i,k)/=0.0d0) then 
                            ib=.true.
                            exit
                        end if
                    end do
                    k=j-1
                    do while (zero.and.ib) 
                        vos(i,k)=1.0d0
                        if (vos(i,k-1)==0.0d0) exit
                        k=k-1
                    end do
                end if

            end if
        end do
    end do

    
    return 
end subroutine ASSIGNINTERNALVOS    

    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HALLOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine HALLOW (problem, imax, jmax, FLAG, vos, C)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !   Subtract vos values from VOF values
    !   
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use interfaces2
    use defs
    implicit none
    character (LEN=30), INTENT(IN) :: problem
    integer, intent(IN):: imax, jmax
    integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
    real(DP), dimension(0:,0:), intent(INOUT) :: vos, C
    
    integer:: i, j
    
    do i=1, imax
        do j=1, jmax
            !substract vos values from vof values for only surface cells
            if (i==87.and.j==113) then 
                continue
            end if
            
            if (.NOT.((IAND(FLAG(i,j),C_E)/=0).OR.(FLAG(i,j)<C_O))) then
                C(i,j)=dmax1(C(i,j)-vos(i,j),0.0d0)
            end if
            
        end do
    end do
   
    return
    end subroutine HALLOW
    
    
    
    
        
    
    

