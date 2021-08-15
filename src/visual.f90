    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%% WR_DCAVITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine WR_DCAVITY(imax, jmax, U, P, z, dz, xlen, um, ro, Re)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  Compares numerical results with Ghia's results
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    implicit none
    integer, intent(IN):: imax, jmax
    real(RP), intent(IN):: xlen, um, ro, Re
    real(RP), dimension(0:), intent(IN):: z, dz
    real(RP), dimension(0:,0:), intent(IN):: U, P
    !
    !-----------------------------------------------------------------------
    !
    integer:: j
    real(RP):: dzz
    real(RP), dimension(0:129):: GU, GP, GZ
    open(3,FILE='ghia.plt',STATUS='UNKNOWN')
    dzz=1.0/128
    !Horizontal velocity component along the centerline
    if (anint(Re)==400) then
        GU(1)=0.000000
        GU(8)=-0.08186
        GU(9)=-0.09266
        GU(10)=-0.10338
        GU(14)=-0.14612
        GU(23)=-0.24299
        GU(37)=-0.32726
        GU(59)=-0.17119
        GU(65)=-0.11477
        GU(80)=0.02135
        GU(95)=0.16256
        GU(110)=0.29093
        GU(123)=0.55892
        GU(124)=0.61756
        GU(125)=0.68439
        GU(126)=0.75837
        GU(129)=1.00000
    else if (anint(Re) == 1000) then
        GU(1)=0.000000
        GU(8)=-0.18109
        GU(9)=-0.20196
        GU(10)=-0.22220
        GU(14)=-0.29730
        GU(23)=-0.38289
        GU(37)=-0.27805
        GU(59)=-0.10648
        GU(65)=-0.06080
        GU(80)=0.05702
        GU(95)=0.18719
        GU(110)=0.33304
        GU(123)=0.46604
        GU(124)=0.51117
        GU(125)=0.57492
        GU(126)=0.65928
        GU(129)=1.00000
        GP(1)=0.110591
        GP(8)=0.109689
        GP(9)=0.109200
        GP(10)=0.108566
        GP(14)=0.104187
        GP(23)=0.081925
        GP(37)=0.040377
        GP(59)=0.004434
        GP(65)=0.00000
        GP(80)=-0.000827
        GP(95)=0.012122
        GP(110)=0.034910
        GP(123)=0.050329
        GP(124)=0.050949
        GP(125)=0.051514
        GP(126)=0.052009
        GP(129)=0.052987
    elseif (anint(Re)==10000) then
        GU(1)=0.000000
        GU(8)=-0.42735
        GU(9)=-0.42537
        GU(10)=-0.41657
        GU(14)=-0.38000
        GU(23)=-0.32709
        GU(37)=-0.23186
        GU(59)=-0.07540
        GU(65)=0.03111
        GU(80)=0.08344
        GU(95)=0.20673
        GU(110)=0.34635
        GU(123)=0.47804
        GU(124)=0.48070
        GU(125)=0.47783
        GU(126)=0.47221
        GU(129)=1.00000
    end if

    do j=1, 129
        GZ(j) = (j-1)*dzz
    end do
    write (3,*) 'VARIABLES= ','"y",','"U",','"P"'
    write (3,*) 'ZONE T= ','GHIA'
    write (3,*) GZ(1),GU(1), GP(1)
    write (3,*) GZ(8),GU(8), GP(8)
    write (3,*) GZ(9),GU(9), GP(9)
    write (3,*) GZ(10),GU(10), GP(10)
    write (3,*) GZ(14),GU(14), GP(14)
    write (3,*) GZ(23),GU(23), GP(23)
    write (3,*) GZ(37),GU(37), GP(37)
    write (3,*) GZ(59),GU(59), GP(59)
    write (3,*) GZ(65),GU(65), GP(65)
    write (3,*) GZ(80),GU(80), GP(80)
    write (3,*) GZ(95),GU(95), GP(95)
    write (3,*) GZ(110),GU(110), GP(110)
    write (3,*) GZ(123),GU(123), GP(123)
    write (3,*) GZ(124),GU(124), GP(124)
    write (3,*) GZ(125),GU(125), GP(125)
    write (3,*) GZ(126),GU(126), GP(126)
    write (3,*) GZ(129),GU(129), GP(129)
    write (3,*)	'ZONE T= ','PRESENT'
    write (3,*) z(0),0.0, (P(imax/2,1)+P(imax/2+1,1))/2.0/(ro*um**2)
    !do j=1,jmax
    !    if (z(j)-dz(j)/2.0>=0.005.AND.z(j)-dz(j)/2.0<=0.015) then
            
    !    write (3,*) (z(j)-dz(j)/2.0-0.005)/0.01, U(imax/2,j)/ um, &
    !    (P(imax/2,j)+P(imax/2+1,j))/2.0/(ro*um**2)
    !    end if
        
    !end do
    !write (3,*) (z(jmax)-0.005)/0.01, (U(imax/2,jmax)+U(imax/2,jmax+1))/2.0/um, &
    !(P(imax/2,jmax)+P(imax/2+1,jmax))/2.0/(ro*um**2)
    do j=1,jmax
                    
        write (3,*) (z(j)-dz(j)/2.0), U(imax/2,j), &
        (P(imax/2,j)+P(imax/2+1,j))/2.0/(ro*um**2)

    end do
    write (3,*) (z(jmax)), (U(imax/2,jmax)+U(imax/2,jmax+1))/2.0, &
    (P(imax/2,jmax)+P(imax/2+1,jmax))/2.0/(ro*um**2)
    return
    end subroutine WR_DCAVITY
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%% WR_PSI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine WR_PSI(problem, imax, jmax, x, z, dz, U, FLAG, t)
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
    character (LEN=30), INTENT(IN) :: problem
    integer, intent(IN):: imax, jmax
    real(RP), intent(IN):: t
    integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
    real(RP), dimension(0:),intent(IN):: x, z, dz
    real(RP), dimension(0:,0:),intent(IN):: U
    !
    !-----------------------------------------------------------------------
    !
    integer i, j
    real(RP), dimension(0:imax+1,0:jmax+1):: PSI
    real(RP), dimension(0:imax):: xp
    real (DP):: A, pi, D, KC, period, disp, vc
    disp=0.0d0
    
    if (t.ge.0.0) then
        if (problem=="circle") then 
            KC=5.0d0
            D=1.0d0
            pi=4.0d0*datan(1.0d0)
            A=KC*D/2.0d0/pi
            period=1.0d0                                        !period of the oscillating cylinder
            disp=-A*dsin(2.0d0*pi*t/period)                     !displacement of the cylinder
            vc=-A*2.0d0*pi/period*dcos(2.0d0*pi*t/period)       !velocity of the cylinder
        end if
        open(4,FILE='stream.plt')
        !
        ! Computation of the stream function at the upper right corner
        ! of cell (i,j) (only if both lower cells are fluid cells)
        !
        do i = 1,imax
            PSI(i,0) = 0.0
            do j = 1,jmax
                if (((IAND(FLAG(i,j),  C_F)/=0).AND.(FLAG(i,  j)<C_E)).OR. &
                ((IAND(FLAG(i+1,j),C_F)/=0).AND.(FLAG(i+1,j)<C_E))) then
                    PSI(i,j) = PSI(i,j-1) + U(i,j)*dz(j)
                else
                    PSI(i,j) = PSI(i,j-1);
                end if
            end do
        end do
        write (4,*) 'VARIABLES= ','"X",','"Z",','"PSI"'
        write (4,*)	'ZONE T="',t,'"',', ', 'I=',imax,',J=',jmax,',   F=POINT'
        do j=1, jmax
            do i=1, imax
                if (problem=="reservoir") then 
                    if (((IAND(FLAG(i,j),  C_F)/=0).AND.(FLAG(i,  j)<C_E)).OR. &
                    ((IAND(FLAG(i+1,j),C_F)/=0).AND.(FLAG(i+1,j)<C_E))) then
                        write(4,*) x(i), z(j), PSI(i,j)
                    else
                        write(4,*) x(i), z(j), 0.0
                    end if
                else
                    write(4,*) x(i)+disp, z(j), PSI(i,j)
                end if
            end do
        end do
    else
        continue
    end if
    
    return
    end subroutine WR_PSI
    !

    !
    !%%%%%%%%%%%%%%%%%%%%%%%% WR_FREESURFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine WR_FREESURFACE(imax, dx, H, t, hi)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !   writes free surface profiles
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    USE nrtype
    implicit none
    integer, intent(IN):: imax
    real(RP), intent(IN):: t, hi
    real(RP), dimension(0:), intent(IN):: dx, H
    !-----------------------------------------------------------------------
    !
    integer :: i
    real(RP):: x
    open(5,FILE='fsurface.plt')
    write (5,'(4hZONE 1x2hT= 1x3h"T= f8.3,1x1h" )') t
    x = dx(1)/2.0
    do i=1, imax
        x  = x + (dx(i)+dx(i-1))/2.0
        !			write (8,50) x, (H(i)-1.0)*hi
        write (5,50) x, H(i)
50      format (f10.5,5x,f12.8)

    end do
    end subroutine WR_FREESURFACE
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine WR_RES2(FLAG, imax, jmax, H, dz, t, h0)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !   writes residual water as a function of time
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    USE nrtype
    USE defs
    implicit none
    integer, intent(IN):: imax, jmax
    integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
    real(RP), dimension(0:), intent(IN):: H, dz
    real(RP), intent(IN):: t, h0
    !
    !-----------------------------------------------------------------------
    !
    integer:: j, js, i
    real(RP):: too, hh
    open(6,FILE='res.plt')
    if (t==0.0) THEN
        write (6,'(4hZONE 1x2hT= 1x3h"h= f8.3,1x1h" )') h0
    end if
    !	  write (4,40) t, (H(1)-h0)*hi
    !write (6,*) t, H(1), H(imax)-h0
    !write (6,*) t*sqrt(9.81/h0), H(1)/h0
    write (6,*) t, H(1)-h0
    !40 format (f7.3,5x,f15.8)

    return
    end subroutine WR_RES2

    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine WR_RES(FLAG, jmax, H, dz, C, t, h0)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !   writes residual water as a function of time
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    USE nrtype
    USE defs
    implicit none
    integer, intent(IN):: jmax
    integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
    real(RP), dimension(0:), intent(IN):: H, dz
    real(DP), dimension(0:,0:), intent(IN) :: C
    real(RP), intent(IN):: t, h0
    !
    !-----------------------------------------------------------------------
    !
    integer:: j, js, i
    real(RP):: too, hh1, hh2
    open(6,FILE='res.plt')
    if(VOLUM_TRACKING) THEN
        if (t==0.0) THEN
            write (6,'(4hZONE 1x2hT= 1x3h"h= f8.3,1x1h" )') h0
        end if
        i = 1
        too = t*sqrt(9.81/h0)
        hh1 = 0.0
        hh2 = 0.0
        do j=1, jmax
            hh1 = hh1 + C(1,j)*dz(j)
            hh2 = hh2 + C(100,j)*dz(j)
        end do
        write (6,*) t, hh1-h0, hh2-h0
40      format (2f12.5)
    else
        write (6,"(2d20.3)") too, H(1)/h0
    end if
    !write (4,*) t, (H(33)-h0)/h0
    !write (4,*) t, (H(10)-h0), (H(60)-h0)
    !write (4,*) t, (H(1)-h0)/0.65

    return
    end subroutine WR_RES
    !
    !%%%%%%%%%%%%%%%%%%%%%%%% WR_FREESURFACE3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine WR_FREESURFACE3(problem, FLAG, C, imax, jmax, dx, dz, t, h0, vxb, vyb, nby, Ls)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !   writes free surface profiles
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    USE nrtype
    use defs
    implicit none
    CHARACTER (LEN=30) :: problem
    integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
    integer, intent(IN):: imax, jmax, nby
    real(RP), intent(IN):: t, h0, Ls
    real(RP), dimension(0:), intent(IN):: dx, dz
    real(DP), dimension(0:,0:), intent(IN) :: C
    real(DP), dimension(:), intent(IN) :: vxb, vyb
    !-----------------------------------------------------------------------
    !
    integer :: i, j, water, ib, k
    real(RP):: x, hh, zb, diff, mindiff
    logical found
    real(RP), dimension(0:imax+1) :: hs
    open(7,FILE='fsurface2.plt',STATUS='UNKNOWN')
    write (7,'(4hZONE 1x2hT= 1x3h"T= f8.3,1x1h" )') t

    x = dx(1)/2.0
    do i=1, imax
        x  = x + (dx(i)+dx(i-1))/2.0
        !find nearest immersed bundary point to assign the y coordinate of the immersed boundary point
        hh=0.0
        do j=1, jmax
            !if ( .NOT.(IAND(FLAG(i,j),  C_E) /= 0) .OR. (FLAG(i,j) < C_O) ) then
            !    do k=1, j-1
            !        hh=hh+dz(j)
            !    end do
                hh=hh+c(i,j)*dz(j)
            !end if
       end do
        
        write (7,50) x, hh
50      format (f10.5,5x,f13.8)
     
        
    end do

    end subroutine WR_FREESURFACE3

    subroutine PLT(U, W, P, C, FLAG, imax, jmax, dx, dz, tt)

    USE nrtype
    use defs
    
    implicit none
    integer, intent(IN) :: imax, jmax
    real(RP), intent(IN) :: tt
    integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
    real(RP), dimension(0:,0:), intent(IN) :: U, W, P
    real(DP), dimension(0:,0:), intent(IN) :: C
    real(RP), dimension(0:), intent(IN):: dx, dz

    integer :: i, j
    real (DP):: x, z
    real(RP), dimension(0:1000,0:1000) :: T

    open (8,FILE='vof.plt')

    !write(18,*) ' VARIABLES = "X","Y","U","V","P","T","F","FLAG" '
    !write(18,*) '    ZONE   I=',imax+1,'   J=',jmax+1,'    '
    write(8,*) ' VARIABLES = "X","Y","F"'
    write (8,*)	'ZONE T="',tt,'"',', ', 'I=',imax,',J=',jmax,',   F=POINT'
    write(8,*)

    x=-dx(0)/2.d0
    z=-dz(0)/2.d0

    do j=1, jmax
        z = z + (dz(j)+dz(j-1))/2.d0
        x=-dx(0)/2.d0
        do i=1, imax
            x=x + (dx(i)+dx(i-1))/2.d0
            !f (FLAG(i,j) /= C_B) THEN
                write(8,40) x, z, c(i,j)
40              format (3f12.5)
            !else
            !    write(8,40) x, z, 5.d0
                !	write(18,*) xX,yY,.5*(u(i,j)+u(i+1,j)),.5*(v(i,j)+v(i,j+1)),p(i,j),T(i,j),30.0,flag(i,j)

            !end if

        end do


    end do

    !CLOSE(UNIT=18)

    return
    end
    
!
!%%%%%%%%%%%%%%%%%%%%%%%% WR_UWP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
	subroutine WR_UWP(problem, imax, jmax, x, z, dx, dz, U, W, P, t)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Writes u, w, p and div into uwp file  for visualization
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
	use nrtype
    use defs
	implicit none
    character (LEN=30), INTENT(IN) :: problem
	integer, intent(IN):: imax, jmax
	real(RP), intent(IN):: t
	real(RP), dimension(0:),intent(IN):: x, z, dx, dz
	real(RP), dimension(0:,0:),intent(IN):: U, W, P
!
!-----------------------------------------------------------------------
!
	integer i, j, ip
	real(RP):: nu, pt, pb, pp, uu, ww, divt, divb, divv, F, xr, xl, xx, zz, ui, uim1, uni, &
                unim1, uun, vcmax, vor, dzn, dzs, dxe,dxw, &
                dwdxn, dwdxs, dudze, dudzw
    real(DP):: KC, D, pi, A, period, disp, vc
    real(RP), dimension(0:imax):: xp
    disp=0.0d0
    if (problem=="circle") then 
        KC=5.0d0
        D=1.0d0
        pi=4.0d0*datan(1.0d0)
        A=KC*D/2.0d0/pi
        period=1.0d0                                    !period of the oscillating cylinder
        disp=-A*dsin(2.0d0*pi*t/period)                 !displacement of the cylinder
        vc=-A*2.0d0*pi/period*dcos(2.0d0*pi*t/period)   !velocity of the cylinder
        vcmax=A*2.0d0*pi/period
        if (movbound) then
            pt=21.2
            !-x coordinates of the inertial system
            do i=1, imax
                xp(i)=x(i)+disp
            end do
            ip=0
            do i=1, imax
                if ((pt.gt.xp(i-1)).AND.(pt.lt.xp(i))) then
                    xr=xp(i)-pt
                    xl=pt-xp(i-1)
                    !velocity of the fluid element in the inertial system
                        ui=U(i,jmax/2)+vc        
                        uim1=U(i-1,jmax/2)+vc
                    !fluid velocities in non-inertial system
                        uni=U(i,jmax/2)
                        unim1=U(i-1,jmax/2)
                    !velocity at pt in the inertial system
                    uu=(ui*xl+uim1*xr)/dx(i)
                    !velocity at pt in the non-inertial system
                    uun=(uni*xl+unim1*xr)/dx(i)
                    exit 
                end if
            end do
            open(40,FILE='uinertial.plt')
            if (t==0.0) write (40,*) 'VARIABLES= ','"t",','"ui",','"uni",','"u0i",','"u0ni",','"uimaxi",','"uimaxni",','"vc"'
            write(40,*) t, uu/vcmax, uun/vcmax, (U(0,jmax/2)+vc)/vcmax, U(0,jmax/2)/vcmax,  &
                        (U(imax,jmax/2)+vc)/vcmax, U(imax,jmax/2)/vcmax, vc/vcmax
        end if
    end if
    if (t.ge.200.0) then 
        open(9,FILE='uwp.plt')
        write (9,*) 'VARIABLES= ','"X",','"Z",','"U",','"W",','"P",','"Un",','"vor",','"i",','"j"'
        write (9,*)	'ZONE T="',T,'"',', ', 'I=',imax,',J=',jmax,',   F=POINT'
        zz = -dz(0)/2.0
        do j=1, jmax
            zz = zz + (dz(j)+dz(j-1))/2.0
            xx = -dx(0)/2.0
            do i=1, imax
                xx = xx + (dx(i)+dx(i-1))/2.0
                !compute vorticity
                dxe=(dx(i)+dx(i+1))/2.0
                dxw=(dx(i)+dx(i-1))/2.0
                dzn=(dz(j)+dz(j+1))/2.0
                dzs=(dz(j)+dz(j-1))/2.0
                dwdxn=0.5*((W(i+1,j)-W(i,j))/dxe+(W(i,j)-W(i-1,j))/dxw)
                dwdxs=0.5*((W(i+1,j-1)-W(i,j-1))/dxe+(W(i,j-1)-W(i-1,j-1))/dxw)
                dudze=0.5*((U(i,j+1)-U(i,j))/dzn+(U(i,j)-U(i,j-1))/dzs)
                dudzw=0.5*((U(i-1,j+1)-U(i-1,j))/dzn+(U(i-1,j)-U(i-1,j-1))/dzs)
                vor=0.5*(dwdxn+dwdxs)-0.5*(dudze+dudzw)
                write(9,*) xx+disp, zz, (U(i,j)+U(i-1,j))/2.0+vc, (W(i,j)+W(i,j-1))/2.0, &
                            P(i,j), (U(i,j)+U(i-1,j))/2.0, vor, i, j
            end do
        end do
    end if
    
	return
    end subroutine WR_UWP
    
     !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%% WR_DAM_BASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine WR_DAM_BASE(imax, jmax, t, ro, GZ, P, h0, dz, H, alpha)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  Writes hydrodynamic pressures on dam face
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    
    implicit none
    integer, intent(IN):: imax, jmax
    real(RP), intent(IN):: t, ro, GZ, h0, alpha
    real(RP), dimension(0:,0:), intent(IN):: P
    real(RP), dimension(0:), intent(IN):: dz, H
    !
    !-----------------------------------------------------------------------
    !
    integer:: j, js
    real(RP):: F, M, ps, F0, M0, z, Ff, dd, P0
    real(RP), DIMENSION(0:jmax):: pd
    open(10,FILE='dam-base.plt',STATUS='UNKNOWN')
    if (t==0) write (10,*)	'ZONE T= ','PRESENT'
    P0=ro*(-GZ)*h0
    F0 = ro*(-GZ)*h0**2/2.0
    z = dz(1)/2.0
    do j=1, jmax
        if (z<h0) then
            ps = (h0-z)*ro*(-GZ)
        else
            ps = 0.0
        end if
        pd(j) = P(1,j)-ps
        z = z + (dz(j)+dz(j+1))/2.0
    end do
    do j=1, jmax-1
        dd = (dz(j)+dz(j+1))/2.0
        F = F + (pd(j)+pd(j+1))*dd/2.0
    end do
    F = F / F0
    ps = (h0-dz(1)/2.0)*ro*(-GZ)
    !write (10,*) t/0.4, pd(1)/(ro*(-GZ)*h0), F
    !write (10,*) t/0.509,  F
    write (10,*) t, pd(1)/alpha/P0
    return
    end subroutine WR_DAM_BASE
    
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%% DAM_FACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   subroutine DAM_FACE(problem, FLAG, imax, jmax, P, H, dz, t, h0, ro, GZ, alpha)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Writes hydrodynamic pressures on dam face
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      use nrtype
      use defs
	  implicit none
      character (LEN=30), INTENT(IN) :: problem
      integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
	  integer, intent(IN):: imax, jmax
	  real(RP), intent(IN):: t, h0, ro, GZ, alpha
	  real(RP), dimension(0:), intent(IN):: dz, H
	  real(RP), dimension(0:,0:), intent(IN):: P
      
!
!-----------------------------------------------------------------------
!
	  integer:: i, j
	  real(RP) :: z, ps, p0
	  open(12,FILE='dam-face.plt',STATUS='UNKNOWN')
      write (12,*) 'ZONE T=', '"',t,'"'
	  !write (12,'(4hZONE 1x2hT= 1x3h"T= f8.6,1x1h" )') t
	  z = -dz(0)/2.0
      p0 = ro*(-GZ)*h0
      do j=1, jmax
          z = z + (dz(j)+dz(j-1))/2.0
          ps = (h0-z)*ro*(-GZ)
          write (12,*) (P(1,j)-ps)/p0/alpha, z/h0
      end do
    
      return
   end subroutine DAM_FACE
   
   !
    !%%%%%%%%%%%%%%%%%%%%%%%% WR_SOLITARY_SHELF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine WR_SOLITARY_SHELF(imax, dx, H, t, h0)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !   writes free surface profiles
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    USE nrtype
    implicit none
    integer, intent(IN):: imax
    real(RP), intent(IN):: t, h0
    real(RP), dimension(0:), intent(IN):: dx, H
    !-----------------------------------------------------------------------
    !
    integer :: i, ia, ib, ic, id
    real(RP):: ha, hb, hc, hd, erra, errb, errc, errd, err, x
    open(44,FILE='sol-shelf.plt')
    ha=1.0e5
    hb=1.0e5
    hc=1.0e5
    hd=1.0e5
    erra=1.0e5
    errb=1.0e5
    errc=1.0e5
    errd=1.0e5
    x=-dx(0)/2.0
    do i=1, imax
        x=x+(dx(i)+dx(i-1))/2.0
        err=abs(x-1.5936)
        if (err.lt.erra) then 
            ha=H(i)
            erra=err
            ia=i
        end if
        err=abs(x-2.762)
        if (err.lt.errb) then 
            hb=H(i)
            errb=err
            ib=i
        end if
        err=abs(x-3.651)
        if (err.lt.errc) then 
            hc=H(i)
            errc=err
            ic=i
        end if
        err=abs(x-4.5146)
        if (err.lt.errd) then 
            hd=H(i)
            errd=err
            id=i
        end if
    end do
    write (44,*) t, (ha-h0)/h0, (hb-h0)/h0, (hc-h0)/h0, (hd-h0)/h0
    
    return
    end subroutine WR_SOLITARY_SHELF
   
