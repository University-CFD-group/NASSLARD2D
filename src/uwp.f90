    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMP_delt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine COMP_delt (problem, delt, imax, jmax, dx, dz, U, W, Re, St, cfl, H, c, uxb, uyb, nm)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    ! Computation of adaptive time stepsize satisfying 
    ! the CFL stability criteria and set the flag "write" if some data
    ! has to be written into a file. 
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use defs
    implicit none
    character (LEN=30) :: problem
    integer, INTENT(IN) :: imax, jmax, nm
    real(RP), INTENT(IN) :: Re, St, cfl, c, uxb, uyb
    real(RP), INTENT(OUT) :: delt
    real(RP), DIMENSION(0:,0:), INTENT(IN) :: U, W
    real(RP), DIMENSION(0:), INTENT(IN) :: dx, dz, H
    !
    !     local variables
    !
    integer  :: i, j
    real(RP) :: uu, ww, deltu, deltw, dxmin, dzmin, deltRePr, &
                deltc, hmax, temp, deltub, deltwb, ubmax, wbmax, deltvof
    !
    ! delt satisfying CFL conditions
    !
    if ( cfl >= 1.0e-10 ) then ! else no time stepsize control
        deltu = 1.0E10
        deltw = 1.0E10
        deltub = 1.0E10
        deltwb = 1.0E10
        dxmin = 1.0E10
        dzmin = 1.0E10
        deltc = 1.0E10
        deltvof=1.0E10

        if (problem /= "dcavity" .AND. problem /= "backstep" .AND. problem /="circle") then
            do i=0, imax
                do j=0, jmax
                    temp = 2.0*dx(i)*dz(j)/(c+1.0)/(dx(i)+dz(j))/St
                    if (temp<deltc)  deltc = temp
                end do
            end do
        end if

        do i = 0, imax
            if (dx(i)<dxmin) dxmin = dx(i)
            do j = 1, jmax
                uu = abs(U(i,j))
                if (uu /= 0.0) then
                    if ((dx(i)+dx(i+1))/2.0/uu < deltu) then 
                        deltu = (dx(i)+dx(i+1))/2.0/uu/St
                    end if
                    
                end if
            end do
        end do
       
        if (uxb/=0.0) deltub=dxmin/uxb/St
        if (uyb/=0.0) deltwb=dzmin/uyb/St
          
        do i = 1,imax
            do j = 0, jmax
                if(dz(j)<dzmin)  dzmin = dz(j)
                ww = ABS(W(i,j))
                if (ww /= 0.0) then
                    if ((dz(j)+dz(j+1))/2.0/ww < deltw) then 
                        deltw = (dz(j)+dz(j+1))/2.0/ww/St
                    end if
                end if
            end do
        end do
        
        if (VOLUM_TRACKING) then
            do i=1, imax
                do j=1, jmax
                    if (dx(i)/abs(U(i,j)-U(i-1,j))<deltvof) deltvof = dx(i)/abs(U(i,j)-U(i-1,j))
                    if (dz(j)/abs(W(i,j)-W(i,j-1))<deltvof) deltvof = dz(j)/abs(W(i,j)-W(i,j-1))
                end do
            end do
        end if
        deltRePr = Re/2.0*(dxmin**2*dzmin**2)/(dxmin**2+dzmin**2)
        if (deltu < deltw) then
            if (deltu < deltRePr) then
                delt = deltu
            else
                delt = deltRePr
            end if
        else
            if (deltw < deltRePr) then
                delt = deltw
            else
                delt = deltRePr
            end if
        end if
        if (deltc<delt) delt = deltc
        if (deltu<delt) delt=deltu
        if (deltw<delt) delt=deltw
        if (deltvof<delt) delt = deltvof
        if (deltub<delt) delt=deltub
        if (deltwb<delt) delt=deltwb
        delt = cfl * delt   ! multiply by safety factor
    end if
    return
    end subroutine COMP_delt
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMP_FG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine COMP_FG (problem, upwind, U, W, UNM1, WNM1, F, G, imax, jmax, gamma, &
                        delt, St, Re, Fr, dx, dz, GZ, axnm1, axn, P, ro, pflag, wE, wW, &
                        induin, indvin, indu, indv, gFIuin, gFIvin, gFIu, gFIv, invaguin, invagvin, invagu, &
                        invagv, ghuin, ghvin, ghu, ghv, ivp, jvp, ppu, ppv, IMu, IMv, uxb, uyb)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  Compute tentative velocity field for horizontal and vertical momentum
    !  equations.
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use interfaces2
    use defs
    implicit none
    character (LEN=30), INTENT(IN) :: problem, upwind
    integer, INTENT(IN) :: imax, jmax, wE, wW, ivp, jvp, ppu, ppv
    real(RP), INTENT(IN) :: gamma, delt, St, Re, Fr, GZ, axnm1, axn, ro, uxb, uyb
    real(RP), DIMENSION(0:), INTENT(IN) :: dx, dz
    real(RP), DIMENSION(0:,0:), INTENT(IN) :: U, W, UNM1, WNM1, P
    real(RP), DIMENSION(0:,0:), INTENT(OUT) :: F, G
    integer(I2B), dimension(0:,0:), intent(IN) :: pflag
    integer, dimension(:,:), intent(IN) :: induin, indvin, indu, indv
    integer, dimension(:,:,:), intent(IN) :: gFIuin, gFIvin, gFIu, gFIv
    real(DP), dimension(:,:,:), intent(IN) :: invaguin, invagvin, invagu, invagv
    real(DP), dimension(:,:), intent(IN) :: ghuin, ghvin, ghu, ghv
    integer, dimension(:), intent(IN) :: IMu, IMv
    !
    !-----------------------------------------------------------------------
    !
    integer :: i, j, ii, jj, i1, j1, i2, j2
    real(RP):: dxu, dzu, dxw, dzw, dzn, dzs, dxe, dxww, conn, connm1, hn, hnm1, press, Q3U, Q3W, xx, zz
    real(DP), dimension(3):: fi
    real(DP) a0, a1, a2, fiu, fiv

    do i = 1, imax-1
        dxu = (dx(i)+dx(i+1))/2.0
        do j = 1, jmax
            dzu = dz(j)
            dzn = (dz(j)+dz(j+1))/2.0
            dzs = (dz(j)+dz(j-1))/2.0
            !
            !  only if both adjacent cells are fluid cells
            !
            if (((IAND(pflag(i,j),  C_F)/=0).AND.(pflag(i,j)  <C_E)).AND.&
                ((IAND(pflag(i+1,j),C_F)/=0).AND.(pflag(i+1,j)<C_E))) then
                if (upwind=="quick") then
                    conn = CONV_UQ(i, j, imax, jmax, U, W, dx, dz)		        !Convective term (QUICK) at n time level
                    connm1 = CONV_UQ(i, j, imax, jmax, UNM1, WNM1, dx, dz)		!Convective term (QUICK) at n-1 time level
                elseif (upwind=="fou") then
                    conn = CONV_U(i, j, U, W, dx, dz, gamma)				    !Convective term (FOU) at n time level
                    connm1 = CONV_U(i, j, UNM1, WNM1, dx, dz, gamma)            !Convective term (FOU) at n-1 time level
                else
                    write (6,*) ' invalid upwind procedure'
                    stop
                end if
                Q3U = 1.0/Fr**2*(P(i,j)-P(i+1,j))*dz(j)/ro
                hn = 1.0/Re*DIFF_U(i, j, U, dx, dz) - conn 
                hnm1 = 1.0/Re*DIFF_U(i, j, UNM1, dx, dz) - connm1 
                F(i,j) = U(i,j) + St*delt/dxu/dzu*(3.0/2.0*hn-1.0/2.0*hnm1+Q3U) - delt*axn
                !F(i,j) = U(i,j) + St*delt/dxu/dzu*hn - delt*axn
            else
                F(i,j)=U(i,j)
            end if
        end do
    end do
    !
    !  compute flux field G
    !
    do i = 1, imax
        dxw = dx(i)
        dxe = (dx(i)+dx(i+1))/2.0
        dxww = (dx(i)+dx(i-1))/2.0
        do j = 1, jmax-1
            dzw = (dz(j)+dz(j+1))/2.0
            !
            !  only if both adjacent cells are fluid cells
            !
            if (((IAND(pflag(i,j),  C_F)/=0).AND.(pflag(i,j)<C_E)).AND.&
                ((IAND(pflag(i,j+1),C_F)/=0).AND.(pflag(i,j+1)<C_E)))then
                Q3W = 1.0/Fr**2*(P(i,j)-P(i,j+1))*dx(i)/ro
                if (upwind=="quick") then
                    conn = CONV_WQ(i, j, imax, jmax, U, W, dx, dz)			    !Convective term (QUICK) at n time level
                    connm1 = CONV_WQ(i, j, imax, jmax, UNM1, WNM1, dx, dz)		!Convective term (QUICK) at n-1 time level
                elseif (upwind=="fou") then
                    conn = CONV_W(i, j, U, W, dx, dz, gamma)                    !Convective term (FOU) at n time level
                    connm1 = CONV_W(i, j, UNM1, WNM1, dx, dz, gamma)			!Convective term (FOU) at n-1 time level
                else
                    write (6,*) ' invalid upwind procedure'
                    stop
                end if
                hn = 1.0/Re*DIFF_W(i, j, W, dx, dz) - conn + 1.0/Fr**2*SRC_W(i, j, GZ, dx, dz)
                hnm1 = 1.0/Re*DIFF_W(i, j, WNM1, dx, dz) - connm1 + 1.0/Fr**2*SRC_W(i, j, GZ, dx, dz)
                G(i,j) = W(i,j) + St*delt/dxw/dzw*(3.0/2.0*hn-1.0/2.0*hnm1+Q3W) 
                !G(i,j) = W(i,j) + St*delt/dxw/dzw*hn
            else
                G(i,j)=W(i,j)
            end if
        end do
    end do
    !interpolate F and G values at inactive cells form the neighboring active cells
    !at ianctive cells which lie in the fluid
    if (immersed) then 
        do i=1, ivp
            ii=induin(i,1)
            jj=induin(i,2)
            i1=gFIuin(i,1,1)          !i index of first interpolation point
            j1=gFIuin(i,1,2)          !j index of first interpolation point
            i2=gFIuin(i,2,1)          !i index of second interpolation point
            j2=gFIuin(i,2,2)          !j index of second interpolation point
            fi(1)=uxb                 !horizontal solid velocity 
            fi(2)=F(i1,j1)            !tentative velocity at first interpolation point
            fi(3)=F(i2,j2)            !tentative velocity at second interpolation point
            !compute linear interpolation constants a0, a1 and a2
            a0=fi(1)+fi(2)*invaguin(i,1,2)+fi(3)*invaguin(i,1,3)
            a1=fi(1)+fi(2)*invaguin(i,2,2)+fi(3)*invaguin(i,2,3)
            a2=fi(1)+fi(2)*invaguin(i,3,2)+fi(3)*invaguin(i,3,3)
            !compute u variable at ghost point using linear interpolation
            F(ii,jj)=a0+a1*ghuin(i,1)+a2*ghuin(i,2)
        end do
        do i=1, jvp
            ii=indvin(i,1)
            jj=indvin(i,2)
            i1=gFIvin(i,1,1)          !i index of first interpolation point
            j1=gFIvin(i,1,2)          !j index of first interpolation point
            i2=gFIvin(i,2,1)          !i index of second interpolation point
            j2=gFIvin(i,2,2)          !j index of second interpolation point
            fi(1)=uyb                 !vertical solid velocity 
            fi(2)=G(i1,j1)            !tentative velocity at first interpolation point
            fi(3)=G(i2,j2)            !tentative velocity at second interpolation point
            !compute linear interpolation constants a0, a1 and a2
            a0=fi(1)+fi(2)*invagvin(i,1,2)+fi(3)*invagvin(i,1,3)
            a1=fi(1)+fi(2)*invagvin(i,2,2)+fi(3)*invagvin(i,2,3)
            a2=fi(1)+fi(2)*invagvin(i,3,2)+fi(3)*invagvin(i,3,3)
            !compute u variable at ghost point using linear interpolation
            G(ii,jj)=a0+a1*ghvin(i,1)+a2*ghvin(i,2)
        end do
        !at the ghost cells
        do i=1, ppu
            ii=indu(i,1)
            jj=indu(i,2)
            i1=gFIu(i,1,1)          !i index of first interpolation point
            j1=gFIu(i,1,2)          !j index of first interpolation point
            i2=gFIu(i,2,1)          !i index of second interpolation point
            j2=gFIu(i,2,2)          !j index of second interpolation point
            fi(1)=0.0d0             !solid velocity (for the case of stationary boundary)
            fi(2)=F(i1,j1)          !velocity at first interpolation point
            fi(3)=F(i2,j2)          !velocity at second interpolation point
            !compute linear interpolation constants a0, a1 and a2
            a0=fi(1)+fi(2)*invagu(i,1,2)+fi(3)*invagu(i,1,3)
            a1=fi(1)+fi(2)*invagu(i,2,2)+fi(3)*invagu(i,2,3)
            a2=fi(1)+fi(2)*invagu(i,3,2)+fi(3)*invagu(i,3,3)
            !compute u variable at ghost point using linear interpolation
            fiu=a0+a1*ghu(i,1)+a2*ghu(i,2)
            !if imagine of ghost point is used then re-compute the ghost point value
            if (IMu(i)==1) fiu=2.0d0*fi(1)-fiu
            F(ii,jj)=fiu
        end do
        do i=1, ppv
            ii=indv(i,1)
            jj=indv(i,2)
            i1=gFIv(i,1,1)          !i index of first interpolation point
            j1=gFIv(i,1,2)          !j index of first interpolation point
            i2=gFIv(i,2,1)          !i index of second interpolation point
            j2=gFIv(i,2,2)          !j index of second interpolation point
            fi(1)=0.0d0             !solid velocity (for the case of stationary boundary)
            fi(2)=G(i1,j1)          !velocity at first interpolation point
            fi(3)=G(i2,j2)          !velocity at second interpolation point
            !compute linear interpolation constants a0, a1 and a2
            a0=fi(1)+fi(2)*invagv(i,1,2)+fi(3)*invagv(i,1,3)
            a1=fi(1)+fi(2)*invagv(i,2,2)+fi(3)*invagv(i,2,3)
            a2=fi(1)+fi(2)*invagv(i,3,2)+fi(3)*invagv(i,3,3)
            !compute u variable at ghost point using linear interpolation
            fiv=a0+a1*ghv(i,1)+a2*ghv(i,2)
            !if imagine of ghost point is used then re-compute the ghost point value
            if (IMv(i)==1) fiv=2.0d0*fi(1)-fiv
            G(ii,jj)=fiv
        end do
    end if
    !
    ! F and G at external boundary
    !
    do j = 1, jmax
        if (wW == 3 .OR. wW == 4) then 
            F(0,j) = U(0,j) - delt*axn
        else
            F(0,j) = U(0,j)
        end if
        if (wE == 3 .OR. wE == 4) then
			F(imax,j) = U(imax,j) - delt*axn
		else
			F(imax,j) = U(imax,j)
		end if
    end do
    do i = 1, imax
        G(i,0) = W(i,0)
        G(i,jmax) = W(i,jmax)
    end do
    
    return
    end subroutine COMP_FG
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADAP_UWP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine ADAP_UW (U, W, F, G, P, pflag, H, imax, jmax, dx, dz, z, Fr, St, ro, delt, c)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  Compute the U, W fields if adjacent cells are fluid cells
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use defs
    implicit none
    integer, INTENT(IN) :: imax, jmax
    real(RP), INTENT(IN) :: delt, St, Fr, ro
    integer(I2B), dimension(0:,0:), intent(IN) :: pflag
    real(RP), DIMENSION(0:), INTENT(IN):: dx, dz, H, z
    real(RP), DIMENSION(0:,0:), INTENT(IN) :: F, G
    real(RP), DIMENSION(0:,0:), INTENT(INOUT) :: P
    real(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, W
    real(DP), dimension(0:,0:), intent(IN) :: c

    !
    !-----------------------------------------------------------------------
    !
    integer :: i, j, ifull, jfull, js
    real(RP):: d, P_N, P_E, P_S,  P_W, dxx, dzz, nu, Q1U, Q3U, Q1W, Q3W, UN, WN
    ifull = 0
    jfull = 0

    do i = 1, imax - 1
        dxx = (dx(i)+dx(i+1))/2.0
        do j = 1, jmax
            !
            ! only if both adjacent cells are fluid cells
            !
            if (((IAND(pflag(i,j),  C_F)/=0).AND.(pflag(i,j)  <C_E)).AND.&
                ((IAND(pflag(i+1,j),C_F)/=0).AND.(pflag(i+1,j)<C_E))) then
                P_E = P(i+1,j)
                Q3U = 1.0/Fr**2*(P(i,j)-P_E)*dz(j)/ro
                U(i,j) = F(i,j)  + St*delt*Q3U/dxx/dz(j)
            end if
        end do
    end do
    do i = 1, imax
        do j = 1, jmax - 1
            dzz = (dz(j)+dz(j+1))/2.0
            !
            ! only if both adjacent cells are fluid cells
            !            
            if (((IAND(pflag(i,j),  C_F)/=0).AND.(pflag(i,j)  <C_E)).AND.&
                ((IAND(pflag(i,j+1),C_F)/=0).AND.(pflag(i,j+1)<C_E))) then
                P_N = P(i,j+1)
                if (.NOT.(VOLUM_TRACKING).AND.(pflag(i,j+1)>=C_O)) then
                    d = H(i) - z(j-1) - dz(j)/2.0
                    nu = (dz(j)+dz(j+1))/2.0/d
                    P_N = nu*P(i,j+1)+(1.0-nu)*P(i,j)
                end if
                Q3W = 1.0/Fr**2*(P(i,j)-P_N)*dx(i)/ro
                W(i,j) = G(i,j) + St*delt*Q3W/dzz/dx(i)
            end if

        end do
    end do
    
    return
    end subroutine ADAP_UW
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMP_RHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine COMP_RHS (F, G, P, RHS, pflag, imax, jmax, delt, Fr, St, Fa, ro, dx, dz)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  Compute right-hand side of Pressure-Poisson equation.
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use defs
    implicit none
    integer, INTENT(IN) :: imax, jmax
    real(RP), INTENT(IN) :: delt, Fr, St, Fa, ro
    real(RP), DIMENSION(0:), INTENT(IN) :: dx, dz
    integer(I2B), DIMENSION(0:,0:), INTENT(IN) :: pflag
    real(RP), DIMENSION(0:,0:), INTENT(IN) :: F, G, P
    real(RP), DIMENSION(0:,0:), INTENT(OUT) :: RHS
    !
    !-----------------------------------------------------------------------
    !
    integer :: i,j
    real(RP):: k
    k = St*Fa*delt/Fr
    do i = 1, imax
        do j = 1, jmax
            !
            !           only for fluid and non-surface cells
            !
            if ((IAND(pflag(i,j),C_F)/=0).AND.(pflag(i,j)<C_O)) then
                RHS(i,j) = ro*Fr**2/St/delt*( ( F(i,j) - F(i-1,j) ) / dx(i) &
                           + ( G(i,j) - G(i,j-1) ) / dz(j)) 
            end if
            
        end do
    end do
    return
    end subroutine COMP_RHS

    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POISSON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine POISSON (problem, P, RHS, pflag, imax, jmax, ppp, dx, dz, z, eps, iter, &
                        itermax, omega, res, ifull, St, Fr, Fa, delt, H, c, ro, indp, gFIp, &
                        invagp, ghp, IMp, ap, bp, cp, q)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  Iteration for Poisson's equation
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use defs
    implicit none
    character (LEN=30), INTENT(IN) :: problem
    integer, INTENT(IN) :: imax, jmax, ppp, itermax, ifull
    integer, INTENT(INOUT) :: iter
    real(RP), INTENT(IN) ::  eps, omega, St, Fr, Fa, delt, ro
    real(RP), INTENT(OUT) :: res, q
    integer(I2B), DIMENSION(0:,0:), INTENT(IN) :: pflag
    real(RP), DIMENSION(0:,0:), INTENT(INOUT) :: P
    real(RP), DIMENSION(0:,0:), INTENT(IN) :: RHS
    real(RP), DIMENSION(0:), INTENT(IN):: dx, dz, z, H
    real(DP), dimension(0:,0:), intent(INOUT) :: c
    integer, dimension(:,:), intent(IN) :: indp
    integer, dimension(:,:,:), intent(IN) :: gFIp
    real(DP), dimension(:,:,:), intent(IN) :: invagp
    real(DP), dimension(:,:), intent(IN) :: ghp
    integer, dimension(:), intent(IN) :: IMp
    real(DP), dimension(:), intent(OUT) :: ap, bp, cp
    !
    !-----------------------------------------------------------------------
    !
    integer i, j, ii, jj, i1, j1,i2, j2, pbound
    real(RP) :: add, beta_2, a, b, cc, d, e, f, k, beta_mod, KK
    real(RP) :: p0, dd, nu, P_N, P_E, P_W, P_S, dxe, dxw, dzn, dzs, zz, xx, alpha, dpdn, pbn
    real(DP), dimension(3):: fi
    real(DP):: fip, a0, a1, a2, a3
    !real(RP), dimension(:,:), allocatable :: PN
    integer eps_E, eps_W, eps_N, eps_S
    k = St*Fa*delt/Fr
    !OPEN(20,FILE='pressure.plt',STATUS='UNKNOWN')
    pbound=2
    p0 = 0.0
    do i = 1, imax
        do j = 1, jmax
            if ( IAND(pflag(i,j),C_F) /= 0 ) then
                p0 = p0 + P(i,j)*P(i,j)
            end if
        end do
    end do
    p0 = sqrt(p0/ifull)
    if (p0 < 0.0001) then
        p0 = 1.0
    end if
    if (absorption) then 
        !allocate (PN(0:imax+1,0:jmax+1))
        alpha=0.95   !absorption coeffcicient
        q=(1.0-alpha)/(1.0+alpha)/St/Fa/Fr
        dzs=(dz(0)+dz(1))/2.0
        KK=Fr**2*q*dzs/delt
    end if
    
    !
    !  SOR iteration
    !
    do iter = 1, itermax
        if (pbound==1) then 
            !
            ! relaxation for fluid cells 
            !
            do i = 1, imax         
                dxe = (dx(i)+dx(i+1))/2.0
                dxw = (dx(i)+dx(i-1))/2.0
                a = 1.0/dxe/dx(i)
                cc = 1.0/dxw/dx(i)
                b = a+cc
                do j = 1, jmax
                    P_N = P(i,j+1)
                    P_S=P(i,j-1)
                    P_E=P(i+1,j)
                    P_W=P(i-1,j)
                    dzn = (dz(j)+dz(j+1))/2.0
                    dzs = (dz(j)+dz(j-1))/2.0
                    d = 1.0/dzn/dz(j)
                    f = 1.0/dzs/dz(j)
                    e = d+f
                    beta_2 = -omega/(b+e)
                    !
                    ! five point star for interior fluid cells
                    !
                    if (pflag(i,j) == C_A) then
                        P(i,j) = (1.0-omega)*P(i,j) -beta_2*(a*P_E+cc*P_W+d*P_N+f*P_S-RHS(i,j))
                    else if ((IAND(pflag(i,j),C_F)/=0).AND.(pflag(i,j)<C_O)) then
                        include 'eps.h'
                        a = eps_E/dxe/dx(i)
                        cc = eps_W/dxw/dx(i)
                        b = a+cc
                        d = eps_N/dzn/dz(j)
                        f = eps_S/dzs/dz(j)
                        e = d+f
                        beta_mod = -omega/(b+e)
                        P(i,j) = (1.0-omega)*P(i,j) - beta_mod*(a*P_E+cc*P_W+d*P_N+f*P_S-RHS(i,j))
                    end if
                end do
            end do

            res = 0.0
            do i = 1,imax
                dxe = (dx(i)+dx(i+1))/2.0
                dxw = (dx(i)+dx(i-1))/2.0
                a = eps_E*1.0/dxe/dx(i)
                cc = eps_W*1.0/dxw/dx(i)
                b = a+cc
                do j = 1,jmax
                    dzn = (dz(j)+dz(j+1))/2.0
                    dzs = (dz(j)+dz(j-1))/2.0
                    d = eps_N*1.0/dzn/dz(j)
                    f = eps_S*1.0/dzs/dz(j)
                    e = d+f
                    !
                    ! only fluid cells
                    !
                    if ((IAND(pflag(i,j),C_F)/=0).AND.(pflag(i,j)< C_O)) then
                        P_N = P(i,j+1)
                        P_S=P(i,j-1)
                        P_E=P(i+1,j)
                        P_W=P(i-1,j)
                        add = a*(P_E-P(i,j))+cc*(P_W-P(i,j))+d*(P_N-P(i,j))+f*(P_S-P(i,j))-RHS(i,j)
                        res = res + add*add
                    end if
                end do
            end do
            res = sqrt(res/ifull)/p0
            if (res < eps) exit
        else
            !
            ! copy values at external boundary
            !
            do i = 1, imax
                if (absorption) then
                    P(i,0)=(1.0+K/2.0)*P(i,1)/(1.0-K/2.0)
                else
                    P(i,0) = P(i,1)
                end if
                P(i,jmax+1) = P(i,jmax)
            end do
            do j = 1,jmax
                P(0,j) = P(1,j)
                P(imax+1,j) = P(imax,j)
            end do

            if (immersed) then 
                do i=1, ppp
                    ii=indp(i,1)
                    jj=indp(i,2)
                    i1=gFIp(i,1,1)          !i index of first interpolation point
                    j1=gFIp(i,1,2)          !j index of first interpolation point
                    i2=gFIp(i,2,1)          !i index of second interpolation point
                    j2=gFIp(i,2,2)          !j index of second interpolation point
                    fi(1)=0.0d0             !dp/dn=0
                    fi(2)=P(i1,j1)          !pressure at first interpolation point
                    fi(3)=P(i2,j2)          !pressure at second interpolation point
                    !compute linear interpolation constants a0, a1 and a2
                    a0=fi(2)*invagp(i,1,2)+fi(3)*invagp(i,1,3)
                    a1=fi(1)+fi(2)*invagp(i,2,2)+fi(3)*invagp(i,2,3)
                    a2=fi(1)+fi(2)*invagp(i,3,2)+fi(3)*invagp(i,3,3)
                    !compute p variable at ghost point using linear interpolation
                    fip=a0+a1*ghp(i,1)+a2*ghp(i,2)
                    P(ii,jj)=fip
                    ap(i)=a0
                    bp(i)=a1
                    cp(i)=a2
                end do
            end if
            !
            !  relaxation method for fluid cells
            !
            do i = 1, imax
                dxe = (dx(i)+dx(i+1))/2.0
                dxw = (dx(i)+dx(i-1))/2.0
                a = 1.0/dxe/dx(i)
                cc = 1.0/dxw/dx(i)
                b = a+cc
                do j = 1, jmax
                    dzn = (dz(j)+dz(j+1))/2.0
                    dzs = (dz(j)+dz(j-1))/2.0
                    d = 1.0/dzn/dz(j)
                    f = 1.0/dzs/dz(j)
                    e = d+f
                    beta_2 = -omega/(b+e+1.0/k**2)
                    !beta_2 = -omega/(b+e)
                    !
                    ! only fluid cells
                    !
                    if ((IAND(pflag(i,j),C_F)/=0).AND.(pflag(i,j)<C_O)) then
                        P_N = P(i,j+1)
                        P_S=P(i,j-1)
                        P_E=P(i+1,j)
                        P_W=P(i-1,j)
                        if (.NOT.(VOLUM_TRACKING).AND.(pflag(i,j+1)>=C_O)) then
                            dd = H(i) - z(j-1) - dz(j)/2.0
                            nu = (dz(j)+dz(j+1))/2.0/dd
                            P_N = nu*P(i,j+1)+(1.0-nu)*P(i,j)
                        end if
                        P(i,j) = (1.0-omega)*P(i,j) -beta_2*(a*P_E+cc*P_W+d*P_N+f*P_S-RHS(i,j))
                    end if
                end do
            end do
            !
            ! computation of residual
            !
            res = 0.0
            do i = 1,imax
                dxe = (dx(i)+dx(i+1))/2.0
                dxw = (dx(i)+dx(i-1))/2.0
                a = 1.0/dxe/dx(i)
                cc = 1.0/dxw/dx(i)
                b = a+cc
                do j = 1,jmax
                    dzn = (dz(j)+dz(j+1))/2.0
                    dzs = (dz(j)+dz(j-1))/2.0
                    d = 1.0/dzn/dz(j)
                    f = 1.0/dzs/dz(j)
                    e = d+f
                    !
                    ! only fluid cells
                    !
                    if ((IAND(pflag(i,j),C_F)/=0).AND.(pflag(i,j)< C_O)) then
                        P_N = P(i,j+1)
                        P_S=P(i,j-1)
                        P_E=P(i+1,j)
                        P_W=P(i-1,j)
                        if (.NOT.(VOLUM_TRACKING).AND.(pflag(i,j+1)>=C_O)) then
                            dd = H(i) - z(j-1) - dz(j)/2.0
                            nu = (dz(j)+dz(j+1))/2.0/dd
                            P_N = nu*P(i,j+1)+(1.0-nu)*P(i,j)
                        end if
                        add = a*P_E-(b+e+1.0/k**2)*P(i,j)+cc*P_W+d*P_N+f*P_S-RHS(i,j)
                        !add = a*P_E-(b+e)*P(i,j)+cc*P_W+d*P_N+f*P_S-RHS(i,j)
                        res = res + add*add
                    end if
                end do
            end do
            res = sqrt(res/ifull)/p0
            if (res < eps) EXIT
        end if

    end do

    return
    end subroutine POISSON
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMP_RES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    LOGICAL FUNCTION COMP_RES (U, W, UT, WT, imax, jmax, ifull)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  Computes residual of Navier-Stokes equations
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    implicit none
    real(RP), DIMENSION(0:,0:), INTENT(IN) :: U, W, UT, WT
    integer, INTENT(IN) :: imax, jmax, ifull
    !
    !-----------------------------------------------------------------------
    !
    integer :: i, j
    real(RP):: eru, erw, erka, erep, umax, wmax, kamax, epmax, epsu, epsw, epska, epsep, eps
    eps = 1.0E-7
    eru = 0.0; erw = 0.0 ; erka = 0.0 ; erep = 0.0
    umax = 0.0 ; wmax = 0.0 ; kamax = 0.0 ; epmax = 0.0
    do i=1, imax
        do j=1, jmax
            epsu = U(i,j) - UT(i,j)
            eru = eru + abs(epsu)
            epsw = W(i,j) - WT(i,j)
            erw = erw + abs(epsw)
            if(abs(U(i,j))>umax) umax = abs(U(i,j))
            if(abs(W(i,j))>wmax) wmax = abs(W(i,j))
        end do
    end do
    eru = eru/ifull/max(umax,wmax)
    erw = erw/ifull/max(umax,wmax)

    if (eru<eps .AND. erw<eps) then  !steady state
        COMP_RES = .TRUE.
    else
        COMP_RES = .FALSE.
    end if
    print *, 'eru= ', eru, 'erw= ', erw

    return
    end FUNCTION COMP_RES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEIGHBORINGPRESSURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    subroutine NEIGHBORINGPRESSURE (ii, jj, P, FLAG, dx, dz, c, P_E, P_W, P_N, P_S)
    use nrtype
    use defs
    implicit none
    integer, INTENT(IN) :: ii, jj
    integer(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
    real(RP), DIMENSION(0:,0:), INTENT(INOUT) :: P
    real(DP), dimension(0:,0:), intent(IN) :: c
    real(RP), DIMENSION(0:), INTENT(IN):: dx, dz
    real(RP), INTENT(INOUT):: P_E, P_W, P_N, P_S
    
    integer:: i, j, hyd
    real(RP):: yip, yim, yi, xjp, xjm, xj, sl1, sl2, h, d, pf, eta
    
    if (FLAG(ii,jj+1)>=C_O) then    !north cell is a surface cell
        i=ii
        j=jj+1
        yip = c(i+1,j-1)*dz(j-1)+c(i+1,j)*dz(j)+c(i+1,j+1)*dz(j+1)
        yim = c(i-1,j-1)*dz(j-1)+c(i-1,j)*dz(j)+c(i-1,j+1)*dz(j+1)
        xjp = c(i-1,j+1)*dx(i-1)+c(i,j+1)*dx(i)+c(i+1,j+1)*dx(i+1)
        xjm = c(i-1,j-1)*dx(i-1)+c(i,j-1)*dx(i)+c(i+1,j-1)*dx(i+1)
        !DY/DX
        sl1 = 2.0*(yip-yim)/(dx(i+1)+2.0*dx(i)+dx(i-1))
        !
        !DX/DY
        sl2 = 2.0*(xjp-xjm)/(dz(j+1)+2.0*dz(j)+dz(j-1))
        if (abs(sl1)<=abs(sl2)) then                !horizontal
            h = (dz(j)+dz(j-1))/2.0
            d = dz(j-1)/2.0 + c(i,j)*dz(j)
            pf = P(i,j-1)
        else                                        !vertical
            if (sl1<0.0) then                       !dy/dx<0 fluid lies on the left of the free surface
                if ((IAND(FLAG(i-1,j),C_F)/=0).AND.(FLAG(i-1,j)<C_O)) then
                    h = (dx(i)+dx(i-1))/2.0
                    d = dx(i-1)/2.0 + c(i,j)*dx(i)
                    pf = P(i-1,j)
                else                                !hydrostatic
                    hyd=1 ! pressure is hydrostatic in the surface cell
                end if
            else                                    !dy/dx>=0 fluid lies on the right of the free surface
                if ((IAND(FLAG(i+1,j),C_F)/=0).AND.(FLAG(i+1,j)<C_O)) then
                    h = (dx(i)+dx(i+1))/2.0
                    d = dx(i+1)/2.0 + c(i,j)*dx(i)
                    pf = P(i+1,j)
                else
                    hyd=1 ! pressure is hydrostatic in the surface cell
                end if
            end if
        end if
        if (hyd==1) then 
            P_N=max(c(i,j)*dz(j)-dz(j)/2.0,0.0) ! pressure is hydrostatic in the surface cell
        else
            eta = h/d
            P_N= eta*P(i,j)+(1.0-eta)*pf
        end if
    elseif (FLAG(ii,jj-1)>=C_O) then    !south cell is a surface cell
        i=ii
        j=jj-1
        yip = c(i+1,j-1)*dz(j-1)+c(i+1,j)*dz(j)+c(i+1,j+1)*dz(j+1)
        yim = c(i-1,j-1)*dz(j-1)+c(i-1,j)*dz(j)+c(i-1,j+1)*dz(j+1)
        xjp = c(i-1,j+1)*dx(i-1)+c(i,j+1)*dx(i)+c(i+1,j+1)*dx(i+1)
        xjm = c(i-1,j-1)*dx(i-1)+c(i,j-1)*dx(i)+c(i+1,j-1)*dx(i+1)
        !DY/DX
        sl1 = 2.0*(yip-yim)/(dx(i+1)+2.0*dx(i)+dx(i-1))
        !
        !DX/DY
        sl2 = 2.0*(xjp-xjm)/(dz(j+1)+2.0*dz(j)+dz(j-1))
        if (abs(sl1)<=abs(sl2)) then                !horizontal
            h = (dz(j)+dz(j+1))/2.0
            d = dz(j+1)/2.0 + c(i,j)*dz(j)
            pf = P(i,j+1)
        else
            if (sl1<0.0) then                       !dy/dx<0 fluid lies on the left of the free surface
                if ((IAND(FLAG(i-1,j),C_F)/=0).AND.(FLAG(i-1,j)<C_O)) then 
                    h = (dx(i)+dx(i-1))/2.0
                    d = dx(i-1)/2.0 + c(i,j)*dx(i)
                    pf = P(i-1,j)
                else
                    hyd=1
                end if
            else                                    !dy/dx>=0 fluid lies on the right of the free surface
                if ((IAND(FLAG(i+1,j),C_F)/=0).AND.(FLAG(i+1,j)<C_O)) then 
                    h = (dx(i)+dx(i+1))/2.0
                    d = dx(i+1)/2.0 + c(i,j)*dx(i)
                    pf = P(i+1,j)
                 else
                    hyd=1
                end if
           end if
        end if
        if (hyd==1) then 
            P_S=max(c(i,j)*dz(j)-dz(j)/2.0,0.0) ! pressure is hydrostatic in the surface cell
        else
            eta = h/d
            P_S= eta*P(i,j)+(1.0-eta)*pf
        end if
    elseif (FLAG(ii+1,jj)>=C_O) then 
        i=ii+1
        j=jj
        yip = c(i+1,j-1)*dz(j-1)+c(i+1,j)*dz(j)+c(i+1,j+1)*dz(j+1)
        yim = c(i-1,j-1)*dz(j-1)+c(i-1,j)*dz(j)+c(i-1,j+1)*dz(j+1)
        xjp = c(i-1,j+1)*dx(i-1)+c(i,j+1)*dx(i)+c(i+1,j+1)*dx(i+1)
        xjm = c(i-1,j-1)*dx(i-1)+c(i,j-1)*dx(i)+c(i+1,j-1)*dx(i+1)
        !DY/DX
        sl1 = 2.0*(yip-yim)/(dx(i+1)+2.0*dx(i)+dx(i-1))
        !
        !DX/DY
        sl2 = 2.0*(xjp-xjm)/(dz(j+1)+2.0*dz(j)+dz(j-1))
        if (abs(sl1)<=abs(sl2)) then        !horizontal
            if (sl2<0.0) then               !dx/dy<0 fluid lies below the free surface
                if ((IAND(FLAG(i,j-1),C_F)/=0).AND.(FLAG(i,j-1)<C_O)) then 
                    h = (dz(j)+dz(j-1))/2.0
                    d = dz(j-1)/2.0 + c(i,j)*dz(j)
                    pf = P(i,j-1)
                else
                    hyd=1
                end if
            else                            !dx/dy<=0 fluid lies above the free surface
                if ((IAND(FLAG(i,j+1),C_F)/=0).AND.(FLAG(i,j+1)<C_O)) then 
                    h = (dz(j)+dz(j+1))/2.0
                    d = dz(j+1)/2.0 + c(i,j)*dz(j)
                    pf = P(i,j+1)
                else
                    hyd=1
                end if
            end if
        else                                !vertical
            h = (dx(i)+dx(i-1))/2.0
            d = dx(i-1)/2.0 + c(i,j)*dx(i)
            pf = P(i-1,j)
        end if
        if (hyd==1) then 
            P_E=max(c(i,j)*dz(j)-dz(j)/2.0,0.0) ! pressure is hydrostatic in the surface cell
        else
            eta = h/d
            P_E= eta*P(i,j)+(1.0-eta)*pf
        end if
    elseif (FLAG(ii-1,jj)>=C_O) then 
        i=ii-1
        j=jj
        yip = c(i+1,j-1)*dz(j-1)+c(i+1,j)*dz(j)+c(i+1,j+1)*dz(j+1)
        yim = c(i-1,j-1)*dz(j-1)+c(i-1,j)*dz(j)+c(i-1,j+1)*dz(j+1)
        xjp = c(i-1,j+1)*dx(i-1)+c(i,j+1)*dx(i)+c(i+1,j+1)*dx(i+1)
        xjm = c(i-1,j-1)*dx(i-1)+c(i,j-1)*dx(i)+c(i+1,j-1)*dx(i+1)
        !DY/DX
        sl1 = 2.0*(yip-yim)/(dx(i+1)+2.0*dx(i)+dx(i-1))
        !
        !DX/DY
        sl2 = 2.0*(xjp-xjm)/(dz(j+1)+2.0*dz(j)+dz(j-1))
        if (abs(sl1)<=abs(sl2)) then        !horizontal
            if (sl2<0.0) then               !dx/dy<0 fluid lies below the free surface
                if ((IAND(FLAG(i,j-1),C_F)/=0).AND.(FLAG(i,j-1)<C_O)) then 
                    h = (dz(j)+dz(j-1))/2.0
                    d = dz(j-1)/2.0 + c(i,j)*dz(j)
                    pf = P(i,j-1)
                else
                    hyd=1
                end if
            else                            !dx/dy<=0 fluid lies above the free surface
                if ((IAND(FLAG(i,j+1),C_F)/=0).AND.(FLAG(i,j+1)<C_O)) then
                    h = (dz(j)+dz(j+1))/2.0
                    d = dz(j+1)/2.0 + c(i,j)*dz(j)
                    pf = P(i,j+1)
                else
                    hyd=1
                end if
            end if
        else                                !vertical
            h = (dx(i)+dx(i+1))/2.0
            d = dx(i+1)/2.0 + c(i,j)*dx(i)
            pf = P(i+1,j)
        end if
        if (hyd==1) then 
            P_W=max(c(i,j)*dz(j)-dz(j)/2.0,0.0) ! pressure is hydrostatic in the surface cell
        else
            eta = h/d
            P_W= eta*P(i,j)+(1.0-eta)*pf
        end if
    elseif (FLAG(ii,jj)>=C_O) then    !north cell is a surface cell
        i=ii
        j=jj
        yip = c(i+1,j-1)*dz(j-1)+c(i+1,j)*dz(j)+c(i+1,j+1)*dz(j+1)
        yim = c(i-1,j-1)*dz(j-1)+c(i-1,j)*dz(j)+c(i-1,j+1)*dz(j+1)
        xjp = c(i-1,j+1)*dx(i-1)+c(i,j+1)*dx(i)+c(i+1,j+1)*dx(i+1)
        xjm = c(i-1,j-1)*dx(i-1)+c(i,j-1)*dx(i)+c(i+1,j-1)*dx(i+1)
        !DY/DX
        sl1 = 2.0*(yip-yim)/(dx(i+1)+2.0*dx(i)+dx(i-1))
        !
        !DX/DY
        sl2 = 2.0*(xjp-xjm)/(dz(j+1)+2.0*dz(j)+dz(j-1))
        if (abs(sl1)<=abs(sl2)) then        !horizontal
            if (sl2<0.0) then               !dx/dy<0 fluid lies below the free surface
                if ((IAND(FLAG(i,j-1),C_F)/=0).AND.(FLAG(i,j-1)<C_O)) then !is there a fluid cell below the surface?
                    h = (dz(j)+dz(j-1))/2.d0
                    d = dz(j-1)/2.d0 + c(i,j)*dz(j)
                    pf = P(i,j-1)
                else
                    hyd=1
                end if
            else
                if ((IAND(FLAG(i,j+1),C_F)/=0).AND.(FLAG(i,j+1)<C_O)) then
                    h = (dz(j)+dz(j+1))/2.0
                    d = dz(j+1)/2.0 + c(i,j)*dz(j)
                    pf = P(i,j+1)
                else
                    hyd=1
                end if
            end if
        else                            !vertical
            if (sl1<0.d0) then              !dy/dx<0 fluid lies on the left of the free surface
                if ((IAND(FLAG(i-1,j),C_F)/=0).AND.(FLAG(i-1,j)<C_O)) then     
                    h = (dx(i)+dx(i-1))/2.d0
                    d = dx(i-1)/2.d0 + c(i,j)*dx(i)
                    pf = P(i-1,j)
                else
                    hyd=1
                end if
            else                            !dy/dx>=0 fluid lies on the right of the free surface
                if ((IAND(FLAG(i+1,j),C_F)/=0).AND.(FLAG(i+1,j)<C_O)) then
                    h = (dx(i)+dx(i+1))/2.d0
                    d = dx(i+1)/2.d0 + c(i,j)*dx(i)
                    pf = P(i+1,j)
                else
                    hyd=1
                end if
            end if
        end if
        if (hyd==1) then 
            P(i,j)=max(c(i,j)*dz(j)-dz(j)/2.0,P(i,j)) ! pressure is hydrostatic in the surface cell  
        else
            eta = h/d
            P(i,j) = eta*P(i,j)+(1.0-eta)*pf
        end if
    else
        return
    end if
 
    return 
    end subroutine NEIGHBORINGPRESSURE
   
