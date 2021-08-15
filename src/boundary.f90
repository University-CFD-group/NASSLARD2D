    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETBCOND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine SETBCOND (problem, U, W, uflag, vflag, pflag, imax, jmax, ppu, ppv, wW, wE, wN, wS, &
                         indu, indv, induin, indvin, gFIu, gFIv, gFIuin, gFIvin, invagu, invagv, &
                            invaguin, invagvin, ghu, ghv, ghuin, ghvin, IMu, IMv, ivp, jvp, uxb, uyb)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  Setting the boundary conditions at the boundary strip.
    !  The flags wW,wE,wNN, and wS can have the values:
    !
    !  1 = slip               2 = no-slip
    !  3 = outflow
    !
    !  Moreover, no-slip conditions are set at internal obstacle cells
    !  by default.1
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use defs
    implicit none
    character (LEN=30) :: problem
    INTEGER, intent(IN) :: imax, jmax, wW, wE, wN, wS, ppu, ppv, ivp, jvp
    INTEGER(I2B), dimension(0:,0:), intent(IN) :: uflag, vflag, pflag
    real(RP), dimension(0:,0:), intent(INOUT) :: U, W
    integer, dimension(:,:), intent(IN) :: indu, indv, induin, indvin
    integer, dimension(:,:,:), intent(IN) :: gFIu, gFIv, gFIuin, gFIvin
    real(DP), dimension(:,:,:), intent(IN) :: invagu, invagv, invaguin, invagvin
    real(DP), dimension(:,:), intent(IN) :: ghu, ghv, ghuin, ghvin
    integer, dimension(:), intent(IN) :: IMu, IMv
    real(RP), intent(IN) :: uxb, uyb
    !
    !-----------------------------------------------------------------------
    !
    integer:: i, j, iter, ii, jj, i1, j1, i2, j2, k, m
    real(RP):: zz, uto, E, kappa, cmu, y1, ypl, x1, xpl, fs, fsd, uto0, uton, vel, res, eps, &
               tow, yp, xp, y2, x2, h1, h2
    
    real(DP) a0, a1, a2, fd, fl, cd, cl, fiu, fiv, um, u1, u2, u3, u4, uxvp, uyvp, vp, ub, xin, yin
    real(DP), dimension(3):: fi
    E = 9.0
    kappa = 0.41
    cmu=0.09
    ub=0.0d0     !solid boundary velocity
    do j = 1, jmax
        !
        ! western and eastern boundary
        !
        if ( wW == 1 ) then                             ! free-slip
            U(0,j) = 0.0			                      ! u = 0
            W(0,j) = W(1,j)                             ! dv/dn = 0
         else if (wW == 2 ) then                         ! no-slip
            U(0,j) = 0.0			                      ! u = 0
            W(0,j) = (-1.0)*W(1,j)                      ! w=0 at the boundary by averaging
        else if (wW == 3) then                          ! outflow
            U(0,j) = U(1,j)
            W(0,j) = W(1,j)
        end if
        if (wE == 1) then                               ! free-slip
            U(imax,j) = 0.0
            W(imax+1,j) = W(imax,j)                     ! dw/dn = 0
        else if (wE == 2 ) THEN                         ! no-slip
            U(imax,j) = 0.0
            W(imax+1,j) = (-1.0)*W(imax,j)              ! w=0 at the boundary by averaging
        else if (wE == 3) THEN                           ! outflow
            U(imax,j) = U(imax-1,j)                     ! du/dn = 0
            W(imax+1,j) = W(imax,j)                     ! dw/dn = 0
        end if
    end do

    do i=1, imax
        !
        ! northern and southern boundary
        !
        if (wN == 1 ) then
            W(i,jmax) = 0.0
            U(i,jmax+1) = U(i,jmax)
        else if (wN == 2 ) then
            W(i,jmax) = 0.0
            U(i,jmax+1) = (-1.0)*U(i,jmax)
        else if (wN == 3) THEN
            W(i,jmax) = W(i,jmax-1)
            U(i,jmax+1) = U(i,jmax)
        end if 

        if (wS == 1 ) then
            W(i,0) = 0.0
            U(i,0) = U(i,1)
        else if (wS == 2 ) then
            W(i,0) = 0.0
            U(i,0) = (-1.0)*U(i,1)
        else if (wS == 3) then
            W(i,0) = W(i,1)
            U(i,0) = U(i,1)
        end if
    end do

    if (immersed) then 
        do i=1, ivp
            ii=induin(i,1)
            jj=induin(i,2)
            i1=gFIuin(i,1,1)          !i index of first interpolation point
            j1=gFIuin(i,1,2)          !j index of first interpolation point
            i2=gFIuin(i,2,1)          !i index of second interpolation point
            j2=gFIuin(i,2,2)          !j index of second interpolation point
            fi(1)=uxb                 !horizontal solid velocity 
            fi(2)=U(i1,j1)            !velocity at first interpolation point
            fi(3)=U(i2,j2)            !velocity at second interpolation point
            !compute linear interpolation constants a0, a1 and a2
            a0=fi(1)+fi(2)*invaguin(i,1,2)+fi(3)*invaguin(i,1,3)
            a1=fi(1)+fi(2)*invaguin(i,2,2)+fi(3)*invaguin(i,2,3)
            a2=fi(1)+fi(2)*invaguin(i,3,2)+fi(3)*invaguin(i,3,3)
            !compute u variable at ghost point using linear interpolation
            U(ii,jj)=a0+a1*ghuin(i,1)+a2*ghuin(i,2)
        end do
        do i=1, jvp
            ii=indvin(i,1)
            jj=indvin(i,2)
            i1=gFIvin(i,1,1)          !i index of first interpolation point
            j1=gFIvin(i,1,2)          !j index of first interpolation point
            i2=gFIvin(i,2,1)          !i index of second interpolation point
            j2=gFIvin(i,2,2)          !j index of second interpolation point
            fi(1)=uyb                 !vertical solid velocity 
            fi(2)=W(i1,j1)            !velocity at first interpolation point
            fi(3)=W(i2,j2)            !velocity at second interpolation point
            !compute linear interpolation constants a0, a1 and a2
            a0=fi(1)+fi(2)*invagvin(i,1,2)+fi(3)*invagvin(i,1,3)
            a1=fi(1)+fi(2)*invagvin(i,2,2)+fi(3)*invagvin(i,2,3)
            a2=fi(1)+fi(2)*invagvin(i,3,2)+fi(3)*invagvin(i,3,3)
            !compute u variable at ghost point using linear interpolation
            W(ii,jj)=a0+a1*ghvin(i,1)+a2*ghvin(i,2)
        end do

        do i=1, ppu
            ii=indu(i,1)
            jj=indu(i,2)
            i1=gFIu(i,1,1)          !i index of first interpolation point
            j1=gFIu(i,1,2)          !j index of first interpolation point
            i2=gFIu(i,2,1)          !i index of second interpolation point
            j2=gFIu(i,2,2)          !j index of second interpolation point
            fi(1)=uxb               !horizontl solid velocity 
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
            U(ii,jj)=fiu
        end do

        do i=1, ppv
            ii=indv(i,1)
            jj=indv(i,2)
            i1=gFIv(i,1,1)          !i index of first interpolation point
            j1=gFIv(i,1,2)          !j index of first interpolation point
            i2=gFIv(i,2,1)          !i index of second interpolation point
            j2=gFIv(i,2,2)          !j index of second interpolation point
            fi(1)=uyb               !vertical solid velocity 
            fi(2)=W(i1,j1)          !velocity at first interpolation point
            fi(3)=W(i2,j2)          !velocity at second interpolation point
            !compute linear interpolation constants a0, a1 and a2
            a0=fi(1)+fi(2)*invagv(i,1,2)+fi(3)*invagv(i,1,3)
            a1=fi(1)+fi(2)*invagv(i,2,2)+fi(3)*invagv(i,2,3)
            a2=fi(1)+fi(2)*invagv(i,3,2)+fi(3)*invagv(i,3,3)
            !compute u variable at ghost point using linear interpolation
            fiv=a0+a1*ghv(i,1)+a2*ghv(i,2)
            !if imagine of ghost point is used then re-compute the ghost point value
            if (IMv(i)==1) fiv=2.0d0*fi(1)-fiv
            W(ii,jj)=fiv
        end do
    end if
    
    return
    end subroutine SETBCOND
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%% SETSPECBCOND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine SETSPECBCOND (problem, dz, U, W, imax, jmax, um, cfreq, t, c, G, P, Fr, delt, q, St)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  Setting specific boundary conditions, depending on "problem" 
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use defs
    implicit none
    character (LEN=30) :: problem
    real(RP), dimension(0:), intent(IN) :: dz
    real(RP), dimension(0:,0:), intent(INOUT) :: U, W
    integer, intent(IN) :: imax, jmax
    real(RP), intent(IN) :: um, cfreq, t, Fr, delt, q, St
    real(DP), dimension(0:,0:), intent(IN) :: c
    real(RP), dimension(0:,0:), intent(IN) :: G, P
    !
    !  local variables
    !
    integer :: i,j
    real(RP):: zz, UI, KAI, EPI, hostep, up, dzz, pbn, pbnp1

    if ( (problem=="drop") .OR. (problem=="dam") ) then
        return
        !-----------------------------------------------------------
        ! Driven Cavity: U = um at the upper boundary              
        !-----------------------------------------------------------

    else if  (problem=="dcavity") then
        do i = 0, imax
            U(i,jmax+1) = 2.0*um - U(i,jmax)
        end do
    elseif (problem=="circle") then
        if (.NOT.movbound) then
            do j=1, jmax
                U(0,j)=um
                W(0,j)=0.0
            end do
        end if
    else if ( problem=="backstep" ) then

        !*------------------*/
        zz = dz(1)/2.0
        do j = 1, jmax
            if (zz>=0.5) then
                U(0,j)=1.0
            end if
            zz = zz + (dz(j)+ dz(j+1))/2.0
        end do
    elseif (problem=="wavef") then 
        up = um*sin(cfreq*t)
        do j=1, jmax
            if (c(1,j)/=0.0d0) then 
                U(0,j) = up
            end if
        end do
    elseif (problem=="reservoir") then 
        if (absorption) then 
            dzz=(dz(0)+dz(1))/2.0
            do i=1, imax
                W(i,0)=G(i,0)-St*delt/Fr**2*(P(i,1)-P(i,0))/dzz-St*q*(P(i,1)+P(i,0))/2.0
            end do
        end if
    end if
    
    return
    end subroutine SETSPECBCOND
    
   
    
   

