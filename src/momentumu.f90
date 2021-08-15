    !
    !%%%%%%%%%%%%%%%%%%%%%%%% CONV_UQ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    real function CONV_UQ(i, j, imax, jmax, U, W, dx, dz)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  computes convective flux of U using QUICK
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use defs
    IMPLICIT NONE
    integer, INTENT (IN):: i, j, imax, jmax
    real(RP), DIMENSION(0:), INTENT (IN):: dx, dz
    real(RP), DIMENSION(0:,0:), INTENT(IN):: U, W
    !
    !-----------------------------------------------------------------------
    !
    integer :: jj
    real(RP):: me, mw, mn, ms, &
    ue, uw, un, us, wn, ws, &
    dzu, dxu, fie, fiw, fin, fis
    dzu = dz(j)
    dxu = (dx(i)+dx(i+1))/2.0
    !
    !   fluid velocities at cell faces
    !
    ue = (U(i,j) + U(i+1,j))/2.0
    uw = (U(i,j) + U(i-1,j))/2.0
    un = (U(i,j+1)*dz(j) + U(i,j)*dz(j+1))/(dz(j)+dz(j+1))
    us = (U(i,j)*dz(j-1) + U(i,j-1)*dz(j))/(dz(j)+dz(j-1))
    wn = (W(i+1,j) + W(i,j))/2.0
    ws = (W(i,j-1) + W(i+1,j-1))/2.0
    !
    !   mass fluxes at cell faces
    !
    me = ue*dzu
    mw = -uw*dzu
    mn = -wn*dxu
    ms = ws*dxu
    !
    !   QUICK
    !
    if (me>=0.0) then
        fie = U(i,j)+(U(i+1,j)-U(i,j))/dx(i+1)*(dx(i+1)/2.0)+ &
        ((U(i-1,j)-U(i,j))/(-dx(i))-(U(i+1,j)-&
        U(i,j))/dx(i+1))*(dx(i+1)/2.0)*(-dx(i+1)/2.0)/(-dx(i)-dx(i+1))
    else
        if (i<imax-1) then
            fie = U(i+1,j)+(U(i,j)-U(i+1,j))/(-dx(i+1))*(-dx(i+1)/2.0)+ &
            ((U(i+2,j)-U(i+1,j))/(dx(i+2))-(U(i,j)-&
            U(i+1,j))/(-dx(i+1)))*(-dx(i+1)/2.0)*(dx(i+1)/2.0)/(dx(i+1)+dx(i+2))
        else
            !FOU
            fie = U(i+1,j)
        end if
    end if

    if (mw<=0.0) then
        if (i>1) then
            fiw = U(i-1,j)+(U(i,j)-U(i-1,j))/(dx(i))*(dx(i)/2.0)+ &
            ((U(i-2,j)-U(i-1,j))/(-dx(i-1))-(U(i,j)-&
            U(i-1,j))/dx(i))*(dx(i)/2.0)*(-dx(i)/2.0)/(-dx(i-1)-dx(i))
        else
            !FOU
            fiw = U(i-1,j)
        end if
    else
        fiw = U(i,j)+(U(i-1,j)-U(i,j))/(-dx(i))*(-dx(i)/2.0)+ &
        ((U(i+1,j)-U(i,j))/(dx(i+1))-(U(i-1,j)-&
        U(i,j))/(-dx(i)))*(dx(i)/2.0)*(-dx(i)/2.0)/(dx(i)+dx(i+1))
    end if
    if (mn<=0.0) then
        fin = U(i,j)+(U(i,j+1)-U(i,j))/((dz(j)+dz(j+1))/2.0)*(dz(j)/2.0)+ &
        ((U(i,j-1)-U(i,j))/((-dz(j)-dz(j-1))/2.0)-(U(i,j+1)-&
        U(i,j))/((dz(j)+dz(j+1))/2.0))*(dz(j)/2.0)*(-dz(j+1)/2.0)/(-dz(j-1)/2.0-dz(j)-dz(j+1)/2.0)
    else
        if (j<jmax) then
            fin = U(i,j+1)+(U(i,j)-U(i,j+1))/((-dz(j)-dz(j+1))/2.0)*(-dz(j)/2.0)+ &
            ((U(i,j+2)-U(i,j+1))/((dz(j+1)+dz(j+2))/2.0)-(U(i,j)-&
            U(i,j+1))/((-dz(j)-dz(j+1))/2.0))*(dz(j)/2.0)*(-dz(j+1)/2.0)/(dz(j+2)/2.0+dz(j+1)+dz(j)/2.0)
        else
            !FOU
            fin = U(i,j+1)
        end if
    end if
    if (ms>=0.0) then
        if (j>1) then
            fis = U(i,j-1)+(U(i,j)-U(i,j-1))/((dz(j)+dz(j-1))/2.0)*(dz(j-1)/2.0)+ &
            ((U(i,j-2)-U(i,j-1))/((-dz(j-1)-dz(j-2))/2.0)-(U(i,j)-&
            U(i,j-1))/((dz(j-1)+dz(j))/2.0))*(-dz(j)/2.0)*(dz(j-1)/2.0)/(-dz(j-2)/2.0-dz(j-1)-dz(j)/2.0)
        else
            !FOU
            fis = U(i,j-1)
        end if
    else
        fis = U(i,j)+(U(i,j-1)-U(i,j))/((-dz(j-1)-dz(j))/2.0)*(-dz(j)/2.0)+ &
        ((U(i,j+1)-U(i,j))/((dz(j)+dz(j+1))/2.0)-(U(i,j-1)-&
        U(i,j))/((-dz(j-1)-dz(j))/2.0))*(-dz(j)/2.0)*(dz(j-1)/2.0)/(dz(j+1)/2.0+dz(j)+dz(j-1)/2.0)
    end if
    !
    !   upwind differencing and central differencing
    !
    CONV_UQ = fie*me + fiw*mw - fin*mn - fis*ms
    return
    end function CONV_UQ
    !
    !%%%%%%%%%%%%%%%%%%%%%%%% DifF_U %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    real function DifF_U(i, j, U, dx, dz)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  computes diffusive flux of U using second order polynomial approximation
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    IMPLICIT NONE
    integer, INTENT(IN):: i, j
    real(RP), DIMENSION(0:), INTENT(IN):: dx, dz
    real(RP), DIMENSION(0:, 0:), INTENT(IN):: U
    !
    !-----------------------------------------------------------------------
    !
    real(RP):: btx, btz, ax, az, dzn
    dzn = (dz(j)+dz(j+1))/2.0
    btx = dx(i)/dx(i+1)
    btz = (dz(j)+dz(j-1))/(dz(j)+dz(j+1))
    ax = (U(i-1,j)-(1.0+btx)*U(i,j)+btx*U(i+1,j))/btx/(1.0+btx)/dx(i+1)**2
    az = (U(i,j-1)-(1.0+btz)*U(i,j)+btz*U(i,j+1))/btz/(1.0+btz)/dzn**2
    DifF_U = 2.0*(dx(i)+dx(i+1))/2.0*dz(j)*(ax+az)
    return
    end function DifF_U
    !
    !%%%%%%%%%%%%%%%%%%%%%%%% CUNV_U %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    real function CONV_U(i, j, U, W, dx, dz, gamma)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  computes convective flux of U
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    IMPLICIT NONE
    integer, INTENT (IN):: i, j
    real(RP), DIMENSION(0:), INTENT (IN):: dx, dz
    real(RP), DIMENSION(0:,0:), INTENT(IN):: U, W
    real(RP), INTENT(IN)::gamma
    !
    !-----------------------------------------------------------------------
    !
    real(RP):: me, mw, mn, ms, &
    ue, uw, un, us, wn, ws, &
    convue, convuw, convun, convus, dzu, dxu
    dzu = dz(j)
    dxu = (dx(i)+dx(i+1))/2.0
    !
    !   fluid velocities at cell faces
    !
    ue = (U(i,j) + U(i+1,j))/2.0
    uw = (U(i,j) + U(i-1,j))/2.0
    un = (U(i,j+1)*dz(j) + U(i,j)*dz(j+1))/(dz(j)+dz(j+1))
    us = (U(i,j)*dz(j-1) + U(i,j-1)*dz(j))/(dz(j)+dz(j-1))
    wn = (W(i+1,j) + W(i,j))/2.0
    ws = (W(i,j-1) + W(i+1,j-1))/2.0
    !
    !   mass fluxes at cell faces
    !
    me = ue*dzu
    mw = -uw*dzu
    mn = -wn*dxu
    ms = ws*dxu
    !
    !   upwind differencing and central differencing
    !
    convue = gamma*(U(i,j)*max(me,0.0) + U(i+1,j)*min(me,0.0)) + &
    (1.0-gamma)*me*ue
    convuw = gamma*(U(i,j)*max(mw,0.0) + U(i-1,j)*min(mw,0.0)) + &
    (1.0-gamma)*mw*uw
    convun = gamma*(U(i,j)*min(mn,0.0) + U(i,j+1)*max(mn,0.0)) + &
    (1.0-gamma)*mn*un
    convus = gamma*(U(i,j)*min(ms,0.0) + U(i,j-1)*max(ms,0.0)) + &
    (1.0-gamma)*ms*us
    CONV_U = convue + convuw - convun - convus
    return
    end function CONV_U
    !
    !%%%%%%%%%%%%%%%%%%%%%%%% COMP_X %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    real function COMP_X(i, j, U, W, dx, dz)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  compute compressibility term for u_momentum equation
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    IMPLICIT NONE
    integer, INTENT(IN):: i, j
    real(RP), DIMENSION(0:), INTENT(IN):: dx, dz
    real(RP), DIMENSION(0:, 0:), INTENT(IN):: U, W
    !
    !-----------------------------------------------------------------------
    !
    real(RP):: DE, DW
    DE = (U(i+1,j)-U(i,j))/dx(i+1) + (W(i+1,j)-W(i+1,j-1))/dz(j)
    DW = (U(i,j)-U(i-1,j))/dx(i) + (W(i,j)-W(i,j-1))/dz(j)
    COMP_X = (DE-DW)*dz(j)
    return
    end function COMP_X
