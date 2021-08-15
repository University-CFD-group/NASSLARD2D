    !
    !%%%%%%%%%%%%%%%%%%%%%%%% CONV_W %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    real function CONV_WQ(i, j, imax, jmax, U, W, dx, dz)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  computes convective flux of W using QUICK
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use defs
    IMPLICIT NONE
    integer, INTENT (IN):: i, j, imax, jmax
    real(RP), DIMENSION(0:), INTENT (IN):: dx, dz
    real(RP), DIMENSION(0:, 0:), INTENT (IN)::U, W
    !
    !-----------------------------------------------------------------------
    !
    integer :: jj
    real(RP):: me, mw, mn, ms, &
    ue, uw, we, ww, wn, ws, &
    dzw, dxw, fie, fiw, fin, fis
    dzw = (dz(j) + dz(j+1))/2.0
    dxw = dx(i)
    !
    !   fluid velocities at cell faces
    !
    uw = (U(i-1,j) + U(i-1,j+1))/2.0
    ue = (U(i,j) + U(i,j+1))/2.0
    ww = (W(i,j)*dx(i-1) + W(i-1,j)*dx(i))/(dx(i)+dx(i-1))
    we = (W(i+1,j)*dx(i) + W(i,j)*dx(i+1))/(dx(i)+dx(i+1))
    wn = (W(i,j) + W(i,j+1))/2.0
    ws = (W(i,j) + W(i,j-1))/2.0
    !
    !   mass fluxes at cell faces
    !
    me = ue*dzw
    mw = -uw*dzw
    mn = -wn*dxw
    ms = ws*dxw
    !
    !   QUICK
    !
    if (me>=0.0) then
        fie = W(i,j)+(W(i+1,j)-W(i,j))/((dx(i)+dx(i+1))/2.0)*(dx(i)/2.0)+ &
        ((W(i-1,j)-W(i,j))/((-dx(i-1)-dx(i))/2.0)-(W(i+1,j)-&
        W(i,j))/((dx(i)+dx(i+1))/2.0))*(-dx(i+1)/2.0)*(dx(i)/2.0)/(-dx(i-1)/2.0-dx(i)-dx(i+1)/2.0)
    else
        if (i<imax) then
            fie = W(i+1,j)+(W(i,j)-W(i+1,j))/((-dx(i)-dx(i+1))/2.0)*(-dx(i+1)/2.0)+ &
            ((W(i+2,j)-W(i+1,j))/((dx(i+1)+dx(i+2))/2.0)-(W(i,j)-&
            W(i+1,j))/((-dx(i)-dx(i+1))/2.0))*(-dx(i+1)/2.0)*(dx(i)/2.0)/(dx(i)/2.0+dx(i+1)+dx(i+2)/2.0)
        else
            !FOU
            fie = W(i+1,j)
        end if
    end if
    if (mw<=0.0) then
        if (i>1) then
            fiw = W(i-1,j)+(W(i,j)-W(i-1,j))/((dx(i)+dx(i-1))/2.0)*(dx(i-1)/2.0)+ &
            ((W(i-2,j)-W(i-1,j))/((-dx(i-2)-dx(i-1))/2.0)-(W(i,j)-&
            W(i-1,j))/((dx(i)+dx(i-1))/2.0))*(dx(i-1)/2.0)*(-dx(i)/2.0)/(-dx(i-2)/2.0-dx(i-1)-dx(i)/2.0)
        else
            !FOU
            fiw = W(i-1,j)
        end if
    else
        fiw = W(i,j)+(W(i-1,j)-W(i,j))/((-dx(i)-dx(i-1))/2.0)*(-dx(i)/2.0)+ &
        ((W(i+1,j)-W(i,j))/((dx(i)+dx(i+1))/2.0)-(W(i-1,j)-&
        W(i,j))/((-dx(i)-dx(i-1))/2.0))*(dx(i-1)/2.0)*(-dx(i)/2.0)/(dx(i-1)/2.0+dx(i)+dx(i+1)/2.0)
    end if
    if (mn<=0.0) then
        fin = W(i,j)+(W(i,j+1)-W(i,j))/dz(j)*(dz(j)/2.0)+ &
        ((W(i,j-1)-W(i,j))/(-dz(j))-(W(i,j+1)-&
        W(i,j))/(dz(j+1)))*(dz(j+1)/2.0)*(-dz(j+1)/2.0)/(-dz(j)-dz(j+1))
    else
        if (j<jmax-1) then
            fin = W(i,j+1)+(W(i,j)-W(i,j+1))/(-dz(j+1))*(-dz(j+1)/2.0)+ &
            ((W(i,j+2)-W(i,j+1))/(dz(j+2))-(W(i,j)-&
            W(i,j+1))/(-dz(j+1)))*(dz(j+1)/2.0)*(-dz(j+1)/2.0)/(dz(j+1)+dz(j+2))
        else
            !FOU
            fin = W(i,j+1)
        end if
    end if
    if (ms>=0.0) then
        if (j>1) then
            fis = W(i,j-1)+(W(i,j)-W(i,j-1))/dz(j)*(dz(j)/2.0)+ &
            ((W(i,j-2)-W(i,j-1))/(-dz(j-1))-(W(i,j)-&
            W(i,j-1))/(dz(j)))*(dz(j)/2.0)*(-dz(j)/2.0)/(-dz(j-1)-dz(j))
        else
            !FOU
            fis = W(i,j-1)
        end if
    else
        fis = W(i,j)+(W(i,j-1)-W(i,j))/(-dz(j))*(-dz(j)/2.0)+ &
        ((W(i,j+1)-W(i,j))/(dz(j+1))-(W(i,j-1)-&
        W(i,j))/(-dz(j)))*(dz(j)/2.0)*(-dz(j)/2.0)/(dz(j+1)+dz(j))
    end if
    !
    !   convective fluxes at cell faces
    !
    CONV_WQ = fie*me + fiw*mw - fin*mn - fis*ms
    return
    end function CONV_WQ
    !
    !%%%%%%%%%%%%%%%%%%%%%%%% DifF_W %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    real function DifF_W(i, j, W, dx, dz)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !   computes diffusive flux of W using second order polynomial approximation
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    IMPLICIT NONE
    integer, INTENT(IN):: i, j
    real(RP), DIMENSION(0:), INTENT(IN):: dx, dz
    real(RP), DIMENSION(0:,0:), INTENT(IN):: W
    !
    !-----------------------------------------------------------------------
    !
    real(RP):: btx, btz, ax, az,dxe
    dxe = (dx(i)+dx(i+1))/2.0
    btx = (dx(i)+dx(i-1))/(dx(i)+dx(i+1))
    btz = dz(j)/dz(j+1)
    ax = (W(i-1,j)-(1.0+btx)*W(i,j)+btx*W(i+1,j))/btx/(1.0+btx)/dxe**2
    az = (W(i,j-1)-(1.0+btz)*W(i,j)+btz*W(i,j+1))/btz/(1.0+btz)/dz(j+1)**2
    DifF_W = 2.0*dx(i)*(dz(j)+dz(j+1))/2.0*(ax+az)
    return
    end function DifF_W
    !
    !%%%%%%%%%%%%%%%%%%%%%%%% SRC_W %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    real function SRC_W(i, j, GZ, dx, dz)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !   computes source of W
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    use nrtype
    IMPLICIT NONE
    integer, INTENT(IN)::i, j
    real(RP), INTENT(IN):: GZ
    real(RP), DIMENSION(0:), INTENT(IN):: dx, dz
    !
    !-----------------------------------------------------------------------
    !
    SRC_W = GZ*dx(i)*(dz(j) + dz(j+1))/2.0
    return
    end function SRC_W
    !
    !%%%%%%%%%%%%%%%%%%%%%%%% CONV_W %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    real function CONV_W(i, j, U, W, dx, dz, gamma)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  computes convective flux of W
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    IMPLICIT NONE
    integer, INTENT (IN):: i, j
    real(RP), DIMENSION(0:), INTENT (IN):: dx, dz
    real(RP), DIMENSION(0:, 0:), INTENT (IN)::U, W
    real(RP), INTENT(IN)::gamma
    !
    !-----------------------------------------------------------------------
    !
    real(RP):: me, mw, mn, ms, &
    ue, uw, we, ww, wn, ws, &
    convwe, convww, convwn, convws, dzw, dxw
    dzw = (dz(j) + dz(j+1))/2.0
    dxw = dx(i)
    !
    !   fluid velocities at cell faces
    !
    ue = (U(i,j) + U(i,j+1))/2.0
    uw = (U(i-1,j) + U(i-1,j+1))/2.0
    ww = (W(i,j)*dx(i-1) + W(i-1,j)*dx(i))/(dx(i)+dx(i-1))
    we = (W(i+1,j)*dx(i) + W(i,j)*dx(i+1))/(dx(i)+dx(i+1))
    wn = (W(i,j) + W(i,j+1))/2.0
    ws = (W(i,j) + W(i,j-1))/2.0
    !
    !   mass fluxes at cell faces
    !
    me = ue*dzw
    mw = -uw*dzw
    mn = -wn*dxw
    ms = ws*dxw
    !
    !   convective fluxes at cell faces
    !
    convwe = gamma*(W(i,j)*max(me,0.0) + W(i+1,j)*min(me,0.0)) + &
    (1.0-gamma)*me*we
    convww = gamma*(W(i,j)*max(mw,0.0) + W(i-1,j)*min(mw,0.0)) + &
    (1.0-gamma)*mw*ww
    convwn = gamma*(W(i,j)*min(mn,0.0) + W(i,j+1)*max(mn,0.0)) + &
    (1.0-gamma)*mn*wn
    convws = gamma*(W(i,j)*min(ms,0.0) + W(i,j-1)*max(ms,0.0)) + &
    (1.0-gamma)*ms*ws
    CONV_W = convwe + convww - convwn - convws
    return
    end function CONV_W
    !
    !%%%%%%%%%%%%%%%%%%%%%%%% COMP_Z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    real function COMP_Z(i, j, U, W, dx, dz)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !   computes diffusive flux of W using second order polynomial approximation
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    IMPLICIT NONE
    integer, INTENT(IN):: i, j
    real(RP), DIMENSION(0:), INTENT(IN):: dx, dz
    real(RP), DIMENSION(0:,0:), INTENT(IN):: U, W
    !
    !-----------------------------------------------------------------------
    !
    real(RP):: DN, DS
    DN = (U(i,j+1)-U(i-1,j+1))/dx(i) + (W(i,j+1)-W(i,j))/dz(j+1)
    DS = (U(i,j)-U(i-1,j))/dx(i) + (W(i,j)-W(i,j-1))/dz(j)
    COMP_Z = (DN-DS)*dx(i)
    return
    end function COMP_Z
