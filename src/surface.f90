    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% J_OF_Z %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    integer function J_OF_Z(jmax, h, dz)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  Gives the j index of given h depth
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    implicit none
    integer, INTENT(IN):: jmax
    real(RP), INTENT(IN):: h
    real(RP), DIMENSION(0:), INTENT(IN):: dz
    !
    !-----------------------------------------------------------------------
    !
    integer:: j
    real(RP):: zS, zN
    zS = 0.0
    if (h == 0.0) then
        J_OF_Z = 0
        return
    end if
    do j=1, jmax
        zN = zS + dz(j)
        if ((h>=zS) .AND. (h<zN)) then
            J_OF_Z = j
            return
        else
            zS = zN
        end if
    end do
    return
    end function J_OF_Z
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MARK_CELLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine MARK_CELLS (FLAG, imax, jmax, dz, ifull, H)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  Mark the cells of the fluid domain
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use interfaces2
    use defs
    implicit none
    integer, INTENT(IN) :: imax, jmax
    integer, INTENT(OUT) :: ifull
    real(RP), DIMENSION(0:), INTENT(IN):: dz
    integer(I2B), DIMENSION(0:,0:), INTENT(INOUT) :: FLAG
    real(RP), DIMENSION(0:), INTENT(IN) :: H
    !
    !-----------------------------------------------------------------------
    !
    integer :: i, j, jj
    !
    !  Set all cells which are not obstacle cells to empty cells
    !
    do i = 0,imax+1
        do j = 0,jmax+1
            if (FLAG(i,j) >= C_F) then
                FLAG(i,j) = IAND(IOR(FLAG(i,j),C_E),NOT(C_NSWO))
            end if
        end do
    end do
    !
    !  Mark cells containing particles as fluid cells (loop over particles)
    !
    do i = 1,imax+1
        jj = J_OF_Z(jmax, H(i), dz)
        do j = 0, jj
            if (.NOT.(FLAG(i,j)<C_F)) then
                FLAG(i,j) = IAND(FLAG(i,j),NOT(C_E))
            end if
        end do
    end do

    !
    ! Mark surface cells
    !
    ifull = 0
    do j = 1,jmax
        do i = 1,imax
            if ( (IAND(FLAG(i,j),C_F) /= 0) .AND. (FLAG(i,j) < C_E) ) then
                if (IAND(FLAG(i-1,j),C_E) /= 0) then
                    FLAG(i,j) = IOR(FLAG(i,j),C_W)
                end if
                if (IAND(FLAG(i+1,j),C_E) /= 0) then
                    FLAG(i,j) = IOR(FLAG(i,j),C_O)
                end if
                if (IAND(FLAG(i,j+1),C_E) /= 0) then
                    FLAG(i,j) = IOR(FLAG(i,j),C_N)
                end if
                if (FLAG(i,j) < C_O) ifull = ifull + 1
            end if
        end do
    end do
    return
    end subroutine MARK_CELLS

    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET_UVP_SURFACE %%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine SET_UVP_SURFACE_TUMMAC (U, W, P, FLAG, imax, jmax, Fr, Re, dx, dz)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    ! Set boundary values at free surface using TUMMAC-V
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use defs
    implicit none
    integer, INTENT(IN) :: imax, jmax
    real(RP), INTENT(IN) :: Fr, Re
    integer(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
    real(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, W, P
    real(RP), DIMENSION(0:), INTENT(IN) :: dx, dz
    !
    !-----------------------------------------------------------------------
    !
    integer :: i, j
    real(RP):: dzn, dzs, dxw, dxe
    !
    ! Set velocity values in empty cells to zero
    !
    do j = 1,jmax
        do i = 1,imax-1
            if ( (IAND(FLAG(i,j),  C_E) /= 0) .AND. &
            (IAND(FLAG(i+1,j),C_E) /= 0) ) then
                U(i,j) = 0.0
            end if
        end do
    end do

    do  j = 1,jmax-1
        do i = 1,imax
            if ( (IAND(FLAG(i,j),  C_E) /= 0) .AND. &
            (IAND(FLAG(i,j+1),C_E) /= 0)) then
                W(i,j) = 0.0
            end if
        end do
    end do

    !
    !Free surface boundary conditions for u_velocity
    !
    do j = 1,jmax
        do i = 1,imax
           
            
            !
            ! treat only surface cells
            !
            if ( .NOT.(IAND(FLAG(i,j),  C_E) /= 0) .OR. (FLAG(i,j) < C_O) ) then
                 
                SELECT CASE (IAND(FLAG(i,j),C_NSWO)) ! mask NSWO_E=0x0f00
                  
                    ! filters surface cells
                CASE (C_N)
                    !Case 1
                    if (((IAND(FLAG(i-1,j),C_F)/=0).AND.(FLAG(i-1,j)<C_O)) .AND. &
                    ((IAND(FLAG(i+1,j),C_F)/=0).AND.(FLAG(i+1,j)<C_O))) then
                        dxw = (dx(i-1)+dx(i))/2.0
                        dxe = (dx(i)+dx(i+1))/2.0
                        W(i,j) = (W(i-1,j)*dxe+W(i+1,j)*dxw)/(dxe+dxw)
                        !Case 2
                    ELSE if (((IAND(FLAG(i-1,j),C_F)/=0).AND.(FLAG(i-1,j)<C_O)) .AND. &
                    (IAND(FLAG(i+1,j+1),C_E) /= 0)) then
                        dxw = (dx(i-1)+dx(i))/2.0
                        W(i,j) = (dz(j)*W(i-1,j)+dxw*W(i,j-1))/(dxw+dz(j))
                        !Case 3
                    ELSE if ((IAND(FLAG(i-1,j+1),C_E) /= 0) .AND. &
                    ((IAND(FLAG(i+1,j),C_F)/=0).AND.(FLAG(i+1,j)<C_O))) then
                        dxe = (dx(i)+dx(i+1))/2.0
                        W(i,j) = (dz(j)*W(i+1,j)+dxe*W(i,j-1))/(dxe+dz(j))
                        !Case 4
                    ELSE
                        W(i,j) = W(i,j-1)
                    end if
                    !Shear velocities
                    if (IAND(FLAG(i-1,j+1),C_E) /= 0) then
                        dzn = (dz(j)+dz(j+1))/2.0
                        dxw = (dx(i)+dx(i-1))/2.0
                        U(i-1,j+1) = U(i-1,j)-dzn/dxw*(W(i,j)-W(i-1,j))
                    end if
                CASE (C_O)
                    !Case 1
                    U(i,j) = U(i-1,j)
                    !Shear velocities
                    if (IAND(FLAG(i+1,j-1),C_E) /= 0) then
                        dzs = (dz(j)+dz(j-1))/2.0
                        dxe = (dx(i)+dx(i+1))/2.0
                        W(i+1,j-1) = W(i,j-1)-dxe/dzs*(U(i,j)-U(i,j-1))
                    end if
                CASE (C_W)
                    U(i-1,j) = U(i,j)
                    !Shear velocities
                    if (IAND(FLAG(i-1,j-1),C_E) /= 0) then
                        dzs = (dz(j)+dz(j-1))/2.0
                        dxw = (dx(i)+dx(i-1))/2.0
                        W(i-1,j-1) = W(i,j-1)+dxw/dzs*(U(i-1,j)-U(i-1,j-1))
                    end if
                CASE (C_WO)
                    !Case 5-6
                    if (((IAND(FLAG(i+1,j-1),  C_F)/=0).AND.(FLAG(i+1,j-1)  <C_E))) then
                        U(i,j) = U(i,j-1)
                    ELSE
                        U(i,j) = U(i,j)
                    end if
                    if (((IAND(FLAG(i-1,j-1),  C_F)/=0).AND.(FLAG(i-1,j-1)  <C_E))) then
                        U(i-1,j) = U(i-1,j-1)
                    ELSE
                        U(i-1,j) = U(i-1,j)
                    end if
                    !Shear velocities
                    if (IAND(FLAG(i-1,j-1),C_E) /= 0) then
                        dzs = (dz(j)+dz(j-1))/2.0
                        dxw = (dx(i)+dx(i-1))/2.0
                        W(i-1,j-1) = W(i,j-1)+dxw/dzs*(U(i-1,j)-U(i-1,j-1))
                    end if
                    if (IAND(FLAG(i+1,j-1),C_E) /= 0) then
                        dzs = (dz(j)+dz(j-1))/2.0
                        dxe = (dx(i)+dx(i+1))/2.0
                        W(i+1,j-1) = W(i,j-1)-dxe/dzs*(U(i,j)-U(i,j-1))
                    end if
                CASE (C_NO)
                    !Case 1
                    if ((IAND(FLAG(i+1,j-1),C_F)/=0).AND.(FLAG(i+1,j-1)<C_E)) then
                        dzs = (dz(j)+dz(j-1))/2.0
                        U(i,j) = (dzs*U(i-1,j)+dx(i)*U(i,j-1))/(dx(i)+dzs)
                        !Case 2
                    ELSE
                        U(i,j) = U(i-1,j)
                    end if
                    !Case 2
                    if ((IAND(FLAG(i-1,j),C_F)/=0).AND.(FLAG(i-1,j)<C_O)) then
                        dxw = (dx(i-1)+dx(i))/2.0
                        W(i,j) = (dz(j)*W(i-1,j)+dxw*W(i,j-1))/(dxw+dz(j))
                        !Case 4
                    ELSE
                        W(i,j) = W(i,j-1)
                    end if
                    !Shear velocities
                    if (IAND(FLAG(i-1,j+1),C_E) /= 0) then
                        dzn = (dz(j)+dz(j+1))/2.0
                        dxw = (dx(i)+dx(i-1))/2.0
                        U(i-1,j+1) = U(i-1,j)-dzn/dxw*(W(i,j)-W(i-1,j))
                    end if
                    if (IAND(FLAG(i+1,j+1),C_E) /= 0) then
                        U(i,j+1) = U(i,j)
                        W(i+1,j) = W(i,j)
                    end if
                    if (IAND(FLAG(i+1,j-1),C_E) /= 0) then
                        dzs = (dz(j)+dz(j-1))/2.0
                        dxe = (dx(i)+dx(i+1))/2.0
                        W(i+1,j-1) = W(i,j-1)-dxe/dzs*(U(i,j)-U(i,j-1))
                    end if

                CASE (C_NW)
                    !Case 3
                    if ((IAND(FLAG(i-1,j-1),C_F)/=0).AND.(FLAG(i-1,j-1)<C_E)) then
                        dzs = (dz(j)+dz(j-1))/2.0
                        U(i-1,j) = (dzs*U(i,j)+dx(i)*U(i-1,j-1))/(dx(i)+dzs)
                        !Case 4
                    ELSE
                        U(i-1,j) = U(i,j)
                    end if
                    !Case 3
                    if ((IAND(FLAG(i+1,j),C_F)/=0).AND.(FLAG(i+1,j)<C_O)) then
                        dxe = (dx(i)+dx(i+1))/2.0
                        W(i,j) = (dz(j)*W(i+1,j)+dxe*W(i,j-1))/(dxe+dz(j))
                        !Case 4
                    ELSE
                        W(i,j)   = W(i,j-1)
                    end if
                    !Shear velocities
                    if (IAND(FLAG(i-1,j+1),C_E) /= 0) then
                        U(i-1,j+1) = U(i-1,j)
                        W(i-1,j)   = W(i,j)
                    end if
                    if (IAND(FLAG(i-1,j-1),C_E) /= 0) then
                        dzs = (dz(j)+dz(j-1))/2.0
                        dxw = (dx(i)+dx(i-1))/2.0
                        W(i-1,j-1) = W(i,j-1)+dxw/dzs*(U(i-1,j)-U(i-1,j-1))
                    end if
                CASE (C_NWO)
                    !Case 5-6
                    if (((IAND(FLAG(i+1,j-1),  C_F)/=0).AND.(FLAG(i+1,j-1)  <C_E))) then
                        U(i,j) = U(i,j-1)
                    ELSE
                        U(i,j) = U(i,j)
                    end if
                    if (((IAND(FLAG(i-1,j-1),  C_F)/=0).AND.(FLAG(i-1,j-1)  <C_E))) then
                        U(i-1,j) = U(i-1,j-1)
                    ELSE
                        U(i-1,j) = U(i-1,j)
                    end if
                    !Case 4
                    W(i,j) = W(i,j-1)
                    !Shear velocities
                    if (IAND(FLAG(i-1,j-1),C_E) /= 0) then
                        dzs = (dz(j)+dz(j-1))/2.0
                        dxw = (dx(i)+dx(i-1))/2.0
                        W(i-1,j-1) = W(i,j-1)+dxw/dzs*(U(i-1,j)-U(i-1,j-1))
                    end if
                    if (IAND(FLAG(i+1,j-1),C_E) /= 0) then
                        dzs = (dz(j)+dz(j-1))/2.0
                        dxe = (dx(i)+dx(i+1))/2.0
                        W(i+1,j-1) = W(i,j-1)-dxe/dzs*(U(i,j)-U(i,j-1))
                    end if
                    if (IAND(FLAG(i-1,j+1),C_E) /= 0) then
                        W(i-1,j)   = W(i,j)
                        U(i-1,j+1) = U(i-1,j)
                    end if
                    if (IAND(FLAG(i+1,j+1),C_E) /= 0) then
                        W(i+1,j) = W(i,j)
                        U(i,j+1) = U(i,j)
                    end if
                CASE DEFAULT
                end SELECT
            end if
        end do
    end do
    !
    ! Second loop do  pressure boundary values
    !

    do j = 1,jmax
        do i = 1,imax
            if (.NOT.((IAND(FLAG(i,j),C_E)/=0).OR.(FLAG(i,j)<C_O))) then
                SELECT CASE (IAND(FLAG(i,j),C_NSWO))
                CASE (C_N)
                    P(i,j) = 2.*Fr**2/Re/dz(j)*(W(i,j)-W(i,j-1))
                CASE (C_O)
                    P(i,j) = 2.*Fr**2/Re/dx(i)*(U(i,j)-U(i-1,j))
                CASE (C_W)
                    P(i,j) = 2.*Fr**2/Re/dx(i)*(U(i,j)-U(i-1,j))
                CASE (C_NO)
                    dzs = (dz(j)+dz(j-1))/2.0
                    dxw = (dx(i)+dx(i-1))/2.0
                    P(i,j) = 1.*Fr**2/Re/2.* &
                    ( (U(i,j)+U(i-1,j)-U(i,j-1)-U(i-1,j-1))/dzs &
                    +   (W(i,j)+W(i,j-1)-W(i-1,j)-W(i-1,j-1))/dxw )
                CASE (C_NW)
                    dzs = (dz(j)+dz(j-1))/2.0
                    dxe = (dx(i)+dx(i+1))/2.0
                    P(i,j) = -1.*Fr**2/Re/2.* &
                    ( (U(i,j)+U(i-1,j)-U(i,j-1)-U(i-1,j-1))/dzs &
                    +   (W(i+1,j)+W(i+1,j-1)-W(i,j)-W(i,j-1))/dxe )
                CASE (C_WO)
                    P(i,j) = 0.
                CASE (C_NWO)
                    P(i,j) = 0.
                CASE DEFAULT
                end SELECT
            end if
        end do
    end do
    return
    end subroutine SET_UVP_SURFACE_TUMMAC

    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DIC %%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine DIC(imax, jmax, FLAG, dx, dz, U,  H, dt, St)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  advance free surface to its new position using depth integrated continuity equation
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use defs
    implicit none
    integer, INTENT(IN):: imax, jmax
    integer(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
    real(RP), DIMENSION(0:), INTENT(IN):: dx, dz

    real(RP), DIMENSION(0:), INTENT(INOUT):: H
    real(RP), DIMENSION(0:,0:), INTENT(IN):: U
    real(RP), INTENT(IN):: dt, St
    !
    !-----------------------------------------------------------------------
    !
    integer:: i, j
    real(RP) :: FR, FL
    do i=1, imax
        FR = 0.0
        FL = 0.0
        do j=1, jmax
            if (((IAND(FLAG(i-1,j),C_E) /= 0) .AND. (IAND(FLAG(i,j),C_E) /= 0)) &
            .OR. ((FLAG(i-1,j) == C_B) .AND. (IAND(FLAG(i,j),C_E) /= 0)))  then
                FL = FL
            ELSE
                FL = FL + U(i-1,j)*dz(j)
            end if
            if (((IAND(FLAG(i,j),C_E) /= 0) .AND. (IAND(FLAG(i+1,j),C_E) /= 0)))  then
                FR = FR
            ELSE
                FR = FR + U(i,j)*dz(j)
            end if
        end do
        H(i) = H(i) - St*dt/dx(i)*(FR-FL)
    end do
    H(0) = H(1)
    H(imax+1) = H(imax)
    return
    end subroutine DIC
    !!    
    !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !!
    !subroutine FCT(imax, jmax, H, cfl, dx)
    !!
    !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !!
    !!  Apply Flux Corrected Tranport method for nonuniform mesh
    !!
    !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !!
    !use nrtype
    !!use interfaces2
    !implicit none
    !integer, INTENT(IN):: imax, jmax
    !real(RP), INTENT(IN):: cfl
    !real(RP), DIMENSION(0:), INTENT(INOUT):: H
    !real(RP), DIMENSION(0:), INTENT(IN):: dx
    !!
    !!-----------------------------------------------------------------------
    !!
    !integer:: i, j
    !real(RP)::  eps1, eps2, AP1, AM1,  SP1, SM1, &
    !LP, LM, k1, k2
    !real(RP), DIMENSION(:), ALLOCATABLE :: HD
    !real(RP):: a, b, c, dxe, dxw, dxee, dxww
    !ALLOCATE(HD(0:imax+1))
    !!eps1 = cfl/(1.0+2.0*cfl)**2*1.5
    !eps1 = cfl/(1.0+2.0*cfl)**2
    !eps2 = cfl/(1.0+2.0*cfl)**2
    !do i=1, imax
    !    dxe = (dx(i)+dx(i+1))/2.0
    !    dxw = (dx(i)+dx(i-1))/2.0
    !    a = eps1*dxw/dx(i)
    !    b = dxe/dx(i)*(-eps1) + dxw/dx(i)*(-eps1)
    !    c = eps1*dxe/dx(i)
    !    HD (i) = H(i) + a*H(i-1)+b*H(i)+c*H(i+1)
    !end do
    !HD(0) = HD(1)
    !HD(imax+1) = HD(imax)
    !do i=1, imax
    !    dxe = (dx(i)+dx(i+1))/2.0
    !    dxw = (dx(i)+dx(i-1))/2.0
    !    AP1 = eps2*dxe*(HD(i+1)-HD(i))
    !    AM1 = eps2*dxw*(HD(i)-HD(i-1))
    !    SP1 = sgn(AP1)
    !    SM1 = sgn(AM1)
    !    if (i==imax) then
    !        k1 = HD(i+1)
    !        dxee = dxe
    !    ELSE
    !        k1 = HD(i+2)
    !        dxee = (dx(i+1)+dx(i+2))/2.0
    !    end if
    !    if (i==1) then
    !        k2 = HD(i-1)
    !        dxww = dxw
    !    ELSE
    !        k2 = HD(i-2)
    !        dxww = (dx(i-1)+dx(i-2))/2.0
    !    end if
    !    LP = SP1*max(0.0,min(abs(AP1), SP1*(k1-HD(i+1))*dxee, SP1*(HD(i)-HD(i-1))*dxw))
    !    LM = SM1*max(0.0,min(abs(AM1), SM1*(HD(i-1)-k2)*dxww, SM1*(HD(i+1)-HD(i))*dxe))
    !    H(i) = HD(i) -(LP-LM)/dx(i)
    !end do
    !H(0) = H(1)
    !H(imax+1) = H(imax)
    !DEALLOCATE(HD)
    !return
    !end subroutine FCT
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SGN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    real function SGN(x)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  find the sign of a given value
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    implicit none
    real(RP), INTENT(IN):: x
    !
    !-----------------------------------------------------------------------
    !
    if (x>0.0) then
        SGN = 1.0
    ELSEif (x<0.0) then
        SGN = -1.0
    ELSE
        SGN = 0.0
    end if
    return
    end function SGN
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET_UVP_SURFACE %%%%%%%%%%%%%%%%%%%%%%%%%
!
   subroutine SET_UWP_SURFACE_VARIABLE (U, W, PS, FLAG, GX, GZ, imax, jmax, &
                                        Fr, Re, dx, dz, delt, c)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Set boundary values at free surface
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype
	  USE defs
      IMPLICIT NONE
      INTEGER, intent(IN) :: imax, jmax
      real(RP), intent(IN) :: GX, GZ, Fr, Re, delt
      INTEGER(I2B), dimension(0:,0:), intent(IN) :: FLAG
      real(RP), dimension(0:,0:), intent(INOUT) :: U, W, PS
      real(RP), dimension(0:), intent(IN):: dx, dz
      real(DP), dimension(0:,0:), intent(IN) :: c
!
! local variables
!
      INTEGER :: i, j
      real :: dxe, dxw, dzn, dzs
      real(DP)::cip1, cim1, cjp1, cjm1, hi, hip1, him1, hj, hjp1, hjm1, hx, hxx, &
              curv, sten, hz, hzz, cangle, rangle, pi, yi, yip, yim, xj, xjp, xjm, &
              sl1, sl2, dhw, dhe, dhn, dhs
      pi=4.d0*datan(1.0d0)
      sten=7.3d-2     !surface tension N/m
      cangle=90.0d0       ! contact angle between the free surface and the solid boundary
!
! Set velocity values in empty cells to zero
!
      DO j = 1,jmax
         DO i = 1,imax-1
            IF ( (IAND(FLAG(i,j),  C_E) /= 0) .AND. &
                 (IAND(FLAG(i+1,j),C_E) /= 0) ) THEN
               U(i,j) = 0.0
            END IF
         END DO
      END DO

      DO  j = 1,jmax-1
         DO i = 1,imax
            PS(i,j) = 0.0
            IF ( (IAND(FLAG(i,j),  C_E) /= 0) .AND. &
                 (IAND(FLAG(i,j+1),C_E) /= 0)) THEN
               W(i,j) = 0.0
            END IF
         END DO
      END DO

      DO j = 1,jmax
         DO i = 1,imax
          !
          ! treat only surface cells
          !
          IF ( (IAND(FLAG(i,j),C_E) == 0) .OR. (FLAG(i,j) < C_O) ) THEN
           dxe = (dx(i)+dx(i+1))/2.0
           dxw = (dx(i)+dx(i-1))/2.0
           dzn = (dz(j)+dz(j+1))/2.0
           dzs = (dz(j)+dz(j-1))/2.0

           SELECT CASE (IAND(FLAG(i,j),C_NSWO)) ! masj NSWO_E=0x0f00
	                                        ! filters surface cells
             CASE (C_N)
               W(i,j) = W(i,j-1)-dz(j)/dx(i)*(U(i,j)-U(i-1,j))
               IF (IAND(FLAG(i-1,j+1),C_E) /= 0) THEN
                 U(i-1,j+1) = U(i-1,j)-dzn/dxw*(W(i,j)-W(i-1,j))
               END IF

             CASE (C_S)
               W(i,j-1) = W(i,j)+dz(j)/dx(i)*(U(i,j)-U(i-1,j))
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 U(i-1,j-1) = U(i-1,j)+dzs/dxw*(W(i,j-1)-W(i-1,j-1))
               END IF

	     CASE (C_O)
               U(i,j) = U(i-1,j)-dx(i)/dz(j)*(W(i,j)-W(i,j-1))
               IF (IAND(FLAG(i+1,j-1),C_E) /= 0) THEN
                 W(i+1,j-1) = W(i,j-1)-dxe/dzs*(U(i,j)-U(i,j-1))
               END IF

	     CASE (C_W)
               U(i-1,j) = U(i,j)+dx(i)/dz(j)*(W(i,j)-W(i,j-1))
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 W(i-1,j-1) = W(i,j-1)+dxw/dzs*(U(i-1,j)-U(i-1,j-1))
               END IF

         CASE (C_NO)
               U(i,j) = U(i-1,j)
               W(i,j) = W(i,j-1)
               IF (IAND(FLAG(i-1,j+1),C_E) /= 0) THEN
                 U(i-1,j+1) = U(i-1,j)-dzn/dxw*(W(i,j)-W(i-1,j))
               END IF
               IF (IAND(FLAG(i+1,j+1),C_E) /= 0) THEN
                 U(i,j+1) = U(i,j)
                 W(i+1,j) = W(i,j)
               END IF
               IF (IAND(FLAG(i+1,j-1),C_E) /= 0) THEN
                W(i+1,j-1) = W(i,j-1)-dxe/dzs*(U(i,j)-U(i,j-1))
               END IF

         CASE (C_NW)
               U(i-1,j) = U(i,j)
               W(i,j)   = W(i,j-1)
               IF (IAND(FLAG(i-1,j+1),C_E) /= 0) THEN
                 U(i-1,j+1) = U(i-1,j)
                 W(i-1,j)   = W(i,j)
               END IF
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 W(i-1,j-1) = W(i,j-1)+dxw/dzs*(U(i-1,j)-U(i-1,j-1))
               END IF

             CASE (C_SW)
               U(i-1,j) = U(i,j)
               W(i,j-1) = W(i,j)
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 U(i-1,j-1) = U(i-1,j)
                 W(i-1,j-1) = W(i,j-1)
               END IF

             CASE (C_SO)
               U(i,j)   = U(i-1,j)
               W(i,j-1) = W(i,j)
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 U(i-1,j-1) = U(i-1,j)+dzs/dxw*(W(i,j-1)-W(i-1,j-1))
               END IF
               IF (IAND(FLAG(i+1,j-1),C_E) /= 0) THEN
                 U(i,j-1)   = U(i,j)
                 W(i+1,j-1) = W(i,j-1)
               END IF

             CASE (C_WO)
               U(i,j)   = U(i,j)   + delt*GX
               U(i-1,j) = U(i-1,j) + delt*GX
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 W(i-1,j-1) = W(i,j-1)+dxw/dzs*(U(i-1,j)-U(i-1,j-1))
               END IF
               IF (IAND(FLAG(i+1,j-1),C_E) /= 0) THEN
                 W(i+1,j-1) = W(i,j-1)-dxe/dzs*(U(i,j)-U(i,j-1))
               END IF

             CASE (C_NS)
               W(i,j)   = W(i,j)   + delt*GZ
               W(i,j-1) = W(i,j-1) + delt*GZ
               IF (IAND(FLAG(i-1,j+1),C_E) /= 0) THEN
                 U(i-1,j+1) = U(i-1,j)-dzn/dxw*(W(i,j)-W(i-1,j))
               END IF
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 U(i-1,j-1) = U(i-1,j)+dzs/dxw*(W(i,j-1)-W(i-1,j-1))
               END IF

             CASE (C_NWO)
               W(i,j) = W(i,j-1)-dz(j)/dx(i)*(U(i,j)-U(i-1,j))
               U(i,j) = U(i,j) + delt*GX
               U(i-1,j) = U(i-1,j) + delt*GX
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 W(i-1,j-1) = W(i,j-1)+dxw/dzs*(U(i-1,j)-U(i-1,j-1))
               END IF
               IF (IAND(FLAG(i+1,j-1),C_E) /= 0) THEN
                 W(i+1,j-1) = W(i,j-1)-dxe/dzs*(U(i,j)-U(i,j-1))
               END IF
               IF (IAND(FLAG(i-1,j+1),C_E) /= 0) THEN
                 W(i-1,j)   = W(i,j)
                 U(i-1,j+1) = U(i-1,j)
               END IF
               IF (IAND(FLAG(i+1,j+1),C_E) /= 0) THEN
                 W(i+1,j) = W(i,j)
                 U(i,j+1) = U(i,j)
               END IF

             CASE (C_NSW)
               U(i-1,j) = U(i,j)+dx(i)/dz(j)*(W(i,j)-W(i,j-1))
               W(i,j)   = W(i,j) + delt*GZ
               W(i,j-1) = W(i,j-1) + delt*GZ
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 W(i-1,j-1)  = W(i,j-1)
                 U(i-1,j-1)  = U(i-1,j)
               END IF
               IF (IAND(FLAG(i-1,j+1),C_E) /= 0) THEN
                 W(i-1,j)   = W(i,j)
                 U(i-1,j+1) = U(i-1,j)
               END IF

             CASE (C_SWO)
               W(i,j-1) = W(i,j)+dz(j)/dx(i)*(U(i,j)-U(i-1,j))
               U(i,j)   = U(i,j) + delt*GX
               U(i-1,j) = U(i-1,j) + delt*GX
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 U(i-1,j-1) = U(i-1,j)
                 W(i-1,j-1) = W(i,j-1)
               END IF
               IF (IAND(FLAG(i+1,j-1),C_E) /= 0) THEN
                 U(i,j-1)    = U(i,j)
                 W(i+1,j-1)  = W(i,j-1)
               END IF

             CASE (C_NSO)
               U(i,j) = U(i-1,j)-dx(i)/dz(j)*(W(i,j)-W(i,j-1))
               W(i,j)   = W(i,j) + delt*GZ
               W(i-1,j) = W(i-1,j) + delt*GZ
               IF (IAND(FLAG(i-1,j+1),C_E) /= 0) THEN
                 U(i-1,j+1) = U(i-1,j)-dzn/dxw*(W(i,j)-W(i-1,j))
               END IF
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 U(i-1,j-1) = U(i-1,j)+dzs/dxw*(W(i,j-1)-W(i-1,j-1))
               END IF
               IF (IAND(FLAG(i+1,j-1),C_E) /= 0) THEN
                 U(i,j-1)   = U(i,j)
                 W(i+1,j-1) = W(i,j-1)
               END IF
               IF (IAND(FLAG(i+1,j+1),C_E) /= 0) THEN
                 U(i,j+1)    = U(i,j)
                 W(i+1,j)    = W(i,j)
               END IF

             CASE (C_NSWO)
               U(i,j)   = U(i,j)   + delt*GX
               U(i-1,j) = U(i-1,j) + delt*GX
               W(i,j)   = W(i,j)   + delt*GZ
               W(i,j-1) = W(i,j-1) + delt*GZ
               IF (IAND(FLAG(i-1,j+1),C_E) /= 0) THEN
                 U(i-1,j+1) = U(i-1,j)
                 W(i-1,j)   = W(i,j)
               END IF
               IF (IAND(FLAG(i+1,j+1),C_E) /= 0) THEN
                 U(i,j+1) = U(i,j)
                 W(i+1,j) = W(i,j)
               END IF
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 U(i-1,j-1)  = U(i-1,j)
                 W(i-1,j-1)  = W(i,j-1)
               END IF
               IF (IAND(FLAG(i+1,j-1),C_E) /= 0) THEN
                 U(i,j-1)    = U(i,j)
                 W(i+1,j-1)  = W(i,j-1)
               END IF

	     CASE DEFAULT

           END SELECT

          END IF

         END DO
      END DO
!
! Second loop DO  pressure boundary values
!
    do j = 1,jmax
        do i = 1,imax
            dxe = (dx(i)+dx(i+1))/2.0
            dxw = (dx(i)+dx(i-1))/2.0
            dzn = (dz(j)+dz(j+1))/2.0
            dzs = (dz(j)+dz(j-1))/2.0
            if (curvature) then
                if (cangle==90.0d0) cangle=cangle-emf
                rangle=cangle*pi/180.0d0
                if (.NOT.((IAND(FLAG(i,j),C_E)/=0).OR.(FLAG(i,j)<C_O))) then
                    yi  = c(i,j-1)*dz(j-1)+c(i,j)*dz(j)+c(i,j+1)*dz(j+1)
                    yip = c(i+1,j-1)*dz(j-1)+c(i+1,j)*dz(j)+c(i+1,j+1)*dz(j+1)
                    yim = c(i-1,j-1)*dz(j-1)+c(i-1,j)*dz(j)+c(i-1,j+1)*dz(j+1)
                    xj =  c(i-1,j)*dx(i-1)+c(i,j)*dx(i)+c(i+1,j)*dx(i+1)
                    xjp = c(i-1,j+1)*dx(i-1)+c(i,j+1)*dx(i)+c(i+1,j+1)*dx(i+1)
                    xjm = c(i-1,j-1)*dx(i-1)+c(i,j-1)*dx(i)+c(i+1,j-1)*dx(i+1)
                    !DY/DX
                    sl1 = 2.0*(yip-yim)/(dx(i+1)+2.0*dx(i)+dx(i-1))
                    !DX/DY
                    sl2 = 2.0*(xjp-xjm)/(dz(j+1)+2.0*dz(j)+dz(j-1))
                    if (abs(sl1)<abs(sl2)) then           !free surface is horizontal, a vertical stencil is used
                                !contact line cell: along the solid boundary, if the interface cell has a neighbor along
                                !the solid boundary that is completely empty cell
                        !left cell is B cell
                        if (iand(FLAG(i-1,j),C_F)/=C_F) then
                            yim=yi+dxw/dtan(rangle)
                        end if
                        !right cell is B cell
                        if (iand(FLAG(i+1,j),C_F) /= C_F) then
                            yip=yi+dxe/dtan(rangle)
                        end if
                        !bottom cell is B cell
                        if (iand(FLAG(i,j-1),C_F) /= C_F) then
                            !if the left cell is an empty cell
                            if (IAND(FLAG(i-1,j),C_E)/=0) then
                                yim=yi-dxw*dtan(rangle)
                            end if
                            !if the right cell is an empty cell
                            if (IAND(FLAG(i+1,j),C_E)/=0) then
                                yip=yi-dxe*dtan(rangle)
                            end if
                        end if
                        if (iand(FLAG(i,j+1),C_F) /= C_F) then
                            !if the left cell is an empty cell
                            if (IAND(FLAG(i-1,j),C_E)/=0) then
                                yim=yi-dxw*dtan(rangle)
                            end if
                            !if the right cell is an empty cell
                            if (IAND(FLAG(i+1,j),C_E)/=0) then
                                yip=yi-dxe*dtan(rangle)
                            end if
                        end if
                        dhe=(yip-yi)/dxe
                        dhw=(yi-yim)/dxw
                        curv=(dhe/dsqrt(1.0d0+dhe**2)-dhw/dsqrt(1.0d0+dhw**2))/dx(i)
                     else       !free surface is vertical, a horizontal stencil is used
                                !contact line cell: along the solid boundary, if the interface cell has a neighbor along
                                !the solid boundary that is completely empty cell
                        !left cell is B cell
                        if (iand(FLAG(i-1,j),C_F)/=C_F) then
                            !if the upper cell is an empty cell
                            if (AND(FLAG(i,j+1),C_E)/=0) then
                                xjp=xj-dzn*dtan(rangle)
                            end if
                            !if the below cell is an empty cell
                            if (IAND(FLAG(i,j-1),C_E)/=0) then
                                xjm=xj-dzs*dtan(rangle)
                            end if
                        end if
                        !right cell is B cell
                        if (iand(FLAG(i+1,j),C_F) /= C_F) then
                            !if the upper cell is an empty cell
                            if (IAND(FLAG(i,j+1),C_E)/=0) then
                                xjp=xj-dzn*dtan(rangle)
                            end if
                            !if the below cell is an empty cell
                            if (IAND(FLAG(i,j-1),C_E)/=0) then
                                xjm=xj-dzs*dtan(rangle)
                            end if
                        end if
                        !bottom cell is B cell
                        if (iand(FLAG(i,j-1),C_F) /= C_F) then
                            xjm=xj+dzs/dtan(rangle)
                        end if
                        if (iand(FLAG(i,j+1),C_F) /= C_F) then
                            xjp=xj+dzn/dtan(rangle)
                        end if
                        dhn=(xjp-xj)/dzn
                        dhs=(xj-xjm)/dzs
                        curv=(dhn/dsqrt(1.0d0+dhn**2)-dhs/dsqrt(1.0d0+dhs**2))/dz(j)
                    end if
                    PS(i,j)=-curv*sten
                end if
            else
                if (.NOT.((IAND(FLAG(i,j),C_E)/=0).OR.(FLAG(i,j)<C_O))) then
                    select case (iand(FLAG(i,j),C_NSWO))
                    case (C_N)
                        PS(i,j) = 2.*Fr**2/Re/dz(j)*(W(i,j)-W(i,j-1))
                    case (C_S)
                        PS(i,j) = 2.*Fr**2/Re/dz(j)*(W(i,j)-W(i,j-1))
                    case (C_O)
                        PS(i,j) = 2.*Fr**2/Re/dx(i)*(U(i,j)-U(i-1,j))
                    case (C_W)
                        PS(i,j) = 2.*Fr**2/Re/dx(i)*(U(i,j)-U(i-1,j))
                    case (C_NO)
                        PS(i,j) = 1.*Fr**2/Re/2.* &
                        ( (U(i,j)+U(i-1,j)-U(i,j-1)-U(i-1,j-1))/dzs &
                        +   (W(i,j)+W(i,j-1)-W(i-1,j)-W(i-1,j-1))/dxw )
                    case (C_NW)
                        PS(i,j) = -1.*Fr**2/Re/2.* &
                        ( (U(i,j)+U(i-1,j)-U(i,j-1)-U(i-1,j-1))/dzs &
                        +   (W(i+1,j)+W(i+1,j-1)-W(i,j)-W(i,j-1))/dxe )
                    case (C_SW)
                        PS(i,j) = 1.*Fr**2/Re/2.* &
                        ( (U(i,j+1)+U(i-1,j+1)-U(i,j)-U(i-1,j))/dzn &
                        +   (W(i+1,j)+W(i+1,j-1)-W(i,j)-W(i,j-1))/dxe)
                    case (C_SO)
                        PS(i,j) = -1.*Fr**2/Re/2.* &
                        ( (U(i,j+1)+U(i-1,j+1)-U(i,j)-U(i-1,j))/dzn &
                        +   (W(i,j)+W(i,j-1)-W(i-1,j)-W(i-1,j-1))/dxw)
                    case (C_WO)
                        PS(i,j) = 0.
                    case (C_NS)
                        PS(i,j) = 0.
                    case (C_NWO)
                        PS(i,j) = 0.
                    case (C_NSW)
                        PS(i,j) = 0.
                    case (C_SWO)
                        PS(i,j) = 0.
                    case (C_NSO)
                        PS(i,j) = 0.
                    case (C_NSWO)
                        PS(i,j) = 0.
                    case default
                    end select
                end if
            end if
        end do
    end do

      RETURN

END SUBROUTINE SET_UWP_SURFACE_VARIABLE
