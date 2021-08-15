!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	subroutine swpx (c, U, W, imax, jmax, wW, wE, delt, dx, dz)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  computes color field
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      use nrtype
	  use defs
      implicit none

      integer, intent(IN) :: imax, jmax, wW, wE
      real(RP), intent(IN) :: delt
      real(RP), dimension(0:), intent (IN):: dx, dz
      real(RP), dimension(0:,0:), intent(IN) :: U, W
      real(DP), dimension(0:,0:), intent(INOUT) :: c
!
!  local variables
!
      integer :: i, j, id
      real(DP), dimension(0:imax+1,0:jmax+1) :: nx, nz, s, lflux, rflux
      real(DP) cc, u1, u2, mrx, mrz, mx, mz
	  real(DP) are1, are2, du1, du2
	  real(DP) dzn, dzs, dxe, dxw, mx1, mx2, mx3, mx4, mz1, mz2, mz3, mz4, &
               len, dx1, dz1, dxc, dzc, sc, Fc, sm, ss, nxx, nzz, delx, delz, dvf0, dv, lf, rf, mm1, mm2

	  intrinsic dmax1, dmin1, dabs, dsqrt
!
!     Solve: dc/dt = -div(uc) + c div(u) (c changes only at cells
!     where grad(c) is not zero)
!
      do j=1, jmax
          dzn = (dz(j)+dz(j+1))/2.0
          dzs = (dz(j)+dz(j-1))/2.0
          do i=1, imax
             
              if (c(i,j) .ne. 0.0d0 .and. c(i,j) .ne. 1.0d0) then
                  dxe = (dx(i)+dx(i+1))/2.0
                  dxw = (dx(i)+dx(i-1))/2.0
                  mx1 = ((c(i+1,j+1)-c(i,j+1))*dz(j) + &
                  (c(i+1,j)-c(i,j))*dz(j+1))/(dz(j)+dz(j+1))/dxe
                  mx2 = ((c(i,j+1)-c(i-1,j+1))*dz(j) + &
                  (c(i,j)-c(i-1,j))*dz(j+1))/(dz(j)+dz(j+1))/dxw
                  mx3 = ((c(i+1,j)-c(i,j))*dz(j-1) + &
                  (c(i+1,j-1)-c(i,j-1))*dz(j))/(dz(j)+dz(j-1))/dxe
                  mx4 = ((c(i,j)-c(i-1,j))*dz(j-1) + &
                  (c(i,j-1)-c(i-1,j-1))*dz(j))/(dz(j)+dz(j-1))/dxw
                  mz1 = ((c(i+1,j+1)-c(i+1,j))*dx(i) + &
                  (c(i,j+1)-c(i,j))*dx(i+1))/(dx(i)+dx(i+1))/dzn
                  mz2 = ((c(i,j+1)-c(i,j))*dx(i-1) + &
                  (c(i-1,j+1)-c(i-1,j))*dx(i))/(dx(i)+dx(i-1))/dzn
                  mz3 = ((c(i+1,j)-c(i+1,j-1))*dx(i) + &
                  (c(i,j)-c(i,j-1))*dx(i+1))/(dx(i)+dx(i+1))/dzs
                  mz4 = ((c(i,j)-c(i,j-1))*dx(i-1) + &
                  (c(i-1,j)-c(i-1,j-1))*dx(i))/(dx(i)+dx(i-1))/dzs

                  mx = 1.0d0/4.0d0*(mx1+mx2+mx3+mx4)
                  mz = 1.0d0/4.0d0*(mz1+mz2+mz3+mz4)

                  nx(i,j) = mx
                  nz(i,j) = mz

                  mrx = dabs(nx(i,j))
                  mrz = dabs(nz(i,j))
                  dx1=mrx*dx(i)+1.d-50
                  dz1=mrz*dz(j)+1.d-50
                  sm=dx1+dz1
                  dxc=dmax1(dx1,dz1)
                  dzc=dmin1(dx1,dz1)
                  are1 = 0.5d0*dzc/dxc
                  Fc=dmin1(c(i,j),1.0d0-c(i,j))
                  if (Fc .lt. are1) then
                      sc = dsqrt(2.0d0*Fc*dxc*dzc)
                  else
                      sc = Fc*dxc + 0.5d0*dzc
                  end if
                  if (c(i,j) .le. 0.5d0) then 
                      s(i,j)=sc
                  else
                      s(i,j)=sm-sc
                  end if
              end if
          end do
      end do

      do j=1, jmax
          delz=dz(j)
          do i=1, imax
              u1=U(i-1,j)*delt
              u2=U(i,j)*delt
              !compute flux on the left face
              if (U(i-1,j) .gt. 0.0) then
                  id=i-1
              else
                  id=i
              end if

              if ((c(id,j) .eq. 1.0d0) .or. (c(id,j) .eq. 0.0d0)) then 
                  lflux(i,j)=c(id,j)
              else        
                  ss=s(id,j)
                  nxx=nx(id,j)
                  nzz=nz(id,j)
                  mrx = dabs(nxx)+1.d-50
                  mrz = dabs(nzz)+1.d-50
                  delx=dabs(u1)+1.d-50
                  dv=delx*delz
                  if (nx(id,j)*u1 .lt. 0.0d0) delx=dx(id)-dabs(u1)
                  dxc=dmax1(mrx*delx,mrz*delz)
                  dzc=dmin1(mrx*delx,mrz*delz)
                  sm=dxc+dzc
                  sc=dmin1(ss,dmax1(sm-ss,0.0d0))
                  if (sc .lt. dzc) then 
                      Fc=0.5d0*sc**2/dxc/dzc
                  else
                      Fc=(sc-0.5d0*dzc)/dxc
                  end if
                  if (ss .le. 0.5d0*sm) then 
                      dvf0=delx*delz*Fc
                  else
                      dvf0=delx*delz*(1.0d0-Fc)
                  end if
                  if (nx(id,j)*u1 .lt. 0.0d0) then 
                      dvf0=c(id,j)*dx(id)*dz(j)-dvf0
                  end if
                  lflux(i,j)=dvf0/dv
              end if
              !compute flux on the right face
              if (U(i,j) .gt. 0.0) then 
                  id=i
              else
                  id=i+1
              end if
              !if the donor cell is single-phase cell where c=0 or 1, we simply obtain dvf0=F*dv
              if ((c(id,j) .eq. 1.0d0) .or. (c(id,j) .eq. 0.0d0)) then 
                  rflux(i,j)=c(id,j)
              else
                  ss=s(id,j)
                  nxx=nx(id,j)
                  nzz=nz(id,j)
                  mrx = dabs(nxx)+1.d-50
                  mrz = dabs(nzz)+1.d-50
                  delx=dabs(u2)+1.d-50
                  dv=delx*delz
                  if (nx(id,j)*u2 .lt. 0.0d0) delx=dx(id)-dabs(u2)
                  dxc=dmax1(mrx*delx,mrz*delz)
                  dzc=dmin1(mrx*delx,mrz*delz)
                  sm=dxc+dzc
                  sc=dmin1(ss,dmax1(sm-ss,0.0d0))
                  if (sc .lt. dzc) then 
                      Fc=0.5d0*sc**2/dxc/dzc
                  else
                      Fc=(sc-0.5d0*dzc)/dxc
                  end if
                  if (ss .le. 0.5d0*sm) then 
                      dvf0=delx*delz*Fc
                  else
                      dvf0=delx*delz*(1.0d0-Fc)
                  end if
                  if (nx(id,j)*u2 .lt. 0.0d0) then 
                      dvf0=c(id,j)*dx(id)*dz(j)-dvf0
                  end if
                  rflux(i,j)=dvf0/dv
              end if
          end do
      end do
!     BC on the fluxes
      do j=1, jmax
          if (wW==3.or.wW==4) then 
              lflux(1,j)=lflux(2,j)
          else
              !lflux(i,j)=0.0d0
              lflux(1,j)=U(0,j)*dz(j)*c(1,j)
          end if
          
          if (wE==3.or.wE==4) then 
              rflux(imax,j)=rflux(imax-1,j)
          else
              rflux(imax,j)=0.0d0
          end if
          
      end do
     
      do i=1, imax
          do j=1, jmax
              u1=U(i-1,j)*delt/dx(i)
              u2=U(i,j)*delt/dx(i)
              c(i,j)=(c(i,j)+u1*lflux(i,j)-u2*rflux(i,j))/(1.0d0+u1-u2)
              c(i,j)=dmin1(1.0d0,dmax1(0.0d0,c(i,j)))
          end do
      end do

      return
    end subroutine swpx
    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	subroutine swpz (c, U, W, imax, jmax, wE, delt, dx, dz)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  computes color field
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      use nrtype
	  use defs
      implicit none

      integer, intent(IN) :: imax, jmax, wE
      real(RP), intent(IN) :: delt
      real(RP), dimension(0:), intent (IN):: dx, dz
      real(RP), dimension(0:,0:), intent(IN) :: U, W
      real(DP), dimension(0:,0:), intent(INOUT) :: c
!
!  local variables
!
      integer :: i, j, jd
      real(DP), dimension(0:imax+1,0:jmax+1) :: nx, nz, s, bflux, tflux
      real(DP) cc, w1, w2, mrx, mrz, mx, mz
	  real(DP) are1, are2, dw1, dw2
	  real(DP) dzn, dzs, dxe, dxw, mx1, mx2, mx3, mx4, mz1, mz2, mz3, mz4, &
               len, dx1, dz1, dxc, dzc, sc, Fc, sm, ss, nxx, nzz, delx, delz, dvf0, dv, bf, tf

	  intrinsic dmax1, dmin1, dabs, dsqrt
!
!     Solve: dc/dt = -div(uc) + c div(u) (c changes only at cells
!     where grad(c) is not zero)
!
      do j=1, jmax
          dzn = (dz(j)+dz(j+1))/2.0
          dzs = (dz(j)+dz(j-1))/2.0
          do i=1,imax
              if (c(i,j) .ne. 0.0d0 .and. c(i,j) .ne. 1.0d0) then
                  dxe = (dx(i)+dx(i+1))/2.0
                  dxw = (dx(i)+dx(i-1))/2.0
                  mx1 = ((c(i+1,j+1)-c(i,j+1))*dz(j) + &
                  (c(i+1,j)-c(i,j))*dz(j+1))/(dz(j)+dz(j+1))/dxe
                  mx2 = ((c(i,j+1)-c(i-1,j+1))*dz(j) + &
                  (c(i,j)-c(i-1,j))*dz(j+1))/(dz(j)+dz(j+1))/dxw
                  mx3 = ((c(i+1,j)-c(i,j))*dz(j-1) + &
                  (c(i+1,j-1)-c(i,j-1))*dz(j))/(dz(j)+dz(j-1))/dxe
                  mx4 = ((c(i,j)-c(i-1,j))*dz(j-1) + &
                  (c(i,j-1)-c(i-1,j-1))*dz(j))/(dz(j)+dz(j-1))/dxw
                  mz1 = ((c(i+1,j+1)-c(i+1,j))*dx(i) + &
                  (c(i,j+1)-c(i,j))*dx(i+1))/(dx(i)+dx(i+1))/dzn
                  mz2 = ((c(i,j+1)-c(i,j))*dx(i-1) + &
                  (c(i-1,j+1)-c(i-1,j))*dx(i))/(dx(i)+dx(i-1))/dzn
                  mz3 = ((c(i+1,j)-c(i+1,j-1))*dx(i) + &
                  (c(i,j)-c(i,j-1))*dx(i+1))/(dx(i)+dx(i+1))/dzs
                  mz4 = ((c(i,j)-c(i,j-1))*dx(i-1) + &
                     (c(i-1,j)-c(i-1,j-1))*dx(i))/(dx(i)+dx(i-1))/dzs
                  mx = 1.0d0/4.0d0*(mx1+mx2+mx3+mx4)
                  mz = 1.0d0/4.0d0*(mz1+mz2+mz3+mz4)

                  nx(i,j) = mx
                  nz(i,j) = mz

                  mrx = dabs(nx(i,j))
                  mrz = dabs(nz(i,j))
                  dx1=mrx*dx(i)+1.d-50
                  dz1=mrz*dz(j)+1.d-50
                  sm=dx1+dz1
                  dxc=dmax1(dx1,dz1)
                  dzc=dmin1(dx1,dz1)
                  are1 = 0.5d0*dzc/dxc
                  Fc=dmin1(c(i,j),1.0d0-c(i,j))
                  if (Fc .lt. are1) then
                      sc = dsqrt(2.0d0*Fc*dxc*dzc)
                  else
                      sc = Fc*dxc + 0.5d0*dzc
                  endif
                  if (c(i,j) .le. 0.5d0) then 
                      s(i,j)=sc
                  else
                      s(i,j)=sm-sc
                  end if
              end if
          end do
      end do
      
      do j=1, jmax
          do i=1, imax
              w1=W(i,j-1)*delt
              w2=W(i,j)*delt
              delx=dx(i)
              !compute flux on bottom face
              if (W(i,j-1) .gt. 0.0) then
                  jd=j-1
              else
                  jd=j
              end if
              if ((c(i,jd) .eq. 1.0d0) .or. (c(i,jd) .eq. 0.0d0)) then 
                  bflux(i,j)=c(i,jd)
              else
                  ss=s(i,jd)
                  nxx=nx(i,jd)
                  nzz=nz(i,jd)
                  mrx = dabs(nxx)+1.d-50
                  mrz = dabs(nzz)+1.d-50
                  delz=dabs(w1)+1.d-50
                  dv=delx*delz
                  if (nz(i,jd)*w1 .lt. 0.0d0)  delz=dz(jd)-dabs(w1)
                  dxc=dmax1(mrx*delx,mrz*delz)
                  dzc=dmin1(mrx*delx,mrz*delz)
                  sm=dxc+dzc
                  sc=dmin1(ss,dmax1(sm-ss,0.0d0))
                  if (sc .lt. dzc) then 
                      Fc=0.5d0*sc**2/dxc/dzc
                  else
                      Fc=(sc-0.5d0*dzc)/dxc
                  end if
                  if (ss .le. 0.5d0*sm) then 
                      dvf0=delx*delz*Fc
                  else
                      dvf0=delx*delz*(1.0d0-Fc)
                  end if
                  if (nz(i,jd)*w1 .lt. 0.0d0) then 
                      dvf0=c(i,jd)*dx(i)*dz(jd)-dvf0
                  end if
                  bflux(i,j)=dvf0/dv
              end if
              !compute flux on top face
              if (W(i,j) .gt. 0.0) then 
                  jd=j
              else
                  jd=j+1
              end if
              !if the donor cell is single-phase cell where c=0 or 1, we simply obtain dvf0=F*dv
              if ((c(i,jd) .eq. 1.0d0) .or. (c(i,jd) .eq. 0.0d0)) then 
                  tflux(i,j)=c(i,jd)
              else
                  ss=s(i,jd)
                  nxx=nx(i,jd)
                  nzz=nz(i,jd)
                  mrx = dabs(nxx)+1.d-50
                  mrz = dabs(nzz)+1.d-50
                  delz=dabs(w2)+1.d-50
                  dv=delx*delz
                  if (nz(i,jd)*w2 .lt. 0.0d0) delz=dz(jd)-dabs(w2) 
                  dxc=dmax1(mrx*delx,mrz*delz)
                  dzc=dmin1(mrx*delx,mrz*delz)
                  sm=dxc+dzc
                  sc=dmin1(ss,dmax1(sm-ss,0.0d0))
                  if (sc .lt. dzc) then 
                      Fc=0.5d0*sc**2/dxc/dzc
                  else
                      Fc=(sc-0.5d0*dzc)/dxc
                  end if
                  if (ss .le. 0.5d0*sm) then 
                      dvf0=delx*delz*Fc
                  else
                      dvf0=delx*delz*(1.0d0-Fc)
                  end if
                  if (nz(i,jd)*w2 .lt. 0.0d0) then 
                      dvf0=c(i,jd)*dx(i)*dz(jd)-dvf0
                  end if
                  tflux(i,j)=dvf0/dv
              end if
          end do
      end do
!     BC on the fluxes
      do i=1, imax
          bflux(i,1)=0.0d0
          tflux(i,jmax)=0.0d0
      end do
      
      do i=1, imax
          do j=1, jmax
              w1=W(i,j-1)*delt/dz(j)
              w2=W(i,j)*delt/dz(j)
              c(i,j)=c(i,j)*(1.0d0-w1+w2)+(w1*bflux(i,j)-w2*tflux(i,j))
              c(i,j)=dmin1(1.0d0,dmax1(0.0d0,c(i,j)))
          end do
      end do
      
      
      return
    end subroutine swpz
    
    subroutine VOF_PLIC (c, U, W, imax, jmax, wW, wE, delt, dx, dz, FLAG)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  computes color field
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      use nrtype
	  use defs
      use interfaces2

      implicit none

      integer, intent(IN) :: imax, jmax, wW, wE
      real(RP), intent(IN) :: delt
      real(RP), dimension(0:), intent (IN):: dx, dz
      real(RP), dimension(0:,0:), intent(IN) :: U, W
      real(DP), dimension(0:,0:), intent(INOUT) :: c
      integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
      
      call swpx (c, U, W, imax, jmax, wW, wE, delt, dx, dz)
      call VOFBOUND (c, FLAG, dx, dz, imax, jmax)
      call swpz (c, U, W, imax, jmax, wE, delt, dx, dz)
      call VOFBOUND (c, FLAG, dx, dz, imax, jmax)
      
      return 
    end subroutine VOF_PLIC
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VOFBOUND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine VOFBOUND (c, FLAG, dx, dz, imax, jmax)
      use nrtype
	  use defs
      implicit none

      integer, intent(IN) :: imax, jmax
      real(RP), dimension(0:), intent (IN):: dx, dz
      real(DP), dimension(0:,0:), intent(INOUT) :: c
      integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
      
      integer:: i, j
      real(DP) dzn, dzs, dxe, dxw
      
      do j=0, jmax+1
          c(0,j)=c(1,j)
          c(imax+1,j) = c(imax,j)
      enddo
      
      do i=0, imax+1
          c(i,0)=c(i,1)
          c(i,jmax+1)=c(i,jmax)
      end do
      
      do i = 1,imax
          do j = 1,jmax
              if ( IAND(FLAG(i,j),C_F) /= C_F) then 
                  c(i,j)=0.0d0 
              end if
          end do
      end do
      !
      do i=1, imax
          do j=1, jmax
              if ( (FLAG(i,j) >= B_N) .AND. (FLAG(i,j) <= B_SWE) ) then
                  dxe=(dx(i)+dx(i+1))/2.0
                  dxw=(dx(i)+dx(i-1))/2.0
                  dzn=(dz(j)+dz(j+1))/2.0
                  dzs=(dz(j)+dz(j-1))/2.0
                  if ( FLAG(i,j) == B_N ) then
                      c(i,j) = c(i,j+1)
                  elseif ( FLAG(i,j) == B_E ) then
                      c(i,j) = c(i+1,j)
                  elseif ( FLAG(i,j) == B_S ) then
                      c(i,j) = c(i,j-1)
                  elseif ( FLAG(i,j) == B_W ) then
                      c(i,j) = c(i-1,j)
                  elseif ( FLAG(i,j) == B_NE ) then
                      c(i,j)=(dxe*c(i,j+1)+dzn*c(i+1,j))/(dxe+dzn)
                  elseif ( FLAG(i,j) == B_SE ) then
                      c(i,j)=(c(i,j-1)*dxe+c(i+1,j)*dzs)/(dxe+dzs)
                  elseif ( FLAG(i,j) == B_SW ) then
                      c(i,j)=(c(i,j-1)*dxw+c(i-1,j)*dzs)/(dxw+dzs)
                  elseif ( FLAG(i,j) == B_NW ) then
                      c(i,j)=(c(i,j+1)*dxw+c(i-1,j)*dzn)/(dxw+dzn)
                  end if
              end if
          end do
      end do
      
      return 
      end subroutine VOFBOUND
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET_UVP_SURFACE_VOF %%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine SET_UVP_SURFACE_VOF (U, W, P, c, FLAG, flagib, imax, jmax, Fr, Re, GX, GZ, delt, dx, dz)
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
    real(RP), INTENT(IN) :: Fr, Re, GX, GZ, delt
    integer(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG, flagib
    real(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, W, P
    real(DP), dimension(0:,0:), intent(INOUT) :: c
    real(RP), DIMENSION(0:), INTENT(IN) :: dx, dz
    !
    !-----------------------------------------------------------------------
    !
    integer :: i, j, k, iminn, imaxx, jminn, jmaxx
    real(RP):: dzn, dzs, dxw, dxe
    real(DP)::cip1, cim1, cjp1, cjm1, hi, hip1, him1, hj, hjp1, hjm1, hx, hxx, &
              curv, sten, hz, hzz, cangle, rangle, pi, yi, yip, yim, xj, xjp, xjm, &
              sl1, sl2, dhw, dhe, dhn, dhs
    pi=4.d0*datan(1.0d0)
    sten=7.3d-2     !surface tension N/m
    cangle=135.0d0       ! contact angle between the free surface and the solid boundary
    !
    ! Set velocity values in empty cells to zero
    !
    do j = 1,jmax
        do i = 1,imax-1
            if ( (IAND(FLAG(i,j),  C_E) /= 0) .AND. (IAND(FLAG(i+1,j),C_E) /= 0) ) then
                U(i,j) = 0.0
            end if
        end do
    end do

    do  j = 1,jmax-1
        do i = 1,imax
            if ( (IAND(FLAG(i,j),  C_E) /= 0) .AND. (IAND(FLAG(i,j+1),C_E) /= 0)) then
                W(i,j) = 0.0
            end if
        end do
    end do

    !
    !Free surface boundary conditions for u_velocity
    !
    do j = 1,jmax
        do i = 1,imax
            dxw = (dx(i-1)+dx(i))/2.0
            dxe = (dx(i)+dx(i+1))/2.0   
            dzn = (dz(j)+dz(j+1))/2.0   
            dzs = (dz(j)+dz(j-1))/2.0
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
                        W(i,j) = (W(i-1,j)*dxe+W(i+1,j)*dxw)/(dxe+dxw)
                        !Case 2
                    ELSE if (((IAND(FLAG(i-1,j),C_F)/=0).AND.(FLAG(i-1,j)<C_O)) .AND. &
                    (IAND(FLAG(i+1,j+1),C_E) /= 0)) then
                        W(i,j) = (dz(j)*W(i-1,j)+dxw*W(i,j-1))/(dxw+dz(j))
                        !Case 3
                    ELSE if ((IAND(FLAG(i-1,j+1),C_E) /= 0) .AND. &
                    ((IAND(FLAG(i+1,j),C_F)/=0).AND.(FLAG(i+1,j)<C_O))) then
                        W(i,j) = (dz(j)*W(i+1,j)+dxe*W(i,j-1))/(dxe+dz(j))
                        !Case 4
                    ELSE
                        W(i,j) = W(i,j-1)
                    end if
                    !Shear velocities
                    if (IAND(FLAG(i-1,j+1),C_E) /= 0) then
                        U(i-1,j+1) = U(i-1,j)-dzn/dxw*(W(i,j)-W(i-1,j))
                    end if
                CASE (C_S)
                    W(i,j-1) = W(i,j)+dz(j)/dx(i)*(U(i,j)-U(i-1,j))
                    if (IAND(FLAG(i-1,j-1),C_E) /= 0) then
                        U(i-1,j-1) = U(i-1,j)+dzs/dxw*(W(i,j-1)-W(i-1,j-1))
                    end if
                CASE (C_O)
                    !Case 1
                    U(i,j) = U(i-1,j)
                    !Shear velocities
                    if (IAND(FLAG(i+1,j-1),C_E) /= 0) then
                        W(i+1,j-1) = W(i,j-1)-dxe/dzs*(U(i,j)-U(i,j-1))
                    end if
                CASE (C_W)
                    U(i-1,j) = U(i,j)
                    !Shear velocities
                    if (IAND(FLAG(i-1,j-1),C_E) /= 0) then
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
                        W(i-1,j-1) = W(i,j-1)+dxw/dzs*(U(i-1,j)-U(i-1,j-1))
                    end if
                    if (IAND(FLAG(i+1,j-1),C_E) /= 0) then
                        W(i+1,j-1) = W(i,j-1)-dxe/dzs*(U(i,j)-U(i,j-1))
                    end if
                CASE (C_NO)
                    !Case 1
                    if ((IAND(FLAG(i+1,j-1),C_F)/=0).AND.(FLAG(i+1,j-1)<C_E)) then
                        U(i,j) = (dzs*U(i-1,j)+dx(i)*U(i,j-1))/(dx(i)+dzs)
                        !Case 2
                    ELSE
                        U(i,j) = U(i-1,j)
                    end if
                    !Case 2
                    if ((IAND(FLAG(i-1,j),C_F)/=0).AND.(FLAG(i-1,j)<C_O)) then
                        W(i,j) = (dz(j)*W(i-1,j)+dxw*W(i,j-1))/(dxw+dz(j))
                        !Case 4
                    ELSE
                        W(i,j) = W(i,j-1)
                    end if
                    !Shear velocities
                    if (IAND(FLAG(i-1,j+1),C_E) /= 0) then
                        U(i-1,j+1) = U(i-1,j)-dzn/dxw*(W(i,j)-W(i-1,j))
                    end if
                    if (IAND(FLAG(i+1,j+1),C_E) /= 0) then
                        U(i,j+1) = U(i,j)
                        W(i+1,j) = W(i,j)
                    end if
                    if (IAND(FLAG(i+1,j-1),C_E) /= 0) then
                        W(i+1,j-1) = W(i,j-1)-dxe/dzs*(U(i,j)-U(i,j-1))
                    end if
                CASE (C_NW)
                    !Case 3
                    if ((IAND(FLAG(i-1,j-1),C_F)/=0).AND.(FLAG(i-1,j-1)<C_E)) then
                        U(i-1,j) = (dzs*U(i,j)+dx(i)*U(i-1,j-1))/(dx(i)+dzs)
                        !Case 4
                    ELSE
                        U(i-1,j) = U(i,j)
                    end if
                    !Case 3
                    if ((IAND(FLAG(i+1,j),C_F)/=0).AND.(FLAG(i+1,j)<C_O)) then
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
                        W(i-1,j-1) = W(i,j-1)+dxw/dzs*(U(i-1,j)-U(i-1,j-1))
                    end if
                CASE (C_SW)
                    U(i-1,j) = U(i,j)
                    W(i,j-1) = W(i,j)
                    if (IAND(FLAG(i-1,j-1),C_E) /= 0) then
                        U(i-1,j-1) = U(i-1,j)
                        W(i-1,j-1) = W(i,j-1)
                    end if
                CASE (C_SO)
                    U(i,j)   = U(i-1,j)
                    W(i,j-1) = W(i,j)
                    if (IAND(FLAG(i-1,j-1),C_E) /= 0) then
                        U(i-1,j-1) = U(i-1,j)+dzs/dxw*(W(i,j-1)-W(i-1,j-1))
                    end if
                    if (IAND(FLAG(i+1,j-1),C_E) /= 0) then
                        U(i,j-1)   = U(i,j)
                        W(i+1,j-1) = W(i,j-1)
                    end if
                CASE (C_NS)
                    W(i,j)   = W(i,j)   
                    W(i,j-1) = W(i,j-1) 
                    if (IAND(FLAG(i-1,j+1),C_E) /= 0) then
                        U(i-1,j+1) = U(i-1,j)-dzn/dxw*(W(i,j)-W(i-1,j))
                    end if
                    if (IAND(FLAG(i-1,j-1),C_E) /= 0) then
                        U(i-1,j-1) = U(i-1,j)+dzs/dxw*(W(i,j-1)-W(i-1,j-1))
                    end if
                CASE (C_NWO)
                    !Case 5-6
                    if (((IAND(FLAG(i+1,j-1),  C_F)/=0).AND.(FLAG(i+1,j-1)  <C_E))) then
                        U(i,j) = U(i,j-1)
                    else
                        U(i,j) = U(i,j)
                    end if
                    if (((IAND(FLAG(i-1,j-1),  C_F)/=0).AND.(FLAG(i-1,j-1)  <C_E))) then
                        U(i-1,j) = U(i-1,j-1)
                    else
                        U(i-1,j) = U(i-1,j)
                    end if
                    !Case 4
                    W(i,j) = W(i,j-1)
                    !Shear velocities
                    if (IAND(FLAG(i-1,j-1),C_E) /= 0) then
                        W(i-1,j-1) = W(i,j-1)+dxw/dzs*(U(i-1,j)-U(i-1,j-1))
                    end if
                    if (IAND(FLAG(i+1,j-1),C_E) /= 0) then
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
                CASE (C_NSW)
                    U(i-1,j) = U(i,j)+dx(i)/dz(j)*(W(i,j)-W(i,j-1))
                    W(i,j)   = W(i,j) 
                    W(i,j-1) = W(i,j-1) 
                    if (IAND(FLAG(i-1,j-1),C_E) /= 0) then
                        W(i-1,j-1)  = W(i,j-1)
                        U(i-1,j-1)  = U(i-1,j)
                    end if
                    if (IAND(FLAG(i-1,j+1),C_E) /= 0) then
                        W(i-1,j)   = W(i,j)
                        U(i-1,j+1) = U(i-1,j)
                    end if
                CASE (C_SWO)
                    W(i,j-1) = W(i,j)+dz(j)/dx(i)*(U(i,j)-U(i-1,j))
                    U(i,j)   = U(i,j)+delt*GX 
                    U(i-1,j) = U(i-1,j)+delt*GX
                    if (IAND(FLAG(i-1,j-1),C_E) /= 0) then
                        U(i-1,j-1) = U(i-1,j)
                        W(i-1,j-1) = W(i,j-1)
                    end if
                    if (IAND(FLAG(i+1,j-1),C_E) /= 0) then
                        U(i,j-1)    = U(i,j)
                        W(i+1,j-1)  = W(i,j-1)
                    end if
                CASE (C_NSO)
                    U(i,j) = U(i-1,j)-dx(i)/dz(j)*(W(i,j)-W(i,j-1))
                    W(i,j)   = W(i,j)+delt*GZ 
                    W(i-1,j) = W(i-1,j)+delt*GZ 
                    if (IAND(FLAG(i-1,j+1),C_E) /= 0) then
                        U(i-1,j+1) = U(i-1,j)-dzn/dxw*(W(i,j)-W(i-1,j))
                    end if
                    if (IAND(FLAG(i-1,j-1),C_E) /= 0) then
                        U(i-1,j-1) = U(i-1,j)+dzs/dxw*(W(i,j-1)-W(i-1,j-1))
                    end if
                    if (IAND(FLAG(i+1,j-1),C_E) /= 0) then
                        U(i,j-1)   = U(i,j)
                        W(i+1,j-1) = W(i,j-1)
                    end if
                    if (IAND(FLAG(i+1,j+1),C_E) /= 0) then
                        U(i,j+1)    = U(i,j)
                        W(i+1,j)    = W(i,j)
                    end if
                CASE (C_NSWO)
                    U(i,j)   = U(i,j)+delt*GX   
                    U(i-1,j) = U(i-1,j)+delt*GX 
                    W(i,j)   = W(i,j)+delt*GZ   
                    W(i,j-1) = W(i,j-1)+delt*GZ 
                    if (IAND(FLAG(i-1,j+1),C_E) /= 0) then
                        U(i-1,j+1) = U(i-1,j)
                        W(i-1,j)   = W(i,j)
                    end if
                    if (IAND(FLAG(i+1,j+1),C_E) /= 0) then
                        U(i,j+1) = U(i,j)
                        W(i+1,j) = W(i,j)
                    end if
                    if (IAND(FLAG(i-1,j-1),C_E) /= 0) then
                        U(i-1,j-1)  = U(i-1,j)
                        W(i-1,j-1)  = W(i,j-1)
                    end if
                    if (IAND(FLAG(i+1,j-1),C_E) /= 0) then
                        U(i,j-1)    = U(i,j)
                        W(i+1,j-1)  = W(i,j-1)
                    end if
                case default
                end select
            end if
        end do
    end do
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
                        !if ((iand(FLAG(i-1,j),C_F)/=C_F).AND.((IAND(FLAG(i,j+1),C_E)/=0).OR.(IAND(FLAG(i,j-1),C_E)/=0))) then
                        if ((iand(FLAG(i-1,j),C_F)/=C_F).OR.(iand(flagib(i-1,j),C_F)/=C_F)) then
                            yim=yi+dxw/dtan(rangle)
                        end if
                        !right cell is B cell
                        if ((iand(FLAG(i+1,j),C_F) /= C_F).OR.(iand(flagib(i+1,j),C_F) /= C_F)) then
                            yip=yi+dxe/dtan(rangle)
                        end if
                        !bottom cell is B cell
                        if ((iand(FLAG(i,j-1),C_F) /= C_F).OR.(iand(flagib(i,j-1),C_F) /= C_F)) then
                            !if the left cell is an empty cell
                            if ((IAND(FLAG(i-1,j),C_E)/=0).OR.(IAND(flagib(i-1,j),C_E)/=0)) then 
                                yim=yi-dxw*dtan(rangle)
                            end if
                            !if the right cell is an empty cell
                            if ((IAND(FLAG(i+1,j),C_E)/=0).OR.(IAND(flagib(i+1,j),C_E)/=0)) then 
                                yip=yi-dxe*dtan(rangle)
                            end if
                        end if
                        if ((iand(FLAG(i,j+1),C_F) /= C_F).OR.(iand(flagib(i,j+1),C_F) /= C_F)) then
                            !if the left cell is an empty cell
                            if ((IAND(FLAG(i-1,j),C_E)/=0).OR.(IAND(flagib(i-1,j),C_E)/=0)) then 
                                yim=yi-dxw*dtan(rangle)
                            end if
                            !if the right cell is an empty cell
                            if ((IAND(FLAG(i+1,j),C_E)/=0).OR.(IAND(flagib(i+1,j),C_E)/=0)) then 
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
                        if ((iand(FLAG(i-1,j),C_F)/=C_F).OR.(iand(flagib(i-1,j),C_F)/=C_F)) then
                            !if the upper cell is an empty cell
                            if ((IAND(FLAG(i,j+1),C_E)/=0).OR.(IAND(flagib(i,j+1),C_E)/=0)) then 
                                xjp=xj-dzn*dtan(rangle)
                            end if
                            !if the below cell is an empty cell
                            if ((IAND(FLAG(i,j-1),C_E)/=0).OR.(IAND(flagib(i,j-1),C_E)/=0)) then 
                                xjm=xj-dzs*dtan(rangle)
                            end if
                        end if
                        !right cell is B cell
                        if ((iand(FLAG(i+1,j),C_F) /= C_F).OR.(iand(flagib(i+1,j),C_F) /= C_F)) then
                            !if the upper cell is an empty cell
                            if ((IAND(FLAG(i,j+1),C_E)/=0).OR.(IAND(flagib(i,j+1),C_E)/=0)) then 
                                xjp=xj-dzn*dtan(rangle)
                            end if
                            !if the below cell is an empty cell
                            if ((IAND(FLAG(i,j-1),C_E)/=0).OR.(IAND(flagib(i,j-1),C_E)/=0)) then 
                                xjm=xj-dzs*dtan(rangle)
                            end if
                        end if
                        !bottom cell is B cell
                        if ((iand(FLAG(i,j-1),C_F) /= C_F).OR.(iand(flagib(i,j-1),C_F) /= C_F)) then
                            xjm=xj+dzs/dtan(rangle)
                        end if
                        if ((iand(FLAG(i,j+1),C_F) /= C_F).OR.(iand(flagib(i,j+1),C_F) /= C_F)) then
                            xjp=xj+dzn/dtan(rangle)
                        end if
                        dhn=(xjp-xj)/dzn
                        dhs=(xj-xjm)/dzs
                        curv=(dhn/dsqrt(1.0d0+dhn**2)-dhs/dsqrt(1.0d0+dhs**2))/dz(j)
                    end if
                    P(i,j)=-curv*sten
                end if
            else
                if (.NOT.((IAND(FLAG(i,j),C_E)/=0).OR.(FLAG(i,j)<C_O))) then
                    select case (iand(FLAG(i,j),C_NSWO))
                    case (C_N)
                        P(i,j) = 2.*Fr**2/Re/dz(j)*(W(i,j)-W(i,j-1))
                    case (C_S)
                        P(i,j) = 2.*Fr**2/Re/dz(j)*(W(i,j)-W(i,j-1))
                    case (C_O)
                        P(i,j) = 2.*Fr**2/Re/dx(i)*(U(i,j)-U(i-1,j))
                    case (C_W)
                        P(i,j) = 2.*Fr**2/Re/dx(i)*(U(i,j)-U(i-1,j))
                    case (C_NO)
                        P(i,j) = 1.*Fr**2/Re/2.* &
                        ( (U(i,j)+U(i-1,j)-U(i,j-1)-U(i-1,j-1))/dzs &
                        +   (W(i,j)+W(i,j-1)-W(i-1,j)-W(i-1,j-1))/dxw )
                    case (C_NW)
                        P(i,j) = -1.*Fr**2/Re/2.* &
                        ( (U(i,j)+U(i-1,j)-U(i,j-1)-U(i-1,j-1))/dzs &
                        +   (W(i+1,j)+W(i+1,j-1)-W(i,j)-W(i,j-1))/dxe )
                    case (C_SW)
                        P(i,j) = 1.*Fr**2/Re/2.* &
                        ( (U(i,j+1)+U(i-1,j+1)-U(i,j)-U(i-1,j))/dzn &
                        +   (W(i+1,j)+W(i+1,j-1)-W(i,j)-W(i,j-1))/dxe)
                    case (C_SO)
                        P(i,j) = -1.*Fr**2/Re/2.* &
                        ( (U(i,j+1)+U(i-1,j+1)-U(i,j)-U(i-1,j))/dzn &
                        +   (W(i,j)+W(i,j-1)-W(i-1,j)-W(i-1,j-1))/dxw)
                    case (C_WO)
                        P(i,j) = 0.
                    case (C_NS)
                        P(i,j) = 0.
                    case (C_NWO)
                        P(i,j) = 0.
                    case (C_NSW)
                        P(i,j) = 0.
                    case (C_SWO)
                        P(i,j) = 0.
                    case (C_NSO)
                        P(i,j) = 0.
                    case (C_NSWO)
                        P(i,j) = 0.
                    case default
                    end select
                end if
            end if
        end do
    end do

    return
    end subroutine SET_UVP_SURFACE_VOF
    
     !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine MARK_CELLS2 (C, FLAG, imax, jmax, ifull)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    ! Marj the cells of the fluid domain
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    use defs
    implicit none
    integer, intent(IN) :: imax, jmax
    integer, intent(OUT) :: ifull
    integer(I2B), dimension(0:,0:), intent(INOUT) :: FLAG
    real(DP), dimension(0:,0:), intent(INOUT) :: C
    !
    !  local variables
    !
    integer :: i, j, k
    real(RP) :: x, y

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
    do i = 0,imax+1
        do j = 0,jmax+1
            if (.NOT.(FLAG(i,j)<C_F).AND.c(i,j)>=emf) then
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

                if (IAND(FLAG(i,j-1),C_E) /= 0) then
                    FLAG(i,j) = IOR(FLAG(i,j),C_S)
                end if

                if (IAND(FLAG(i,j+1),C_E) /= 0) then
                    FLAG(i,j) = IOR(FLAG(i,j),C_N)
                end if

                if (FLAG(i,j) < C_O) then
                    ifull = ifull + 1
                end if
            end if
        end do
    end do
    

    return
    end subroutine MARK_CELLS2


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine INIT_VOLUME(problem, C, FLAG, imax,jmax, dx, dz, H)

    use nrtype
    use defs

    implicit none
    character (LEN=30), INTENT(IN) :: problem
    integer, intent(IN) :: imax, jmax
    real(RP), dimension(0:), intent(IN):: dx, dz, H
    integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
    real(DP), dimension(0:,0:), intent(INOUT) :: C

    integer i, j
    real(DP) :: zs, zn, x, z, a, b, r
    real(DP):: frac
    x=-dx(0)/2.0d0
    do i=1,imax
        x=x+(dx(i-1)+dx(i))/2.0d0
        z=-dz(0)/2.0d0
        zs = 0.0
        do j=1,jmax
            z=z+(dz(j)+dz(j+1))/2.0d0
            if (problem=="cangle") then 
                a=1.0d0
                b=0.0d0
                r=0.4d0
                if ((x-a)**2+(z-b)**2.le.r**2) then 
                    c(i,j)=1.0d0
                else
                    c(i,j)=0.0d0
                end if
            elseif (problem=="wavef") then 
                zn = zs + dz(j)
                if (immersed) then 
                    if (x.le.5.2d0) then 
                        if (H(i)>zn) then
                            c(i,j)=1.0d0
                        elseif ((zn>=H(i)) .AND. (H(i)>zs)) then
                            frac = (H(i)-zs)/(zn-zs)
                            c(i,j) = frac
                        else
                            c(i,j)=0.0d0
                        end if
                    else
                        c(i,j) = 0.0d0
                    end if
                else
                    if (H(i)>zn) then
                        c(i,j)=1.0d0
                    elseif ((zn>=H(i)) .AND. (H(i)>zs)) then
                        frac = (H(i)-zs)/(zn-zs)
                        c(i,j) = frac
                    else
                        c(i,j)=0.0d0
                    end if
                end if
                zs=zn
            else
                zn = zs + dz(j)
                if (H(i)>zn) then
                    c(i,j)=1.0d0
                elseif ((zn>=H(i)) .AND. (H(i)>zs)) then
                    frac = (H(i)-zs)/(zn-zs)
                    c(i,j) = frac
                else
                    c(i,j) = 0.0d0
                end if
                !if (immersed) then 
                !    if (IAND(flagib(i,j),C_F)==0) c(i,j)=0.0d0 
                !end if

                zs = zn
            end if

        enddo
    enddo
    
    !do i=1, imax
    !    do j=1, jmax/2
    !        C(i,j)=1.d0
    !    end do
    !end do


    do i=0, imax+1
        c(i,0) = c(i,1)
        c(i,jmax+1) = c(i,jmax)
    end do

    do j=1, jmax
        c(0,j) = c(1,j)
        c(imax+1,j) = c(imax,j)
    end do

    do i=imax/3, imax/3+4
        do j=jmax-6, jmax-2
            !C(i,j)=1.0D0
        end do
    end do
        
    return
    end subroutine
   
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRESSURE_INTERPOLATE %%%%%%%%%%%%%%%%%%%%%%%%%
!
   subroutine PRESSURE_INTERPOLATE (PS, P, c, FLAG, imax, jmax, dx, dz)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Set boundary values at free surface
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      use nrtype
	  use defs
      IMPLICIT NONE
      integer, intent(IN) :: imax, jmax
      real(RP), dimension(0:,0:), intent(IN) :: PS
      real(RP), dimension(0:,0:), intent(INOUT) :: P
      real(DP), dimension(0:,0:), intent(IN) :: c
      integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
      real(RP), dimension(0:), intent(IN):: dx, dz
!
! local variables
!
      integer :: i, j, ihv, hyd
      real(DP) :: xjp, xjm, xj, yip, yim, yi, sl1, sl2, h, d, eta, pf
      hyd=0
      do j=1, jmax
          do i=1, imax
              if (.not.((iand(FLAG(i,j),C_E)/=0).OR.(FLAG(i,j)<C_O))) then
                  yip = c(i+1,j-1)*dz(j-1)+c(i+1,j)*dz(j)+c(i+1,j+1)*dz(j+1)
                  yim = c(i-1,j-1)*dz(j-1)+c(i-1,j)*dz(j)+c(i-1,j+1)*dz(j+1)
                  xjp = c(i-1,j+1)*dx(i-1)+c(i,j+1)*dx(i)+c(i+1,j+1)*dx(i+1)
                  xjm = c(i-1,j-1)*dx(i-1)+c(i,j-1)*dx(i)+c(i+1,j-1)*dx(i+1)
                  !DY/DX
                  sl1 = 2.0*(yip-yim)/(dx(i+1)+2.0*dx(i)+dx(i-1))
                  !
                  !DX/DY
                  sl2 = 2.0*(xjp-xjm)/(dz(j+1)+2.0*dz(j)+dz(j-1))
                  if (abs(sl1)<abs(sl2)) then          !horizontal
                      if (sl2<0.d0) then               !dx/dy<0 fluid lies below the free surface
                          if ((IAND(FLAG(i,j-1),C_F)/=0).AND.(FLAG(i,j-1)<C_O)) then !is there a fluid cell below the surface?
                              h = (dz(j)+dz(j-1))/2.d0
                              d = dz(j-1)/2.d0 + c(i,j)*dz(j)
                              pf = P(i,j-1)
                          else                                                         ! else pressure is hydrostatic 
                              hyd=1
                          end if
                      else                            !dx/dy>=0 fluid lies above the free surface
                          if ((IAND(FLAG(i,j+1),C_F)/=0).AND.(FLAG(i,j+1)<C_O)) then !is there a fluid cell above the surface?
                              h = (dz(j)+dz(j+1))/2.d0
                              d = dz(j+1)/2.d0 + c(i,j)*dz(j)
                              pf = P(i,j+1)
                          else                                                         ! else pressure is hydrostatic 
                              hyd=1
                          end if
                      end if
                  else                                 !vertical
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
                      P(i,j)=max(c(i,j)*dz(j)-dz(j)/2.0,0.0) ! pressure is hydrostatic in the surface cell
                  else
                      eta = h/d
                      P(i,j) = eta*PS(i,j)+(1.d0-eta)*pf
                  end if
              end if
          end do
      end do

      return
   end subroutine PRESSURE_INTERPOLATE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	SUBROUTINE YOUNG (c0, u, v, imax, jmax, dt, dxx, dyy, wE)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  computes color field
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype
	  use defs

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: imax, jmax, wE
      REAL(RP), INTENT(IN) :: dt
      REAL(RP), DIMENSION(0:), INTENT (IN):: dxx, dyy
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: u, v
      REAL(DP), DIMENSION(0:,0:), INTENT(INOUT) :: c0
!
! local variables
!
      INTEGER :: i, j, itype
      REAL(DP), DIMENSION(0:imax+1,0:jmax+1):: c1
      REAL(RP) :: c, ft, fb, fl, fr, ut, ub, ul, ur, rnx, rny, emk2, u1, u2, u3, &
                    u4, f1, f2, f3, f4, emikro, rnx1, rny1, &
                  tanbeta, cotbeta, tanalfa, cotalfa, dx, dy
      REAL(RP), DIMENSION(0:imax+1,0:jmax+1) :: flxm, flxp, flym, flyp, div
      emikro = 1.0E-10
      emk2 = 0.0   
      dx = dxx(1)
      dy = dyy(1)
      do i=0, imax+1
          do j=0, jmax+1
              c1(i,j) = c0(i,j)
          end do
      end do
      
      do i=1, imax
          do j=1, jmax
              ft=0.0
              fb=0.0
              fl=0.0
              fr=0.0
              c=c0(i,j)
              ut=v(i,j)
              ub=v(i,j-1)
              ul=u(i-1,j)
              ur=u(i,j)
              if (c.ge.1.0-emikro) then
                  if(ut.gt.0.0) then
                      ft=ut*dt*dx*c
                  end if
                  if(ub.lt.0.0) then
                      fb=abs(ub)*dt*dx*c
                  endif
                  if(ul.lt.0.0) then
					  fl=abs(ul)*dt*dy*c
                  endif
                  if(ur.gt.0.0) then
                      fr=ur*dt*dy*c
                  endif
              else if (c.gt.0.0) then
                  rnx=(c0(i+1,j+1)+2.0*c0(i+1,j)+c0(i+1,j-1)-c0(i-1,j+1)-2.0*c0(i-1,j)-c0(i-1,j-1))/dx
                  rny=(c0(i+1,j+1)+2.0*c0(i,j+1)+c0(i-1,j+1)-c0(i+1,j-1)-2.0*c0(i,j-1)-c0(i-1,j-1))/dy
                  if(abs(rny).le.emk2) then
                      if(rnx.ge.0.0) then
                          u1=ut
                          u2=ur
                          u3=ub
                          u4=ul
                      else
                          u1=ut
                    	  u2=-ul
                       	  u3=ub
                          u4=-ur
                      end if
                      f1=0.0
                      f2=0.0
                      f3=0.0
                      f4=0.0
                      if(u1.gt.0.0) then
                          f1=abs(u1)*dt*c*dx
                      endif
                      if(u3.lt.0.0) then
                          f3=abs(u3)*dt*c*dx
                      endif
                      if(u2.gt.0.0) then
                          if(abs(u2)*dt.le.c*dx)then
                              f2=abs(u2)*dt*dy
                          else
							  f2=c*dx*dy
                          endif
                      end if
                      if(u4.lt.0.0)then
                          if(abs(u4)*dt.le.(1.0-c)*dx)then
                              f4=0.0
                          else
                              f4=(abs(u4)*dt-(1.0-c)*dx)*dy
                          endif
                      end if
                      if(rnx.ge.0.0)then
                          ft=f1
                          fb=f3
                          fr=f2
                          fl=f4
                      else
                          ft=f1
                          fb=f3
                          fr=f4
                          fl=f2
                      end if
                  else if(abs(rnx).le.emk2)then
                      if(rny.ge.0.0)then
                          u1=ut
                          u2=ur
                          u3=ub
                          u4=ul
                      else
                          u1=-ub
                          u2=ur
                          u3=-ut
                          u4=ul
                      end if
                      f1=0.0
                      f2=0.0
                      f3=0.0
                      f4=0.0
                      if(u1.gt.0.0)then
                          if(u1*dt.le.c*dy)then
                              f1=u1*dt*dx
                          else
                              f1=c*dx*dy
                          endif
                      end if
                      if(u3.lt.0.0)then
                          if(abs(u3)*dt.le.(1.0-c)*dy)then
                              f3=0.0
                          else
                              f3=(abs(u3)*dt-(1.0-c)*dy)*dx
                          end if
                      end if
                      if(u2.gt.0.0)then
                          f2=u2*dt*c*dy
                      end if
                      if(u4.lt.0.0)then
                          f4=abs(u4)*dt*c*dy
                      end if
                      if(rny.ge.0)then
                          ft=f1
                          fb=f3
                          fr=f2
                          fl=f4
                      else
                          ft=f3
                          fb=f1
                          fr=f2
                          fl=f4
                      end if
                  else
                      if(rnx.gt.0.0.and.rny.lt.0.0)then
                          u1=ut
                          u2=ur
                          u3=ub
                          u4=ul
                      else if(rnx.lt.0.0.and.rny.gt.0.0)then
                          u1=-ub
                          u2=-ul
                          u3=-ut
                          u4=-ur
                      else if(rnx.lt.0.0.and.rny.lt.0.0)then
                          u1=ut
                          u2=-ul
                          u3=ub
                          u4=-ur
                      else
                          u1=-ub
                          u2=ur
                          u3=-ut
                          u4=ul
                      endif
                      
                          rnx1=abs(rnx)
                          rny1=-abs(rny)
                          tanbeta=-rnx1/rny1
                          cotbeta=1.0/tanbeta
                          tanalfa=dx/dy*tanbeta
                          cotalfa=1.0/tanalfa
                          if(tanalfa.le.1.0)then
                              if(c.le.0.5*tanalfa)then
                                  itype=1
                              else if(c.le.1.0-0.5*tanalfa)then
                                  itype=2
                              else
                                  itype=4
                              endif
                          else
                              if(c.le.0.5*cotalfa)then
                                  itype=1
                              else if(c.le.1.0-0.5*cotalfa)then
                                  itype=3
                              else
                                  itype=4
                              endif
                          endif
						 call transport(c,itype, tanalfa, tanbeta, u1, u3, u4, u2, dx, dy, dt, f1, f3, f4, f2)
                         if(rnx.gt.0.0.and.rny.lt.0.0)then
                             ft=f1
                             fb=f3
                             fl=f4
                             fr=f2
                         else if(rnx.lt.0.0.and.rny.gt.0.0)then
                             ft=f3
							 fb=f1
							 fl=f2
							 fr=f4
                         else if(rnx.lt.0.0.and.rny.lt.0.0)then
                             ft=f1
							 fb=f3
							 fl=f2
							 fr=f4
                         else
                             ft=f3
							 fb=f1
							 fl=f4
							 fr=f2
                         endif
                  endif
              endif
              c1(i,j+1)=c1(i,j+1)+ft/dx/dy
              c1(i,j-1)=c1(i,j-1)+fb/dx/dy
              c1(i+1,j)=c1(i+1,j)+fr/dx/dy
              c1(i-1,j)=c1(i-1,j)+fl/dx/dy 
              c1(i,j)=c1(i,j)-(ft+fb+fl+fr)/dx/dy
          enddo
      enddo
   
      do i=1, imax
          do j=1, jmax
              c0(i,j)=max(min(c1(i,j),1.0d0),0.d0)
          enddo
      enddo

      do i=1,imax
         c0(i,0)      = c0(i,1)
         c0(i,jmax+1) = c0(i,jmax)
      enddo
      do j=1,jmax
         c0(0,j)=c0(1,j)
          !if (wE==3 .OR. wE==4) c0(imax,j)=c0(imax-1,j)
         c0(imax+1,j) = c0(imax,j)
      enddo

      
      
    END SUBROUTINE YOUNG
    
 SUBROUTINE transport (c, itype, tanalfa, tanbeta, u1, u3, u4, u2, dx, dy, dt, f1, f3, f4, f2)
      USE nrtype
	  use defs
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: itype
      REAL(RP), INTENT(IN) :: c, tanalfa, tanbeta, u1, u3, u4, u2, dx, dy, dt
      REAL(RP), INTENT(OUT) :: f1, f3, f4, f2
      
      REAL(RP) :: cotalfa, cotbeta, s1, s2, s3, s4, temp1
      cotalfa=1.0/tanalfa
	  cotbeta=1.0/tanbeta
	  f1=0.0
	  f3=0.0
	  f4=0.0
	  f2=0.0
	  if(itype.eq.1)then
		s1=0.0
		s3=sqrt(2.0*c*cotalfa)
		s4=0.0
		s2=sqrt(2.0*c*tanalfa)
		if(u1.gt.0.0)then
            if(u1*dt.le.(1.0-s2)*dy)then
                f1=0.0
			else
				temp1=u1*dt-(1.0-s2)*dy
				f1=0.5*temp1*temp1*cotbeta
			endif
		endif
		if(u2.gt.0.0)then
			if(u2*dt.ge.s3*dx)then
				f2=c*dx*dy
			else
				f2=0.5*u2*dt*(2.0-u2*dt/(s3*dx))*s2*dy
			endif
		endif
		if(u3.lt.0.0)then
			if((abs(u3)*dt.ge.s2*dy))then
				f3=c*dx*dy
			else
				f3=0.5*abs(u3)*dt*(2.0-abs(u3)*dt/(s2*dy))*s3*dx
			endif
		endif
		if(u4.lt.0.0)then
			if(abs(u4)*dt.le.(1.0-s3)*dx)then
				f4=0.0
			else
				temp1=abs(u4)*dt-(1.0-s3)*dx
				f4=0.5*temp1*temp1*tanbeta
			endif
		endif
	else if(itype.eq.2)then
		s1=0.0
		s3=1.0
		s4=c-0.5*tanalfa
		s2=c+0.5*tanalfa
		if(u1.gt.0.0)then
			if(u1*dt.le.(1.0-s2)*dy)then
				f1=0.0
			else if(u1*dt.le.(1.0-s4)*dy)then
				temp1=u1*dt-(1.0-s2)*dy
				f1=0.5*temp1*temp1*cotbeta
			else
				f1=u1*dt*dx-(1.0-c)*dx*dy
			endif

		endif
		if(u2.gt.0.0)then
			f2=u2*dt*(s2*dy-0.5*u2*dt*tanbeta)
		endif
		if(u3.lt.0.0)then
			if((abs(u3)*dt.le.s4*dy))then
				f3=abs(u3)*dt*dx
			else if((abs(u3)*dt.le.s2*dy))then
				temp1=abs(u3)*dt-s4*dy
				f3=abs(u3)*dt*dx-0.5*temp1*temp1*cotbeta
			else
				f3=c*dx*dy
			endif
		endif
		if(u4.lt.0.0)then
			f4=abs(u4)*dt*(s4*dy+0.5*abs(u4)*dt*tanbeta)
		endif
	else if(itype.eq.3)then
		s1=c-0.5*cotalfa
		s3=c+0.5*cotalfa
		s4=0.0
		s2=1.0
		if(u1.gt.0.0)then
			f1=u1*dt*(s1*dx+0.5*u1*dt*cotbeta)
		endif
		if(u2.gt.0.0)then
			if(u2*dt.le.s1*dx)then
				f2=u2*dt*dy
			else if(u2*dt.le.s3*dx)then
				temp1=u2*dt-s1*dx
				f2=u2*dt*dy-0.5*temp1*temp1*tanbeta
			else
				f2=c*dx*dy
			endif
		endif
		if(u3.lt.0.0)then
			f3=abs(u3)*dt*(s3*dx-0.5*abs(u3)*dt*cotbeta)
		endif
		if(u4.lt.0.0)then
			if((abs(u4)*dt.le.(1.0-s3)*dx))then
    			f4=0.0
			else if((abs(u4)*dt.le.(1.0-s1)*dx))then
				temp1=abs(u4)*dt-(1.0-s3)*dx
				f4=0.5*temp1*temp1*tanbeta
			else
				f4=abs(u4)*dt*dy-(1.0-c)*dx*dy
			endif
		endif
	else if(itype.eq.4)then
		s1=1.0-sqrt(2.0*(1.0-c)*cotalfa)
    	s3=1.0
		s4=1.0-sqrt(2.0*(1.0-c)*tanalfa)
		s2=1.0
		if(u1.gt.0.0)then
			if(u1*dt.ge.(1.0-s4)*dy)then
				f1=u1*dt*dx-(1.0-c)*dx*dy
			else
				f1=u1*dt*(s1*dx+0.5*u1*dt*cotbeta)
			endif
		endif
		if(u2.gt.0)then
			if(u2*dt.le.s1*dx)then
				f2=u2*dt*dy
			else
				temp1=u2*dt-s1*dx
				f2=u2*dt*dy-0.5*temp1*temp1*tanbeta
			endif
		endif
		if(u3.lt.0.0)then
			if(abs(u3)*dt.le.s4*dy)then
				f3=abs(u3)*dt*dx
			else
				temp1=abs(u3)*dt-s4*dy
				f3=abs(u3)*dt*dx-0.5*temp1*temp1*cotbeta
			endif
		endif
		if(u4.lt.0.0)then
			if(abs(u4)*dt.ge.(1.0-s1)*dx)then
				f4=abs(u4)*dt*dy-(1.0-c)*dx*dy
			else
				f4=abs(u4)*dt*(s4*dy+0.5*abs(u4)*dt*tanbeta)
			endif
		endif
	else
		write(*,*) 'error in subroutine transport'
		stop
	endif
    END SUBROUTINE transport
    
    subroutine VOF_COMP_DEPTH (FLAG, c, imax, jmax, dz, H)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  computes color field
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      use nrtype
	  use defs
      implicit none
      INTEGER(I2B), dimension(0:,0:), intent(IN) :: FLAG
      real(DP), dimension(0:,0:), intent(IN) :: c
      integer, intent(IN) :: imax, jmax
      real(RP), dimension(0:), intent (IN):: dz
      real(RP), dimension(0:), intent(OUT) :: H
      
      integer:: i,j
      do i=1, imax
          H(i)=0.0
          do j=1, jmax
              if ( IAND(FLAG(i,j),C_F) /= C_F) then 
                  H(i)=H(i)+dz(j)
              else
                  H(i)=H(i)+c(i,j)*dz(j)
              end if
          end do
      end do
      
      return 
      end subroutine VOF_COMP_DEPTH
             
      
    
    