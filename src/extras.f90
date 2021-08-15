    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WAVE_PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine WAVE_PARAMETERS (problem, Fr, St, GZ, TT, h0, c, lambda)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    ! Compute wave speed and wave length
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    implicit none
    character (LEN=30) :: problem
    real(RP), intent(IN) :: Fr, St, GZ, TT, h0
    real(RP), intent(OUT) :: c, lambda
    !
    !-----------------------------------------------------------------------
    !
    real(RP) :: pi
    pi = 4.0*atan(1.0)
    lambda = (St/Fr)**2*(-GZ)*TT**2/2.0/pi
    ! wave speed computation
    !
    if (problem=="damb" .OR. problem == "tank" .OR. problem=="hjump" .OR. problem=="solitary") then
        c = sqrt(-GZ*h0)/Fr
    else if (problem == "wavef" .OR. problem == "reservoir") then
        if (h0/lambda<0.05) then
            !
            !Shallow water
            !
            c = sqrt(-GZ*h0)/Fr
        else
            !
            !deep water
            !
            c = St/Fr**2*(-GZ)*TT/2.0/pi
        end if
    end if

    return 
    end subroutine WAVE_PARAMETERS


!
!%%%%%%%%%%%%%%%%%%%%%%%% GRID_GEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
	subroutine GRID_GEN( problem, mesh_file, imax, jmax, xlen, zlen, dx, dz, h0, lambda, xtoe)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  generates non-uniform rectangular grid
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
	use nrtype
    use defs
	implicit none
	character (LEN=30), intent(IN) :: problem, mesh_file
	integer, intent (OUT) :: imax, jmax
	real(RP), intent(IN):: h0, lambda
    real(RP), intent(INOUT):: xlen, zlen
	real(RP), dimension(0:), intent(OUT):: dx, dz
    real(RP), intent(OUT):: xtoe
!
!-----------------------------------------------------------------------
!
	integer:: i, j, k, Nx, Nz, np, temp, nxx, nzz, ii
	real(RP):: betax, betaz, sum, sumx, sumz, dxmin, dymin, sume, Ls, dxp, dzp
    real(RP), dimension(0:1000) ::dxx, dzz
    real(RP):: k1, dzmin, dxmax, dzmax, alpha, betas, betab, wave_height, dzw, dxl
    real(RP):: sumz_surf, dzs, cut, xx, zz, zu, left, top, bottom, center
    
    open(11,FILE='grid.plt',STATUS='UNKNOWN')
    betax=0.90
    betaz=0.90
    Nx=25   !number of grid points in the clustered region 
    Nz=30   !number of grid points in the clustered region 
    dxmin=0.0239
    dymin=0.0239
    !dxmin=0.02
    !dymin=0.02
    center=2.0
    left=8.0
    top=13.0
	if (problem == "circle") then
        if (movbound) then 
            left=xlen/2.0-center
            top=zlen/2.0-center
        end if
        nxx=int(center/dxmin)
        sumx=0.0
        do i=1, nxx
            dxx(i)=center/nxx
            sumx=sumx+dxx(i)
        end do
        sum=0
        do i=0, Nx-2
            sum=sum+betax**i
        end do
        dxx(nxx+Nx-1)=left/sum
        temp=nxx
        sumx=sumx+dxx(nxx+Nx-1)
        nxx=nxx+1
        do i=temp+Nx-2, nxx, -1
            dxx(i)=dxx(i+1)*betax
            sumx=sumx+dxx(i)
            nxx=nxx+1
        end do
        sumx=0.0
        do i=1, nxx
            dx(i)=dxx(nxx+1-i)
            sumx=sumx+dx(i)
        end do
        dx(0)=dx(1)
        temp=0
        do i=nxx+1, 2*nxx-1
            temp=temp+1
            dx(i)=dx(nxx-temp)
            sumx=sumx+dx(i)
        end do
        do while(sumx<=40.0) 
            dx(i)=dx(i-1)/betax
            sumx=sumx+dx(i)
            i=i+1
        end do
        imax=i-2
        nzz=int(center/dymin)
        sumz=0.0
        do i=1, nzz
            dzz(i)=center/nzz
            sumz=sumz+dzz(i)
        end do   
        sum=0
        do i=0, Nz-2
            sum=sum+betaz**i
        end do
        dzz(nzz+Nz-1)=top/sum
        temp=nzz
        sumz=sumz+dzz(nzz+Nz-1)
        nzz=nzz+1
        do i=temp+Nz-2, nzz, -1
            dzz(i)=dzz(i+1)*betaz
            sumz=sumz+dzz(i)
            nzz=nzz+1
        end do
        sumz=0.0
        temp=0
        do i=1, nzz
            dz(i)=dzz(nzz-temp)
            sumz=sumz+dz(i)
            temp=temp+1
        end do
        dz(0)=dz(1)
        temp=0
        sume=0.0
        do j=nzz+1, 2*nzz
            dz(j)=dz(nzz-temp)
            sumz=sumz+dz(j)
            temp=temp+1
            sume=sume+dz(j)
        end do
        jmax=j-1
        dz(jmax+1)=dz(jmax)
	!if (problem == "circle") then
 !       center=2.0
 !       left=8.0
 !       top=13.0
 !       if (movbound) then 
 !           center=center
 !           left=xlen/2.0-center
 !           top=zlen/2.0-center
 !       end if
 !       betax=0.85
 !       betaz=0.85
 !       dxmin=0.024
 !       dymin=0.024
 !       nxx=int(center/dxmin)
 !       sumx=0.0
 !       do i=1, nxx
 !           dxx(i)=center/nxx
 !           sumx=sumx+dxx(i)
 !       end do
 !       do while(sumx<left+center) 
 !           dxx(i)=dxx(i-1)/betax
 !           sumx=sumx+dxx(i)
 !           i=i+1
 !       end do
 !       Nx=i-1
 !       if (sumx.gt.left+center) then 
 !           sumx=sumx-dxx(Nx)
 !           dxx(Nx)=left+center-sumx
 !       end if
 !       sumx=0.0
 !       do i=1, Nx
 !           dx(i)=dxx(Nx-i+1)
 !           sumx=sumx+dx(i)
 !       end do
 !       do while(sumx.lt.left+2.0*center) 
 !           dx(i)=dx(i-1)
 !           sumx=sumx+dx(i)
 !           i=i+1
 !       end do
 !       do while(sumx.lt.xlen) 
 !           dx(i)=dx(i-1)/betax
 !           sumx=sumx+dx(i)
 !           i=i+1
 !       end do
 !       imax=i-1
 !       dx(0)=dx(1)
 !       dx(imax+1)=dx(imax)
 !       nzz=int(center/dymin)
 !       sumz=0.0
 !       do i=1, nzz
 !           dzz(i)=center/nzz
 !           sumz=sumz+dzz(i)
 !       end do 
 !       do while(sumz.le.top+center) 
 !           dzz(i)=dzz(i-1)/betaz
 !           sumz=sumz+dzz(i)
 !           i=i+1
 !       end do
 !       NZ=i-1
 !       if (sumz.gt.top+center) then 
 !           sumz=sumz-dzz(Nz)
 !           dzz(Nz)=top+center-sumz
 !       end if
 !       sumz=0.0
 !       do i=1, Nz
 !           dz(i)=dzz(Nz-i+1)
 !           sumz=sumz+dz(i)
 !       end do
 !       temp=0
 !       do i=Nz+1, 2*Nz
 !           dz(i)=dz(Nz-temp)
 !           sumz=sumz+dz(i)
 !           temp=temp+1
 !       end do
 !       jmax=i-1
 !       dz(jmax+1)=dz(jmax)
    elseif (problem=="tank") then 
        !grid in x direction
        Nx=40
        betax=0.95
        sumx=0.0
        do i=0, Nx-1
            sumx=sumx+betax**(-i)
        end do
        dx(1)=xlen/2.0/sumx
        dx(0)=dx(1)
        sumx=dx(1)
        do i=2, Nx
            dx(i)=dx(i-1)/betax
            sumx=sumx+dx(i)
        end do
        temp=0
        do i=Nx+1, 2*Nx
            temp=temp+1
            dx(i)=dx(Nx+1-temp)
            sumx=sumx+dx(i)
        end do
        dx(2*Nx+1)=dx(2*Nx)
        imax=2*Nx
        !grid in z direction
        Nz=40
        betaz=0.98
        sumz=0.0
        do i=0, Nz-1
            sumz=sumz+betaz**(-i)
        end do
        dz(1)=zlen/4.0/sumz
        dz(0)=dz(1)
        sumz=dz(1)
        do i=2, Nz
            dz(i)=dz(i-1)/betaz
            sumz=sumz+dz(i)
        end do
        temp=0
        do i=Nz+1, 2*Nz
            temp=temp+1
            dz(i)=dz(Nz+1-temp)
            sumz=sumz+dz(i)
        end do
        temp=0
        do i=2*Nz+1, 4*Nz
            temp=temp+1
            dz(i)=dz(2*Nz+1-temp)
            sumz=sumz+dz(i)
        end do
        
        dz(4*Nz+1)=dz(4*Nz)    
        jmax=4*Nz 
    else if (problem=="reservoir") then 
        open (3, FILE = mesh_file, STATUS = 'OLD', FORM = 'FORMATTED')
        read (3,*) k1
        read (3,*) dxmin
        read (3,*) dxmax
        read (3,*) dzmin
        read (3,*) dzmax
        read (3,*) alpha
        read (3,*) betas
        read (3,*) betab
        close(3)
        if (immersed) then
            zu=h0-0.1*h0
            nz=int(zu/dzmin)
            dzp=zu/nz      !uniform in z direction
            j=0
            zz=0.0
            do while (zz.le.zu) 
                j=j+1
                dz(j)=dzp
                zz=zz+dzp
            end do
            do while(zz.le.h0)
                j=j+1
                dz(j)=dz(j-1)*betas
                zz=zz+dz(j)
            end do
            do while(zz.le.zlen) 
                j=j+1
                dz(j)=dz(j-1)/betas
                zz=zz+dz(j)
            end do
            jmax=j
            dz(0)=dz(1)
            dz(jmax+1)=dz(j)
            Ls=19.0*h0     !Location of the dam-toe
            xx=dxmin
            dxp=dxmin
            do while(xx.le.Ls)
                dxp=dxp/alpha
                xx=xx+dxp
            end do
            xx=dxp
            i=1
            dx(i)=dxp
            do while (dxp.ge.dxmin)
                i=i+1
                dxp=dxp*alpha
                dx(i)=dxp
                xx=xx+dxp
            end do
            xtoe=xx+dx(i)/3.0/h0   ! shift the dam on the grid to overcome the possible ambigiuties on the grid
            !do while (i.le.ii+10)
            !    i=i+1
            !    dx(i)=dx(i-1)/alpha
            !    xx=xx+dx(i)
            !end do
            do while (xx.lt.xlen)
                i=i+1
                dx(i)=dx(i-1)
                xx=xx+dx(i)
            end do
            xlen=xx
            imax=i 
            dx(0)=dx(1)
            dx(imax+1)=dx(imax)
        else
            nx = 10
            nz = 40
            !wave_height = wh(casen)
            !wave_height = 25.0
            wave_height = 0.01
            dzw = wave_height/nz
            !
            !	uniform region in 2lambda
            !
            dxl = lambda/nx
            !dxx(0) = dxl
            !dxx(1) = dxl
            dxx(0) = dxmin
            dxx(1) = dxmin
            sumx = dxl
            i = 1
            loopx1: do
                if (sumx>=2.0*lambda) exit loopx1
                i = i + 1
                !dxx(i) = dxl
                dxx(i) = dxmin
                sumx = sumx + dxx(i)
            end do loopx1
            !
            !	cluster along x_axis
            !
            loopx: do
                if (sumx>=xlen) then
                    exit loopx
                else
                    i = i + 1
                end if
                if (dxx(i-1)<dxmax) then
                    dxx(i) = dxx(i-1)/alpha
                ELSE
                    dxx(i) = dxmax
                end if
                sumx = sumx + dxx(i)
            end do loopx
            imax = i
            dxx(imax+1) = dxx(imax)
            !
            !   cluster along z_axis
            !
            sumz_surf = wave_height
            dzs = dzmax
            k = 0
            loop_surf: do
                if (dzs<=dzw) then
                    exit loop_surf
                else
                    dzs = dzs*betas
                    k = k + 1
                    sumz_surf = sumz_surf + dzs
                end if
            end do loop_surf
            dzz(0) = dzmin
            dzz(1) = dzmin
            sumz = dzz(1)
            j = 1
            cut = 1.0
            loopz : do
                if (dzz(j)>=dzmax .OR. sumz>h0-sumz_surf) then
                    if (sumz>h0-sumz_surf) then
                        cut = 1.1 ! erken çýkarsa önlem al
                    end if
                    exit loopz
                else
                    j = j + 1
                    dzz(j) = dzz(j-1)/ betab
                    sumz =sumz + dzz(j)
                end if
            end do loopz

            loopz2 : do
                if (sumz>h0-sumz_surf .OR. cut /= 1.0) then
                    exit loopz2
                else
                    j = j + 1
                    dzz(j) = dzz(j-1)
                    sumz =sumz + dzz(j)
                end if
            end do loopz2
            loopz3 : do
                if (sumz>h0-wave_height) then
                    exit loopz3
                else
                    j = j + 1
                    dzz(j) = dzz(j-1)*betas
                    sumz =sumz + dzz(j)
                end if
            end do loopz3
            loopz4 : do
                if (sumz>h0+wave_height) then
                    exit loopz4
                else
                    j = j + 1
                    dzz(j) = dzw
                    sumz =sumz + dzz(j)
                end if
            end do loopz4
            jmax = j
            dzz(jmax+1) = dzz(jmax)
            do i=0, imax+1
                dx(i)=dxx(i)
            end do
            do j=0, jmax+1
                dz(j)=dzz(j)
            end do
        end if
        
        write (11,*) 'VARIABLES= ','"X",','"Z"'
        write (11,*) 'ZONE i=',imax,',j=',jmax
        zz=0.0
        do j=1, jmax
            xx=0.0
            do i=1, imax
                write(11,*) xx, zz
                xx=xx+dx(i)
            end do
            zz = zz+dz(j)
        end do
        	
    end if

	return
    end subroutine GRID_GEN
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHOPRA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine CHOPRA (jmax, h0, a, t, T0, alpha, ro, GZ, z)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  Compute hydrodynamic pressure according to Chopra approximation
    !
    !  M	: Number of nodes on the dam interface
    !  h0	: Water depth at the reservoir
    !  a	: Sonic velocity
    !  t	: Time level to compute the hydrodynamic pressure
    !  TT	: Period of the harmonic motion
    !  alpha: Horizontal ground acceleration, in units of G
    !  ro	: Density of the fluid
    !  GZ	: Gravitational acceleration
    !  PC	: Computed hydrodynamic pressures on the dam face
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    implicit none
    integer, intent(IN) :: jmax
    real(RP), intent(IN) :: h0, a, t, T0, alpha, ro, GZ
    real(RP), dimension(0:), INTENT(IN) :: z
    !
    !-----------------------------------------------------------------------
    !
    integer :: j, n, n1, tstep
    real(RP):: dz, sum, lambda, pi, freq, part1, part2, tt, delt, t_end, p, pc, ps, zz
    pi = 4.0*atan(1.0)
    freq = 2.0*pi/T0
    tstep = 500
    t_end = 20.0*T0
    delt = t_end/tstep
    ps = 0.4*ro*(-GZ)*h0
    OPEN(12,FILE='dam-face.plt',STATUS='UNKNOWN')
    OPEN(10,FILE='dam-base.plt',STATUS='UNKNOWN')
    write (12,*)	'ZONE T= ','CHOPRA'
    write (10,*)	'ZONE T= ','CHOPRA'
    do j=1, jmax
        part1 = 0.0
        part2 = 0.0
        if (z(j)>h0) exit
        n1 = 1
        findn1:	DO
            lambda = (2*n1-1)*pi/2.0/h0
            if (lambda**2>freq**2/a**2) then
                exit findn1
            else
                n1 = n1 + 1
            end if
        end do findn1
        
        do n=1, n1-1
            lambda = (2*n-1)*pi/2.0/h0
            part1 = part1 + (-1)**(n-1)/(2*n-1)/sqrt(freq**2/a**2-lambda**2)*cos(lambda*z(j))
        end do
        do n=n1, 1000
            lambda = (2*n-1)*pi/2.0/h0
            part2 = part2 + (-1.0)**(n-1)/(2*n-1)/sqrt(lambda**2-freq**2/a**2)*cos(lambda*z(j))
        end do
        pc = 4.0*alpha*ro*(-GZ)/pi*(part1*sin(freq*t)+part2*cos(freq*t))/ps
        write (12,*) pc, z(j)/h0
    end do
    part1 = 0.0
    part2 = 0.0
    tt = 0.0
    zz=0.0
    do while (tt <= t_end)
        part1=0.0
        part2=0.0
        n1 = 1
        findn12:	do
            lambda = (2*n1-1)*pi/2.0/h0
            if (lambda**2>freq**2/a**2) then
                exit findn12
            else
                n1 = n1 + 1
            end if
        end do findn12
        do n=1, n1-1
            lambda = (2*n-1)*pi/2.0/h0
            part1 = part1 + (-1)**(n-1)/(2*n-1)/sqrt(freq**2/a**2-lambda**2)
        end do
        do n=n1, 1000
            lambda = (2*n-1)*pi/2.0/h0
            part2 = part2 + (-1.0)**(n-1)/(2*n-1)/sqrt(lambda**2-freq**2/a**2)
        end do
        pc = 4.0*alpha*ro*(-GZ)/pi*(part1*sin(freq*tt)+part2*cos(freq*tt))/ps
        write (10,*) tt, pc
        tt = tt + delt
    end do

    return
    end subroutine CHOPRA
    
!
!%%%%%%%%%%%%%%%%%%%%%%%% EL_CENTRO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
	SUBROUTINE READ_EARTHQUAKE_ACCEL(EC, delt)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Read El Centro earthquake acceleration record and store to EA array
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
	  USE nrtype
	  IMPLICIT NONE
	  REAL(RP), DIMENSION(0:), INTENT(OUT):: EC
	  REAL(RP), INTENT(OUT):: delt
	  INTEGER :: i, j, k
	  REAL(RP) :: t, a
      OPEN (18, FILE = 'DZC180.AT2', STATUS = 'OLD', FORM = 'FORMATTED')
      !OPEN (18, FILE = 'elcentro_EW.dat', STATUS = 'OLD', FORM = 'FORMATTED')

	  delt = 0.005
      !READ(18,*)
	  !READ(18,*)
	  !READ(18,*)
	  !READ(18,*) 
      !READ(18,*) 
	  100 FORMAT(5E15.6)
	  EC(0) = 0.0
	  READ(18,*, END=20) (EC(i), i=0,20000)
      !READ(18,100, END=20) (EC(i), i=0,20000)
	  20 CONTINUE
	  OPEN(25,FILE='accel.plt',STATUS='UNKNOWN')
	  t = 0.0
   50 format (f8.5,1x,f10.8)

	  DO j=0, i-1
		WRITE(25,*) t, ',',EC(j)
		t = t + delt
	  END DO
	  CLOSE(18)
	  RETURN
      END SUBROUTINE READ_EARTHQUAKE_ACCEL
      
!
!%%%%%%%%%%%%%%%%%%%%%%%% SEISMIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
	REAL FUNCTION SEISMIC(t, dt, EC)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Computes the acceleration at time t using an interpolation technique
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
	  USE nrtype
	  IMPLICIT NONE
	  REAL(RP), INTENT(IN):: t, dt
	  REAL(RP), DIMENSION(0:), INTENT(IN):: EC
	  INTEGER :: i, ii
	  REAL :: tt
	  tt = 0.02
	  i = t/dt
	  ii = tt/dt
	  SEISMIC = ((t-i*dt)*EC(i+1)+((i+1)*dt-t)*EC(i))/dt
	  RETURN
	END FUNCTION SEISMIC