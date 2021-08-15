    !
    !%%%%%%%%%%%%%%%%%%%%%%%% READ_PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine READ_PARAMETERS (FILE, problem, atype, mesh_file, meshtype, upwind, &
                                xlength, zlength, h0, TT, um, wm, imax, jmax, &
                                t_end, cfl, del_out, gamma, f_c, d_c, &
                                delx, delz, &
                                itermax, eps, omega, ro, dvis, a, GZ, &
                                wW, wE, wN, wS)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  reads parameter data and returns to main program
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    use nrtype
    implicit none
    integer, intent(OUT) ::  imax, jmax, &
    itermax, &
    wW, wE, wN, wS, f_c, d_c
    real(RP), intent(OUT) :: xlength, zlength, h0, TT, um ,wm, &
    t_end, cfl, del_out, gamma, &
    eps, omega, ro, dvis, a, GZ, delx, delz
    character (LEN=12), intent(IN) :: FILE
    character (LEN=30), intent(OUT) :: problem, atype, mesh_file, meshtype, upwind
    !
    !-----------------------------------------------------------------------
    !
    character (LEN=30) :: temp
    character (LEN=4) :: temp2
    temp2 =FILE(1:4)
    open (2, FILE = FILE, STATUS = 'OLD', FORM = 'FORMATTED')

    read (2,*) problem
    read (2,*) atype
    read (2,*) temp
    mesh_file = temp
    read (2,*) meshtype
    read (2,*) upwind
    read (2,'(1a1)')
    read (2,*) xlength
    read (2,*) zlength
    read (2,*) h0
    read (2,*) TT
    read (2,*) um
    read (2,*) wm
    read (2,*) imax
    read (2,*) jmax
    delx = xlength/imax
    delz = zlength/jmax
    read (2,'(1a1)')
    read (2,*) t_end
    read (2,*) cfl
    read (2,*) del_out
    read (2,*) gamma
    read (2,*) f_c
    read (2,*) d_c
    read (2,'(1a1)')
    read (2,*) itermax
    read (2,*) eps
    read (2,*) omega
    read (2,'(1a1)')
    read (2,*) ro
    read (2,*) dvis
    read (2,*) a
    read (2,*) GZ
    read (2,'(1a1)')
    read (2,*) wW
    read (2,*) wE
    read (2,*) wN
    read (2,*) wS

    close(2)
    return

    end subroutine READ_PARAMETERS


    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%% INIT_UWP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    SUBROUTINE INIT_UWP(problem, U, W, P, x, z, dx, dz, imax, jmax, ro, GZ, h0, &
                        H, Ls, nm, xm, ym)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  initializes U, W, P and H
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    USE nrtype
    use defs
    use interfaces2
    IMPLICIT NONE
    CHARACTER (LEN=30), INTENT(IN) :: problem
    INTEGER, INTENT(IN) :: imax, jmax, nm
    REAL(RP), INTENT(IN) :: ro, GZ, h0
    REAL(RP), INTENT(OUT) :: Ls
    REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, W, P
    REAL(RP), DIMENSION(0:), INTENT(IN) :: x, z, dx, dz
    REAL(RP), DIMENSION(0:), INTENT(INOUT) :: H
    real(DP), dimension(:), allocatable, intent(IN) :: xm, ym
    !
    !-----------------------------------------------------------------------
    !
    INTEGER :: i, j, ic, ii
    REAL(RP):: xx, zz, UI, xcrest, js, hs, d, c, k, L, t, dt, KAI, EPI, hostep, dist, distmin
    UI=1.0
    KAI=0.003
    EPI=0.0005  !
    if (problem=="circle") then 
        UI=1.0
        if (VOLUM_TRACKING) then 
            do i=1, imax
                H(i)=h0
            end do
        end if
    end if
    
    
    do i = 0,imax+1
        do j = 0,jmax+1
            U(i,j) = 0.0
            W(i,j) = 0.0
            P(i,j) = 0.0
        end do
    end do
            
    if (problem == "backstep") then

        do i=0, imax+1
            zz = -dz(0)/2.0
            do j = 0,jmax
                if (zz<0.5) U(i,j) = 0.0
                zz=zz+(dz(j)+dz(j+1))/2.0
            end do
        end do

    else if (problem=="damb") then
        do i=0, imax
            if (x(i)-dx(i)/2.0<=h0) then   ! The column width is 1.0
                !if (x(i)-dx(i)/2.0>=1.5.AND.x(i)-dx(i)/2.0<=3.5) then   ! The column width is 1.0
                !if (x(i)-dx(i)/2.0>=3.0) then   ! The column width is 1.0
                H(i) = h0
            else
                H(i) = 0.0
            endif
        end do
    else if(problem=="tank".or.problem=="reservoir".or.problem=="wavef") then
        do i=0, imax+1
            H(i)=h0
        end do
    elseif (problem=="solitary") then 
        d = h0			!initial water depth
        !hs = 0.3        !wave height A 
        hs = 0.12*d        !wave height A 
        k = sqrt(3.0/4.0*hs/d**3)
        L = 2.0/k*acosh(sqrt(1.0/0.05))
        Ls=L
        !xcrest = 2.0*L-L/2.0	!distance from wave crest to the left solid boundary
        !xcrest = 15	!distance from wave crest to the left solid boundary
        xcrest = 1.5936	!distance from wave crest to the left solid boundary
        c = sqrt(-GZ*d)*(1.0+0.5*hs/d-3.0/20.0*(hs/d)**2)
        do i=1, imax
            H(i) = d+hs/(cosh(k*(x(i)-dx(i)/2.0-xcrest)))**2
            js = J_OF_Z(jmax, H(i), dz)
            zz = -dz(0)/2.0
            do j=0, jmax
                U(i,j) = sqrt(-GZ*d)*hs/d/(cosh(k*(x(i)-xcrest)))**2
                W(i,j) = sqrt(-GZ*3.0/d)*sqrt((hs/d)**3)*z(j)/ &
                (cosh(k*(x(i)-dx(i)/2.0-xcrest)))**2* &
                tanh(k*(x(i)-dx(i)/2.0-xcrest))
                if (zz<=H(i)) then
                    P(i,j) = (H(i)-zz)*ro*(-GZ)
                end if
                zz = zz + (dz(j)+dz(j+1))/2.0
            end do
        end do
        H(0) = H(1)
        H(imax+1) = H(imax)

        open(21,FILE='fsurface2.plt',STATUS='UNKNOWN')
        t=0.0
        dt=1.0
        do while (t<20.0)
            write (21,'(4hZONE 1x2hT= 1x3h"T= f8.3,1x1h" )') t
            do i=1, imax
                zz = d+hs/(cosh(k*(x(i)-dx(i)/2.0-xcrest-c*t)))**2
                write (21,50) x(i)-dx(i)/2.0, zz
50              format (f10.5,5x,f13.8)
            end do
            t = t + dt
        end do

    end if
    if (problem /= "dcavity" .AND. problem /= "backstep" .AND. problem/="circle") then
        xx=-dx(0)/2.0
        do i=1, imax
            xx=xx+(dx(i)+dx(i-1))/2.0
            zz = -dz(0)/2.0
            do j=0, jmax
                !if (immersed) then 
                !    distmin=1.0e8
                !    do ii=1, nm
                !        dist=dsqrt((xx-xm(ii))**2+(zz-ym(ii))**2)
                !        if (dist<distmin) then 
                !            distmin=dist
                !            ic=ii    !arclength index of the nearest marker point
                !        end if
                !    end do
                !    if (zz.le.H(i).AND.zz.ge.ym(ic).AND.flagib(i,j)/=C_B) then
                !        P(i,j) = (H(i)-zz)*ro*(-GZ)
                !    end if
                !else
                    if (zz<=H(i)) then
                        P(i,j) = (H(i)-zz)*ro*(-GZ)
                    endif
                !end if
                zz = zz + (dz(j)+dz(j+1))/2.0
            end do
        end do
    end if
    !
    !   left and right boundaries
    !
    do j=0, jmax+1
        P(0,j) = P(1,j)
        P(imax+1,j) = P(imax,j)
    end do
    H(0)=H(1)
    H(imax+1)=H(imax)
    
    RETURN
    END SUBROUTINE INIT_UWP

    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INIT_FLAG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    subroutine INITFLAG (uflag, vflag, pflag, imax, jmax, ibound, vos)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !   Defines the computational geometry
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    USE nrtype
    use defs
    IMPLICIT NONE
    integer, intent(IN) :: imax, jmax
    integer(I2B), dimension(0:,0:), intent(INOUT) :: uflag, vflag, pflag
    integer, intent(OUT) :: ibound
    real(DP), dimension(0:,0:), intent(INOUT) :: vos
    !
    !-----------------------------------------------------------------------
    !
    INTEGER :: i, j
    REAL(RP):: x, z, zz, hostep, xb, zb, xt, zt, dist, m, alpha, d, hs, k, L, Lb, y
    logical inclined
    integer(I2B), dimension(0:imax+1,0:jmax+1) :: FLAG
    !
    !  initialize boundary cells
    !
    do i = 0,imax+1
        FLAG(i,0) = C_B
        FLAG(i,jmax+1) = C_B
    end do
    do j = 1,jmax
        FLAG(0,j) = C_B
        FLAG(imax+1,j) = C_B
    end do

    do i = 1,imax
        do j = 1,jmax
            FLAG(i,j) = C_F
        end do
    end do
    !do i=1, imax
    !    do j=1, jmax
    !        if (i>139.AND.i<151) then 
    !            if (j<11) then 
    !                FLAG(i,j) = C_B
    !            end if
    !        end if
    !    end do
    !end do
    
    
    ibound=0
    do i = 1, imax
        do j = 1, jmax
            if ( IAND(FLAG(i,j),C_F) /= C_F) then
                ibound = ibound + 1
            end if
            pflag(i,j) = FLAG(i,j) + ( IAND(FLAG(i-1,j),C_F) * B_W  &
            +   IAND(FLAG(i+1,j),C_F) * B_E  &
            +   IAND(FLAG(i,j-1),C_F) * B_S  &
            +   IAND(FLAG(i,j+1),C_F) * B_N ) / C_F
        end do
    end do   
    
    if (immersed) then 
        do i=0, imax+1
            do j=0, jmax+1
                uflag(i,j)=FLAG(i,j)
                vflag(i,j)=FLAG(i,j)
                pflag(i,j)=FLAG(i,j)
                vos(i,j)=0.0d0
            end do
        end do
    end if
               
    return
    end subroutine INITFLAG

    !
    !%%%%%%%%%%%%%%%%%%%%%%%% DIMENSIONLESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    SUBROUTINE DIMENSIONLESS (problem, atype, imax, jmax, xlen, zlen, delx, delz, um, wm, vm, dvis, &
                                a, ro, GZ, Re, Fr, St, Fa, t_end, TT, h0, cfreq)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  Compute dimensionless numbers
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    USE nrtype
    use defs
    IMPLICIT NONE
    CHARACTER (LEN=30), INTENT(IN) :: problem, atype
    INTEGER, INTENT(IN) :: imax, jmax
    REAL(RP), INTENT(IN) :: dvis, a
    REAL(RP), INTENT(INOUT) :: GZ, xlen, zlen, delx, delz, um, wm, vm, ro, h0
    REAL(RP), INTENT(OUT) :: Re, Fr, St, Fa, t_end, TT, cfreq
    !
    !-----------------------------------------------------------------------
    !
    REAL(RP):: pi, d, hs, c, AA, KC
    pi = 4.0*atan(1.0)

    if (atype == "dimensional") then
        St = 1.0
        Fr = 1.0
        Fa = a
        Re = 1.0/(dvis/ro)
    else
        if (problem == "damb") then
            vm = sqrt(-GZ*h0)
        else
            vm = max(um,wm)
        end if
        if (problem == "damb" .OR. problem == "backstep" .OR. problem == &
            "dcavity" .OR. problem=="hjump") then
            Fr = 1.0
            St = 1.0
            Fa = a/sqrt((-GZ)*h0)
            Re = ro*vm*h0/dvis
        elseif (problem == "solitary") then
            d = h0
            hs = 0.2*d
            c = sqrt(-GZ*d)*(1.0+0.5*hs/d-3.0/20.0*(hs/d)**2)
            St = 1.0
            Fr = 1.0
            Fa = a/sqrt((-GZ)*h0)
            Re = ro*c*h0/dvis
        elseif (problem=="circle") then
            if (movbound) then 
                KC=5.0      !Keulegan-Carpenter Number
                AA=KC/2.0/pi
                vm=2.0*pi*AA/TT
            else
                vm = max(um,wm)
            end if
            
            Fr = 1.0
            St = Vm*TT/h0
            Fa = a/vm
            Re = ro*vm*h0/dvis
        else
            Fa = a/sqrt((-GZ)*h0)
            St = vm*TT/h0
            Fr = vm/sqrt((-GZ)*h0)
            Re = ro*vm*h0/dvis
        endif

        xlen = xlen/h0
        zlen = zlen/h0
        t_end = t_end/TT
        delx = delx/h0
        delz = delz/h0
        h0=h0/h0
        um = um/vm
        wm = wm/vm
        vm	= vm/vm
        TT = TT/TT
        if (GZ/=0.0) then
            GZ = GZ/(-GZ)
        else
            GZ=0.0
        endif
        ro = ro/ro
        if (problem == "backstep") then
		    Fa = 1.0
		    GZ = 0.0
		    Fr = 1.0
		end if
    endif
    cfreq = 2.0*pi/TT
    RETURN
    END SUBROUTINE DIMENSIONLESS
