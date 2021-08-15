    !  NASSLARDBM.f90 
    !
    !  FUNCTIONS:
    !  NASSLARDBM - Entry point of console application.
    !

    !****************************************************************************
    !
    !  PROGRAM: NASSLARDBM
    !
    !  PURPOSE:  Entry point for the console application.
    !
    !****************************************************************************

    program NASSLARDBM
    

    use nrtype
    use interfaces
    use defs
    implicit none
    integer :: imax, jmax, f_c, d_c, itermax, wW, wE, wNN, wS, i, j, ibound, ifull, iter, nm, nbx, nbz, n, &
                step, nmax, ppu, ppv, ppp, ivp, jvp, k, zn
    integer(I2B), dimension(:,:), allocatable :: uflag, vflag, pflag
    real(RP), dimension(:,:), allocatable :: UNP1, WNP1, UN, WN, UNM1, WNM1, F, G, RHS, PS, dpp, P, ph
    REAL(DP), dimension(:,:), allocatable :: C
    real(RP), dimension(:), allocatable   :: dx, dz, dxx, dzz, x, z, H, EC
    
    real(RP) :: xlen, zlen, h0, TT, um, wm, vm, t_end, cfl, del_out, gamma, delx, delz, eps, omega, ro, dvis, a, GZ, &
                t, delt, Re, St, Fr, Fa, res, cfreq, cc, lambda, ts, axn, axnm1, hhi, hhip1, Ls, uu, vv, tstr,  xtoe, &
                dam_height, alpha, sw, delt_rec, pi, q
    real(DP), dimension(:), allocatable   :: xm, ym, sc, nxm, nym, ux, uy, vx, vy, xc, yc, au, bu, cu, av, bv, cv, ap, bp, cp
    real(DP), dimension(:,:,:), allocatable :: invagu, invagv, invagp, invaguin, invagvin
    real(DP), dimension(:,:), allocatable :: ghu, ghv, ghp, ghuin, ghvin, vos
    integer, dimension(:,:,:), allocatable :: gFIu, gFIv, gFIp, gFIuin, gFIvin
    integer, dimension(:,:), allocatable :: indu, indv, indp, ua, va, induin, indvin
    integer, dimension(:), allocatable :: IMu, IMv, IMp
    real(RP) :: uxb, uyb
    real(DP) :: arclen
    
    character (LEN=30) :: problem, FNAME0, mesh_file, upwind, atype, mesh_type
    logical found

    OPEN(35,FILE='uv.plt',STATUS='UNKNOWN')
    OPEN(36,FILE='pc.plt',STATUS='UNKNOWN')
    allocate(dxx(0:10000), dzz(0:10000))
    VOLUM_TRACKING =.false.         ! choose vof method for the free-surface displacement
    immersed=.false.                ! choose immersed-boundary method for the treatment of irregular solid boundaries
    curvature=.false.               ! calculate the total curvature of free-surface for the pressure at the free-surface
    movbound=.false.                ! consider moving solid boundary
    absorption=.TRUE.              !absorption effects on the reservoir bottom
    hf=.false.                      ! use height function method to track the free-surface 
    emf = 1.0d-6
    emf1 = 1.d0-emf
    !FNAME0= "DOBS1INP.DAT"
    !FNAME0= "DCAV1INP.DAT"
    FNAME0= "RESERINP.DAT"
    !FNAME0= "DAMB1INP.DAT"
    !FNAME0= "SOLITINP.DAT"
    call READ_PARAMETERS (FNAME0, problem, atype, mesh_file, mesh_type, upwind, &
                            xlen, zlen, h0, TT, um, wm, imax, jmax, &
                            t_end, cfl, del_out, gamma, f_c, d_c, &
                            delx, delz, &
                            itermax, eps, omega, ro, dvis, a, GZ, &
                            wW, wE, wNN, wS)
    call DIMENSIONLESS (problem, atype, imax, jmax, xlen, zlen, delx, delz, um, wm, vm, dvis, a, ro, GZ, Re, Fr, St, Fa, &
    t_end, TT, h0, cfreq)
    if (problem /= "dcavity" .AND. problem /= "backstep" .AND. problem /= "circle") then
        call WAVE_PARAMETERS (problem, Fr, St, GZ, TT, h0, cc, lambda)
    end if
    if (mesh_type=='uniform') then
        do i=0, imax+1
            dxx(i) = delx
        end do
        do j=0, jmax+1
            dzz(j) = delz
        end do
    else
        call GRID_GEN( problem, mesh_file, imax, jmax, xlen, zlen, dxx, dzz, h0, lambda, xtoe)
    end if
    !
    ! allocate arrays
    !
    allocate (UNP1(0:imax+1,0:jmax+1))
    allocate (WNP1(0:imax+1,0:jmax+1))
    allocate (UN(0:imax+1,0:jmax+1))
    allocate (WN(0:imax+1,0:jmax+1))
    allocate (UNM1(0:imax+1,0:jmax+1))
    allocate (WNM1(0:imax+1,0:jmax+1))
    allocate (F(0:imax+1,0:jmax+1))
    allocate (G(0:imax+1,0:jmax+1))
    allocate (P(0:imax+1,0:jmax+1))
    allocate (dpp(0:imax+1,0:jmax+1))
    allocate (RHS(0:imax+1,0:jmax+1))
    allocate (PS(0:imax+1,0:jmax+1))
    allocate (uflag(0:imax+1,0:jmax+1))
    allocate (vflag(0:imax+1,0:jmax+1))
    allocate (pflag(0:imax+1,0:jmax+1))
    if (VOLUM_TRACKING) then 
        allocate (C(0:imax+1,0:jmax+1))
    end if
    if (immersed) then         
        allocate (ua(imax,jmax))
        allocate (va(imax,jmax))
        allocate (ux(0:imax+1))
        allocate (uy(0:jmax+1))
        allocate (vx(0:imax+1))
        allocate (vy(0:jmax+1))
        allocate (xc(0:imax+1))
        allocate (yc(0:jmax+1))
        allocate (vos(0:imax+1,0:jmax+1))
    end if
    if (problem=="reservoir") then 
        allocate (EC(0:20000))
    end if
    allocate (dx(0:imax+1))
    allocate (dz(0:jmax+1))
    allocate (x(0:imax))
    allocate (z(0:jmax))
    allocate (H(0:imax+1))
    
    do i=0, imax+1
        dx(i) = dxx(i)
    end do

    do j=0, jmax+1
        dz(j) = dzz(j)
    end do

    deallocate(dxx, dzz)
    x(0) = 0.0
    do i = 1, imax
        x(i) = dx(i) + x(i-1)
    end do
    z(0) = 0.0
    do j = 1, jmax
        z(j) = dz(j) + z(j-1)
    end do
    t=0.0
    uxb=0.0
    uyb=0.0
    call INITFLAG (uflag, vflag, pflag, imax, jmax, ibound, vos)
    call WR_UWP(problem, imax, jmax, x, z, dx, dz, UN, WN, P, t)
    
    
    if (immersed) then
        !generate immersed boundary on Eulerian mesh
        call ARCGEN(problem, imax, jmax, dx, dz, xm, ym, sc, nm, t, h0, xlen, zlen, arclen, xtoe, dam_height, sw)
        allocate(nxm(nm))
        allocate(nym(nm))
        call IBM(problem, imax, jmax, ppu, ppv, ppp, ibound, t, x, z, dx, dz, nm, xm, ym, sc, nxm, nym, uflag, vflag, &
                pflag, ua, va, ivp, jvp, ux, vx, uy, vy, xc, yc, invagu, invagv, invagp, invaguin, invagvin, &
                gFIu, gFIv, gFIp, gFIuin, gFIvin, indu, indv, indp, induin, indvin, ghu, ghv, ghp, ghuin, &
                ghvin, IMu, IMv, IMp, arclen, vos)
        allocate(au(ppu))
        allocate(bu(ppu))
        allocate(cu(ppu))
        allocate(av(ppv))
        allocate(bv(ppv))
        allocate(cv(ppv))
        allocate(ap(ppp))
        allocate(bp(ppp))
        allocate(cp(ppp))
    end if
    call INIT_UWP(problem, UNM1, WNM1, P, x, z, dx, dz, imax, jmax, ro, GZ, h0, H, Ls, nm, xm, ym)  
    if(VOLUM_TRACKING) then
        call INIT_VOLUME(problem, C, pflag, imax, jmax, dx, dz, H)
    end if
    do i=0, imax+1
        do j=0, jmax+1
            dpp(i,j)=0.0
        end do
    end do
    call SETBCOND (problem, UNM1, WNM1, uflag, vflag, pflag, imax, jmax, ppu, ppv, wW, wE, wNN, wS,  indu, &
                    indv, induin, indvin, gFIu, gFIv, gFIuin, gFIvin, invagu, invagv, invaguin, invagvin, &
                    ghu, ghv, ghuin, ghvin, IMu, IMv, ivp, jvp, uxb, uyb)
    !call SETSPECBCOND (problem, dz, UNM1, WNM1, imax, jmax, um, cfreq, t, c, G, P, Fr, delt, q, St, PN)    
    do i=0, imax+1
        do j=0, jmax+1
            UN(i,j) = UNM1(i,j)
            WN(i,j) = WNM1(i,j)
            UNP1(i,j)=UN(i,j)
            WNP1(i,j)=WN(i,j)
            dpp(i,j)=0.0
        end do
    end do
    ts=0.0
    step=0
    zn=0
    tstr=0.0
    pi=4.0*atan(1.0)
    !if (problem =="reservoir")  call READ_EARTHQUAKE_ACCEL(EC, delt_rec)
    if (immersed.and.movbound) then 
        call WR_PSI(problem, imax, jmax, x, z, dz, UNP1, uflag, t)
        call WR_UWP(problem, imax, jmax, x, z, dx, dz, UN, WN, P, t)
        if (movbound) call WR_CIRCLE_POS(t, nm, xlen, zlen, zn)
    end if
    !call FIND_INTERS_OBS (problem, imax, jmax, nm, x, z, xm, ym, sc, nxm, nym, vos)
    !call PLT(UNP1, WNP1, P, C, pflag, imax, jmax,dx, dz, t)
    do while (t < t_end)
        if (problem=="solitary") call WR_SOLITARY_SHELF(imax, dx, H, t, h0)
        if (problem=="reservoir") then
            alpha=1.0
            axn = um*cfreq*cos(cfreq*(t+delt/2.0))
            !axn=-alpha*(-GZ)
            !axn = -um*cfreq*cos(cfreq*(t+delt/2.0))
            !axn = SEISMIC((t+delt/2.0), delt_rec, EC)*(-GZ)
        elseif (problem=="circle".and.movbound) then
            !axn=4.0*(5.0/2.0/pi)*pi**2*sin(2.0*pi*t)
            axn=2.0*pi*sin(2.0*pi*t)
        end if
        step=step+1
        call COMP_delt (problem, delt, imax, jmax, dx, dz, UN, WN, Re, St, cfl, H, cc, uxb, uyb, nm)
        delt=0.001
        !if (problem /= "dcavity" .AND. problem /= "backstep" .AND. problem /= "circle") then
        if (problem /= "dcavity" .AND. problem /= "backstep" ) then
            if (VOLUM_TRACKING) then
                call MARK_CELLS2 (C, pflag, imax, jmax, ifull)
				call SET_UWP_SURFACE_VARIABLE (UN, WN, PS, pflag, -axn, GZ, imax, jmax, Fr, Re, dx, dz, delt, c)
                call PRESSURE_INTERPOLATE (PS, P, C, pflag, imax, jmax, dx, dz)
                if (immersed) call HALLOW (problem, imax, jmax, pflag, vos, C)
            else
                call MARK_CELLS (pflag, imax, jmax, dz, ifull, H)
				call SET_UVP_SURFACE_TUMMAC (UN, WN, P, pflag, imax, jmax, Fr, Re, dx, dz)
			end if
        else
            ifull = imax*jmax - ibound
        end if
        
        call COMP_FG (problem, upwind, UN, WN, UNM1, WNM1, F, G, imax, jmax, gamma, &
                      delt, St, Re, Fr, dx, dz, GZ, axnm1, axn, P, ro, pflag, wE, wW, &
                      induin, indvin, indu, indv, gFIuin, gFIvin, gFIu, gFIv, invaguin, invagvin, invagu, &
                      invagv, ghuin, ghvin, ghu, ghv, ivp, jvp, ppu, ppv, IMu, IMv, uxb, uyb)
        if (immersed) then 
            call FORCING(imax, jmax, ppu, ppv, UN, WN, F, G, dx, dz, delt, t, indu, indv, gFIu, gFIv, &
                            invagu, invagv, ghu, ghv, IMu, IMv, au, bu, cu, av, bv, cv, uxb, uyb)
        end if
        call COMP_RHS (F, G, P, RHS, pflag, imax, jmax, delt, Fr, St, Fa, ro, dx, dz)
        if (ifull > 0) then
            call POISSON (problem, dpp, RHS, pflag, imax, jmax, ppp, dx, dz, z, eps, iter, &
                          itermax, omega, res, ifull, St, Fr, Fa, delt, H, c, ro, indp, gFIp, &
                          invagp, ghp, IMp, ap, bp, cp, q)
            print *,'t = ', t, 'res= ', res, 'iter= ', iter
        end if
        !advance the intermediate velocity
        call ADAP_UW (UNP1, WNP1, F, G, dpp, pflag, H, imax, jmax, dx, dz, z, Fr, St, ro, delt, c)
        do i=0, imax+1
            do j=0, jmax+1
                P(i,j)=P(i,j)+dpp(i,j)
                !dpp(i,j)=0.0
            end do
        end do
        
        CALL SETBCOND (problem, UNP1, WNP1, uflag, vflag, pflag, imax, jmax, ppu, ppv, wW, wE, wNN, wS,  &
                        indu, indv, induin, indvin, gFIu, gFIv, gFIuin, gFIvin, invagu, invagv, invaguin, &
                        invagvin, ghu, ghv, ghuin, ghvin, IMu, IMv, ivp, jvp, uxb, uyb)
        
        call SETSPECBCOND (problem, dz, UNP1, WNP1, imax, jmax, um, cfreq, t, c, G, dpp, Fr, delt, q, St)
         do i=0, imax+1
            do j=0, jmax+1
                dpp(i,j)=0.0
            end do
        end do
     
        if (problem /= "dcavity" .AND. problem /= "backstep" .AND. problem/="circle") then
            call SET_UWP_SURFACE_VARIABLE (UNP1, WNP1, PS, pflag, -axn, GZ, imax, jmax, Fr, Re, dx, dz, delt, c)
            if (VOLUM_TRACKING) then 
                call VOF_PLIC (c, UNP1, WNP1, imax, jmax, wW, wE, delt, dx, dz, pflag)
                call VOF_COMP_DEPTH (pflag, c, imax, jmax, dz, H)
            end if
        end if
        if (problem == "dcavity" .OR. problem == "backstep" .OR. problem=="circle") then
            if (COMP_RES (UNP1, WNP1, UN, WN, imax, jmax, ifull)) exit
        end if
        do i=0, imax+1
            do j=0, jmax+1
                UNM1(i,j)=UN(i,j)
                WNM1(i,j)=WN(i,j)
                UN(i,j)=UNP1(i,j)
                WN(i,j)=WNP1(i,j)
            end do
        end do
        if (t.ge.0.0.and.ts.ge.0.1) then
            if (immersed) then 
                call COMP_SURFACE_FORCE(problem, imax, jmax, nm, ppp, ppu, ppv, invagp, invagu, invagv, gFIp, gFIu, gFIv, &
                                        indp, indu, indv, ghp, ghu, ghv, P, UNP1, WNP1, sc, Re, xm, ym, t, au, bu, cu, av, &
                                        bv, cv, ap, bp, cp, uxb, uyb, ro, GZ, h0, alpha, sw, pflag, dx, dz, xc, yc)
            end if
            !call WR_PSI(problem, imax, jmax, x, z, dz, UNP1, uflag, t)
            !call WR_UWP(problem, imax, jmax, x, z, dx, dz, UNP1, WNP1, P, t)
            !call PLT(UNP1, WNP1, P, C, pflag, imax, jmax,dx, dz, t)
            call WR_FREESURFACE(imax, dx, H, t, h0)
            
            if (movbound) call WR_CIRCLE_POS(t, nm, xlen, zlen, zn)
            if (.not.(immersed).and.problem=="reservoir") then 
                call DAM_FACE(problem, pflag, imax, jmax, P, H, dz, t, h0, ro, GZ, alpha)
            end if
            ts=0.0
        end if
        axnm1=axn
        t = t + delt
        ts = ts + delt
        if (t.ge.200.0) then 
            tstr=tstr+delt
            if (tstr.ge.0.1) then
                call WR_PSI(problem, imax, jmax, x, z, dz, UNP1, uflag, t)
                tstr=0.0
            end if
        end if
    end do
    
    call WR_DCAVITY(imax, jmax, UNP1, P, z, dz, xlen, um, ro, Re)
    call WR_PSI(problem, imax, jmax, x, z, dz, UNP1, pflag, t)
    if (problem=="reservoir".and.immersed)  then
        !compute hydrodynamic pressure at a point in the reservoir
        allocate (ph(0:imax+1,0:jmax+1))
        do i=1, imax
            do j=1, jmax
                if (yc(j).le.h0) then 
                    ph(i,j)=(P(i,j)-ro*(-GZ)*(h0-yc(j)))/(ro*alpha*(-GZ)*h0**2)
                else
                    ph(i,j)=0.0
                end if
            end do
        end do
        call WR_UWP(problem, imax, jmax, x, z, dx, dz, UNP1, WNP1, ph, t)
    else
        call WR_UWP(problem, imax, jmax, x, z, dx, dz, UNP1, WNP1, P, t)
    end if                
    if (.not.(immersed).and.problem=="reservoir") then 
        call DAM_FACE(problem, pflag, imax, jmax, P, H, dz, t, h0, ro, GZ, alpha)
    end if
    if (immersed) then 
        call COMP_SURFACE_FORCE(problem, imax, jmax, nm, ppp, ppu, ppv, invagp, invagu, invagv, gFIp, gFIu, gFIv, &
                                indp, indu, indv, ghp, ghu, ghv, P, UNP1, WNP1, sc, Re, xm, ym, t, au, bu, cu, av, &
                                bv, cv, ap, bp, cp, uxb, uyb, ro, GZ, h0, alpha, sw, pflag, dx, dz, xc, yc)
        do j=1, jmax
            if (uy(j).ge.15.0d0) then 
                print *,'geldi'
                exit 
            end if
        end do
        do i=1, imax
            if (ux(i).ge.10.5d0) then 
                uu=UNP1(i,j)
                write (35,*) ux(i)-10.5d0, uu
            end if
        end do 
    end if
    

    stop ' normal'

    end program NASSLARDBM

