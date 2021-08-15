    module interfaces
    !
    ! init.f90 interfaces
    !

    interface
    subroutine READ_PARAMETERS (FILE, problem, atype, mesh_file, meshtype, upwind, &
    xlength, zlength, h0, TT, um, wm, imax, jmax, &
    t_end , cfl, del_out, gamma, f_c, d_c, &
    delx, delz, &
    itermax, eps, omega, ro, dvis, a, GZ, &
    wW, wE, wN, wS)
    use nrtype
    implicit none
    integer, INTENT(OUT) :: imax, jmax, &
    itermax, &
    wW, wE, wN, wS, f_c, d_c
    real(RP), INTENT(OUT) :: xlength, zlength, h0, TT, um ,wm, &
    t_end , cfl, del_out, gamma, &
    eps, omega, ro, dvis, a, GZ, delx, delz
    character (LEN=12), INTENT(IN) :: FILE
    character (LEN=30), INTENT(OUT) :: problem, atype, mesh_file, meshtype, upwind
    end subroutine READ_PARAMETERS
    end interface

    interface
    subroutine INIT_UWP(problem, U, W, P, x, z, dx, dz, imax, jmax, ro, GZ, h0, H, Ls, nm, xm, ym)
    use nrtype
    implicit none
    character (LEN=30), INTENT(IN) :: problem
    integer, INTENT(IN) :: imax, jmax, nm
    real(RP), INTENT(IN) :: ro, GZ, h0
    REAL(RP), INTENT(OUT) :: Ls
    real(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, W, P
    real(RP), DIMENSION(0:), INTENT(IN) :: x, z, dx, dz
    real(RP), DIMENSION(0:), INTENT(INOUT) :: H
    real(DP), dimension(:), allocatable, intent(IN) :: xm, ym
    end subroutine INIT_UWP
    end interface

    interface
    subroutine INITFLAG (uflag, vflag, pflag, imax, jmax, ibound, vos)
    use nrtype
    use defs
    implicit none
    integer, INTENT(IN) :: imax, jmax
    integer(I2B), dimension(0:,0:), intent(INOUT) :: uflag, vflag, pflag
    integer, intent(OUT) :: ibound
    real(DP), dimension(:,:), intent(INOUT) :: vos
    end subroutine INITFLAG
    end interface

    interface
    subroutine DIMENSIONLESS (problem, atype, imax, jmax, xlen, zlen, delx, delz, um, wm, vm, &
                                dvis, a, ro, GZ, Re, Fr, St, Fa, &
                                t_end , TT, h0, cfreq)
    use nrtype
    implicit none
    character (LEN=30), INTENT(IN) :: problem, atype
    integer, INTENT(IN) :: imax, jmax
    real(RP), INTENT(IN) :: dvis, a
    real(RP), INTENT(INOUT) :: GZ, xlen, zlen, delx, delz, um, wm, vm, ro, h0
    real(RP), INTENT(OUT) :: Re, Fr, St, Fa, t_end , TT, cfreq
    end subroutine DIMENSIONLESS
    end interface

    interface
    subroutine SETBCOND (problem, U, W, uflag, vflag, pflag, imax, jmax, ppu, ppv, wW, wE, wN, wS, &
                          indu, indv, induin, indvin, gFIu, gFIv, gFIuin, gFIvin, invagu, invagv, invaguin, &
                          invagvin, ghu, ghv, ghuin, ghvin, IMu, IMv, ivp, jvp, uxb, uyb)
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
    end subroutine SETBCOND
    end interface
    
    
    interface
    subroutine SETSPECBCOND (problem, dz, U, W, imax, jmax, um, cfreq, t, c, G, P, Fr, delt, q, St)
    use nrtype
    implicit none
    character (LEN=30) :: problem
    real(RP), dimension(0:), intent(IN) :: dz
    real(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, W
    integer, INTENT(IN) :: imax, jmax
    real(RP), intent(IN) :: um, cfreq, t, Fr, delt, q, St
    real(DP), dimension(0:,0:), intent(IN) :: c
    real(RP), dimension(0:,0:), intent(IN) :: G, P
    end subroutine SETSPECBCOND
    end interface

    interface
    subroutine COMP_delt (problem, delt, imax, jmax, dx, dz, U, W, Re, St, cfl, H, c, uxb, uyb, nm)
    use nrtype
    implicit none
    character (LEN=30) :: problem
    integer, INTENT(IN) :: imax, jmax, nm
    real(RP), INTENT(IN) :: Re, St, cfl, c, uxb, uyb
    real(RP), INTENT(OUT) :: delt
    real(RP), DIMENSION(0:,0:), INTENT(IN) :: U, W
    real(RP), DIMENSION(0:), INTENT(IN) :: dx, dz, H
    end subroutine COMP_delt
    end interface
    !
    ! momentumu.f90 interfaces
    !

    interface
    subroutine COMP_FG (problem, upwind, U, W, UNM1, WNM1, F, G, imax, jmax, gamma, &
    delt, St, Re, Fr, dx, dz, GZ, axnm1, axn, P, ro, pflag, wE, wW, &
    induin, indvin, indu, indv, gFIuin, gFIvin, gFIu, gFIv, invaguin, invagvin, invagu, &
    invagv, ghuin, ghvin, ghu, ghv, ivp, jvp, ppu, ppv, IMu, IMv, uxb, uyb)
    use nrtype
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
    end subroutine COMP_FG
    end interface

    interface
    subroutine COMP_RHS (F, G, P, RHS, pflag, imax, jmax, delt, Fr, St, Fa, ro, dx, dz)
    use nrtype
    implicit none
    integer, INTENT(IN) :: imax, jmax
    real(RP), INTENT(IN) :: delt, Fr, St, Fa, ro
    real(RP), DIMENSION(0:), INTENT(IN) :: dx, dz
    integer(I2B), DIMENSION(0:,0:), INTENT(IN) :: pflag
    real(RP), DIMENSION(0:,0:), INTENT(IN) :: F, G, P
    real(RP), DIMENSION(0:,0:), INTENT(OUT) :: RHS
    end subroutine COMP_RHS
    end interface

    interface
    subroutine POISSON (problem, P, RHS, pflag, imax, jmax, ppp, dx, dz, z, eps, iter, &
    itermax, omega, res, ifull, St, Fr, Fa, delt, H, c, ro, indp, gFIp, invagp, ghp, IMp, ap, bp, cp, q)
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
    end subroutine POISSON   
    end interface   

    interface
    subroutine ADAP_UW (U, W, F, G, P, pflag, H, imax, jmax, dx, dz, z, Fr, St, ro, delt, c)
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
    end subroutine ADAP_UW
    end interface

    interface
    LOGICAL FUNCTION COMP_RES (U, W, UT, WT, imax, jmax, ifull)
    use nrtype
    implicit none
    real(RP), DIMENSION(0:,0:), INTENT(IN) :: U, W, UT, WT
    integer, INTENT(IN) :: imax, jmax, ifull
    end FUNCTION COMP_RES
    end interface

    interface
    subroutine WR_DCAVITY(imax, jmax, U, P, z, dz, xlen, um, ro, Re)
    use nrtype
    implicit none
    integer, INTENT(IN):: imax, jmax
    real(RP), INTENT(IN):: xlen, um, ro, Re
    real(RP), DIMENSION(0:), INTENT(IN):: z, dz
    real(RP), DIMENSION(0:,0:), INTENT(IN):: U, P
    end subroutine WR_DCAVITY
    end interface

    interface
    subroutine WR_PSI(problem, imax, jmax, x, z, dz, U, FLAG, t)
    use nrtype
    implicit none
    character (LEN=30), INTENT(IN) :: problem
    integer, INTENT(IN):: imax, jmax
    real(RP), intent(IN):: t
    integer(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
    real(RP), DIMENSION(0:),INTENT(IN):: x, z, dz
    real(RP), DIMENSION(0:,0:),INTENT(IN):: U
    end subroutine WR_PSI
    end interface

    interface
    subroutine WAVE_PARAMETERS (problem, Fr, St, GZ, TT, h0, c, lambda)
    use nrtype
    implicit none
    character (LEN=30) :: problem
    real(RP), INTENT(IN) :: Fr, St, GZ, TT, h0
    real(RP), INTENT(OUT) :: c, lambda
    end subroutine
    end interface



    interface
    subroutine MARK_CELLS (FLAG, imax, jmax, dz, ifull, H)
    use nrtype
    implicit none
    integer, INTENT(IN) :: imax, jmax
    integer, INTENT(OUT) :: ifull
    real(RP), DIMENSION(0:), INTENT(IN):: dz
    integer(I2B), DIMENSION(0:,0:), INTENT(INOUT) :: FLAG
    real(RP), DIMENSION(0:), INTENT(IN) :: H
    end subroutine MARK_CELLS
    end interface

    interface
    subroutine SET_UVP_SURFACE_TUMMAC (U, W, P, FLAG, imax, jmax, Fr, Re, dx, dz)
    use nrtype
    implicit none
    integer, INTENT(IN) :: imax, jmax
    real(RP), INTENT(IN) :: Fr, Re
    integer(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
    real(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, W, P
    real(RP), DIMENSION(0:), INTENT(IN) :: dx, dz
    end subroutine SET_UVP_SURFACE_TUMMAC
    end interface

    interface
    subroutine DIC(imax, jmax, FLAG, dx, dz, U,  H, dt, St)
    use nrtype
    implicit none
    integer, INTENT(IN):: imax, jmax
    integer(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
    real(RP), DIMENSION(0:), INTENT(IN):: dx, dz
    real(RP), DIMENSION(0:), INTENT(INOUT):: H
    real(RP), DIMENSION(0:,0:), INTENT(IN):: U
    real(RP), INTENT(IN):: dt, St
    end subroutine DIC
    end interface

    

   

    interface
    subroutine WR_FREESURFACE(imax, dx, H, t, hi)
    use nrtype
    implicit none
    integer, INTENT(IN):: imax
    real(RP), INTENT(IN):: t, hi
    real(RP), DIMENSION(0:), INTENT(IN):: dx, H
    end subroutine WR_FREESURFACE
    end interface
    INTERFACE
      SUBROUTINE WR_FREESURFACE3(problem, FLAG, C, imax, jmax, dx, dz, &
                                    t, h0, vxb, vyb, ny, Ls)
	     USE nrtype
	     IMPLICIT NONE
        CHARACTER (LEN=30) :: problem
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
	     INTEGER, INTENT(IN):: imax, jmax, ny
	     REAL(RP), INTENT(IN):: t, h0, Ls
	     REAL(RP), DIMENSION(0:), INTENT(IN):: dx, dz
	     REAL(DP), DIMENSION(0:,0:), INTENT(IN) :: C
         real(DP), dimension(:), intent(IN) :: vxb, vyb
	  END SUBROUTINE WR_FREESURFACE3
   END INTERFACE

    interface
    subroutine WR_RES2(FLAG, imax, jmax, H, dz, t, h0)
    use nrtype
    use defs
    implicit none
    integer, INTENT(IN):: imax, jmax
    integer(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
    real(RP), DIMENSION(0:), INTENT(IN):: H, dz
    real(RP), INTENT(IN):: t, h0
    end subroutine WR_RES2
    end interface
    
    interface
        subroutine WR_RES(FLAG, jmax, H, dz, C, t, h0)
       USE nrtype
    USE defs
    implicit none
    integer, intent(IN):: jmax
    integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
    real(RP), dimension(0:), intent(IN):: H, dz
    real(DP), dimension(0:,0:), intent(IN) :: C
    real(RP), intent(IN):: t, h0
    end subroutine WR_RES
    end interface


 INTERFACE
      SUBROUTINE INIT_VOLUME(problem, C, FLAG, imax,jmax, dx, dz, H)
        USE nrtype
        IMPLICIT NONE
        character (LEN=30), INTENT(IN) :: problem
        INTEGER, INTENT(IN) :: imax, jmax
		REAL(RP), DIMENSION(0:), INTENT(IN):: dx, dz, H
		INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
        REAL(DP), DIMENSION(0:,0:), INTENT(INOUT) :: C
      END SUBROUTINE INIT_VOLUME
   END INTERFACE
INTERFACE
      SUBROUTINE MARK_CELLS2 (C, FLAG, imax, jmax, ifull)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imax, jmax
         INTEGER, INTENT(OUT) :: ifull
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(INOUT) :: FLAG
    	 REAL(DP), DIMENSION(0:,0:), INTENT(INOUT) :: C
      END SUBROUTINE MARK_CELLS2
   END INTERFACE

 INTERFACE
      SUBROUTINE PLT(U, W, P, C, FLAG, imax, jmax, dx, dz, tt)
        USE nrtype
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: imax, jmax
		REAL(RP), INTENT(IN) :: tt
		INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
		REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, W, P
		REAL(DP), DIMENSION(0:,0:), INTENT(IN) :: C
		REAL(RP), DIMENSION(0:), INTENT(IN):: dx, dz
	  END SUBROUTINE PLT
   END INTERFACE
INTERFACE
      SUBROUTINE SET_UWP_SURFACE_VARIABLE (U, W, PS, FLAG, GX, GZ, imax, jmax, &
                               Fr, Re, dx, dz, delt, c)
      USE nrtype
	  USE defs
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imax, jmax
      REAL(RP), INTENT(IN) :: GX, GZ, Fr, Re, delt
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, W, PS
      REAL(RP), DIMENSION(0:), INTENT(IN):: dx, dz
      real(DP), dimension(0:,0:), intent(IN) :: c
      END SUBROUTINE SET_UWP_SURFACE_VARIABLE
   END INTERFACE


       

     
     interface
        subroutine VOF_PLIC (c, U, W, imax, jmax, wW, wE, delt, dx, dz, FLAG)
      use nrtype
	  use defs
      implicit none
      integer, intent(IN) :: imax, jmax, wW, wE
      real(RP), intent(IN) :: delt
      real(RP), dimension(0:), intent (IN):: dx, dz
      real(RP), dimension(0:,0:), intent(IN) :: U, W
      real(DP), dimension(0:,0:), intent(INOUT) :: c
      integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
      end subroutine VOF_PLIC
     end interface
     

     
     interface
        subroutine VOF_COMP_DEPTH (FLAG, c, imax, jmax, dz, H)
            use nrtype
	        use defs
            implicit none
            integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
            real(DP), dimension(0:,0:), intent(IN) :: c
            integer, intent(IN) :: imax, jmax
            real(RP), dimension(0:), intent (IN):: dz
            real(RP), dimension(0:), intent(OUT) :: H
            end subroutine VOF_COMP_DEPTH
     end interface
     
     INTERFACE
        SUBROUTINE PRESSURE_INTERPOLATE (PS, P, c, FLAG, imax, jmax, dx, dz)
            USE nrtype
	        USE defs
             IMPLICIT NONE
            INTEGER, INTENT(IN) :: imax, jmax
            REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: PS
            REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: P
            REAL(DP), DIMENSION(0:,0:), INTENT(IN) :: c
            INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
            REAL(RP), DIMENSION(0:), INTENT(IN):: dx, dz
       END SUBROUTINE PRESSURE_INTERPOLATE
   END INTERFACE
       
     interface
     subroutine IBM(problem, imax, jmax, ppu, ppv, ppp, ibound, t, xg, yg, dx, dy, N, x, y, sc, nxm, nym, uflag, vflag, &
                    pflag, ua, va, ivp, jvp, ux, vx, uy, vy, xc, yc, &
                    invagu, invagv, invagp, invaguin, invagvin, gFIu, gFIv, gFIp, gFIuin, gFIvin, &
                    indu, indv, indp, induin, indvin, ghu, ghv, ghp, ghuin, ghvin, IMu, IMv, IMp, P, vos)
    use defs
    implicit none
    character (LEN=30), INTENT(IN) :: problem
    integer, intent(IN) :: imax, jmax, N
    real(RP), intent(IN):: t
    real(RP), dimension(0:), intent(IN) :: xg, yg, dx, dy    
    real(DP), dimension(:), intent(IN) :: x, y, sc
    real(DP), dimension(:), intent(OUT) :: nxm, nym
    integer, intent(OUT) :: ibound, ppu, ppv, ppp, ivp, jvp
    integer(I2B), dimension(0:,0:), intent(INOUT) :: uflag, vflag, pflag
    real(DP), dimension(:,:,:), allocatable, intent(OUT) :: invagu, invagv, invagp, invaguin, invagvin
    integer, dimension(:,:,:), allocatable, intent(OUT) :: gFIu, gFIv, gFIp, gFIuin, gFIvin
    integer, dimension(:,:), allocatable, intent(OUT) :: indu, indv, indp, induin, indvin
    real(DP), dimension(:,:), allocatable, intent(OUT) :: ghu, ghv, ghp, ghuin, ghvin
    integer, dimension(:), allocatable, intent(OUT) :: IMu, IMv, IMp
    integer, dimension(:,:), intent(OUT) :: ua, va
    real(DP), dimension(0:), intent(OUT) :: ux, vx, uy, vy, xc, yc
    real(DP), intent(IN) :: P
    real(DP), dimension(0:,0:), intent(INOUT) :: vos
     end subroutine IBM
     end interface
     
     
     
     interface
        subroutine FIND_INTERS_OBS (problem, imax, jmax, nm, x, z, xm, ym, sm, nxm, nym, vos)
        use nrtype
        use interfaces2
        use defs
        implicit none
        character (LEN=30), INTENT(IN) :: problem
        integer, intent(IN):: imax, jmax, nm
        real(RP), dimension(0:), intent(IN):: x, z
        real(DP), dimension(:), intent(IN) :: xm, ym, sm
        real(DP), dimension(:), intent(IN) :: nxm, nym
        real(DP), dimension(0:,0:), intent(INOUT) :: vos
        end subroutine
     end interface
     
     interface
    subroutine COMP_SURFACE_FORCE(problem, imax, jmax, nm, ppp, ppu, ppv, invagp, invagu, invagv, &
                gFIp, gFIu, gFIv, indp, indu, indv, ghp, ghu, ghv, P, U, V, s, Re, xm, ym, t, au, bu, cu, &
                av, bv, cv, ap, bp, cp, uxb, uyb, &
                ro, GZ, h0, alphag, sw, pflag, dx, dy, xc, yc)
    use defs
    implicit none
    character (LEN=30), intent(IN) :: problem
    integer, INTENT(IN) :: imax, jmax, nm, ppp, ppu, ppv
    real(DP), dimension(:,:,:), intent(IN) :: invagp, invagu, invagv
    integer, dimension(:,:,:), intent(IN) :: gFIp, gFIu, gFIv
    integer, dimension(:,:), intent(IN) :: indp, indu, indv
    real(DP), dimension(:,:), intent(INOUT) :: ghp, ghu, ghv
    real(RP), dimension(0:,0:), intent(IN) :: P, U, V
    real(DP), dimension(:), intent(IN) :: s, xm, ym
    real(RP), intent(IN):: Re, t, uxb, uyb, ro, GZ, h0, alphag, sw
    real(DP), dimension(:), intent(IN) :: au, bu, cu, av, bv, cv, ap, bp, cp
    integer(I2B), dimension(0:,0:), intent(IN) :: pflag
    real(RP), dimension(0:), intent(IN) :: dx, dy
    real(DP), dimension(0:), intent(IN) :: xc, yc
     end subroutine COMP_SURFACE_FORCE
     end interface
     
    ! INTERFACE
    !    subroutine ADAP_UW (imax, jmax, U, W, UINT, WINT, delt)
    !   use nrtype
    !use defs
    !implicit none
    !integer, INTENT(IN) :: imax, jmax
    !real(RP), INTENT(IN) :: delt
    !real(RP), DIMENSION(0:,0:), INTENT(IN) :: UINT, WINT
    !real(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, W
    !END subroutine  ADAP_UW
    !END INTERFACE
    
   
    
    
    
   interface
    subroutine ARCGEN(problem, Mx, My, dx, dy, x, y, sc, N, t, h0, xlen, zlen, P, xtoe, dam_height, sw)
    use nrtype
    use defs
    implicit none
    character (LEN=30), INTENT(IN) :: problem
    integer, intent(IN) :: Mx, My
    real(RP), dimension(0:), intent(IN) :: dx, dy 
    real(DP), dimension(:), allocatable, intent(OUT) :: x, y, sc
    integer, intent(OUT) :: N
    real(RP), intent(IN) :: t, h0, xlen, zlen, xtoe
    real(DP), intent(OUT) :: P
    real(RP), intent(OUT) :: dam_height, sw
    end subroutine ARCGEN
    end interface
    
    interface
    subroutine FORCING(imax, jmax, ppu, ppv, U, V, F, G, dx, dy, delt, t, indu, indv, gFIu, gFIv, invagu, invagv, &
                        ghu, ghv, IMu, IMv, au, bu, cu, av, bv, cv, uxb, uyb)
    use nrtype
    use defs
    implicit none
    integer, intent(IN) :: imax, jmax, ppu, ppv     
    real(RP), dimension(0:,0:), intent(INOUT) :: U, V
    real(RP), dimension(0:,0:), intent(INOUT) :: F, G
    real(RP), dimension(0:), intent(IN) :: dx, dy
    real(RP), intent(IN):: delt, t, uxb, uyb
    integer, dimension(:,:), intent(IN) :: indu, indv
    integer, dimension(:,:,:), intent(IN) :: gFIu, gFIv
    real(DP), dimension(:,:,:), intent(IN) :: invagu, invagv
    real(DP), dimension(:,:), intent(IN) :: ghu, ghv
    integer, dimension(:), intent(IN) :: IMu, IMv
    real(DP), dimension(:), intent(OUT) :: au, bu, cu, av, bv, cv
    end subroutine FORCING
    end interface
    
    

    interface
    subroutine WR_CIRCLE_POS(t, N, xlen, zlen, step)
    use nrtype
    use defs
    implicit none
    integer, intent(IN)::  N
    integer, intent(INOUT):: step
    real(RP), intent(IN):: t, xlen, zlen
    end subroutine WR_CIRCLE_POS
    end interface
    
    interface
        integer function FINDPOINT(x, vp)
       use nrtype
    use defs
    implicit none
    real(DP), INTENT(IN):: vp
    real(DP), DIMENSION(0:), INTENT(IN):: x
    end function FINDPOINT 
    end interface
    
    interface
    subroutine GRID_GEN( problem, mesh_file, imax, jmax, xlen, zlen, dx, dz, h0, lambda, xtoe)
	use nrtype
	implicit none
	character (LEN=30), intent(IN) :: problem, mesh_file
	integer, intent (OUT) :: imax, jmax
	real(RP), intent(IN):: h0, lambda
    real(RP), intent(INOUT):: xlen, zlen
	real(RP), dimension(0:), intent(OUT):: dx, dz
    real(RP), intent(OUT):: xtoe
    end subroutine GRID_GEN
    end interface
    
    interface
    subroutine COMP_DRAGF (problem, upwind, U, W, UNM1, WNM1, F, G, &
                            imax, jmax, gamma, delt, St, Re, Fr, dx, dz, GZ, P, ro, &
                            uflag, vflag, x, z, t, xm, ym, nm, nx, ny)
    use nrtype
    use defs
    implicit none
    character (LEN=30), INTENT(IN) :: problem, upwind
    integer, INTENT(IN) :: imax, jmax, nm
    real(RP), INTENT(IN) :: gamma, delt, St, Re, Fr, GZ, ro, t
    real(RP), DIMENSION(0:), INTENT(IN) :: dx, dz, x, z
    real(RP), DIMENSION(0:,0:), INTENT(IN) :: U, W, UNM1, WNM1, P, F, G
    integer, dimension(0:,0:), intent(IN) :: uflag, vflag
    real(DP), dimension(:), intent(IN) :: xm, ym
    real(DP), dimension(:), intent(IN) :: nx, ny
    end subroutine COMP_DRAGF
    end interface
    
    INTERFACE
	  SUBROUTINE WR_UWP(problem, imax, jmax, x, z, dx, dz, U, W, P, t)
		 USE nrtype
		 IMPLICIT NONE
         character (LEN=30), INTENT(IN) :: problem
		 INTEGER, INTENT(IN):: imax, jmax
		 REAL(RP), INTENT(IN):: t
		 REAL(RP), DIMENSION(0:),INTENT(IN):: x, z, dx, dz
		 REAL(RP), DIMENSION(0:,0:),INTENT(IN):: U, W, P
	  END SUBROUTINE WR_UWP
      END INTERFACE
      
       INTERFACE
      SUBROUTINE HIRT_VOF (C, U, W, FLAG, dx, dz, imax, jmax, dt)
      USE nrtype
	  use defs
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imax, jmax
      REAL(RP), INTENT(IN) :: dt
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: u, w
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:), INTENT(IN) :: dx, dz
      REAL(DP), DIMENSION(0:,0:), INTENT(INOUT) :: c
     END SUBROUTINE HIRT_VOF
     END INTERFACE
     
     INTERFACE

 	SUBROUTINE swpxz (c,u,w,flag,flxm,flxp,flzm,flzp,divx,divz,imax,kmax, delt,dell)
      USE nrtype
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imax, kmax
      REAL(RP), INTENT(IN) :: delt, dell
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: u,w
      REAL(DP), DIMENSION(0:,0:), INTENT(INOUT) :: c,flxm,flxp,flzm,flzp,divx,divz
	  END SUBROUTINE swpxz
      END INTERFACE 
      
       
       
       INTERFACE
   SUBROUTINE SWPzz (C, V, FLAG, VOF1, VOF2, VOF3, imax, kmax, delt,dell)
      USE nrtype
      USE defs
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imax, kmax
      REAL(RP), INTENT(IN) :: delt, dell
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: V
      REAL(Dp), DIMENSION(0:,0:), INTENT(INOUT) :: C, VOF1, VOF2, VOF3
	  END SUBROUTINE SWPzz
   END INTERFACE

 INTERFACE
   SUBROUTINE SWPxx (C, U, FLAG, VOF1, VOF2, VOF3, imax, kmax, delt,dell)
      USE nrtype
      USE defs
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imax, kmax
      REAL(RP), INTENT(IN) :: delt, dell
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U
      REAL(Dp), DIMENSION(0:,0:), INTENT(INOUT) :: C, VOF1, VOF2, VOF3
	  END SUBROUTINE SWPxx
   END INTERFACE

INTERFACE
  SUBROUTINE YOUNG (c0, u, v, imax, jmax, dt, dxx, dyy, wE)
      USE nrtype
	  use defs
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imax, jmax, wE
      REAL(RP), INTENT(IN) :: dt
      REAL(RP), DIMENSION(0:), INTENT (IN):: dxx, dyy
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: u, v
      REAL(DP), DIMENSION(0:,0:), INTENT(INOUT) :: c0
  END SUBROUTINE YOUNG
END INTERFACE 

INTERFACE
    SUBROUTINE transport (c, itype, tanalfa, tanbeta, u1, u3, u4, u2, dx, dy, dt, f1, f3, f4, f2)
      USE nrtype
	  use defs
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: itype
      REAL(DP), INTENT(IN):: c
      REAL(RP), INTENT(IN) :: tanalfa, tanbeta, u1, u3, u4, u2, dx, dy, dt
      REAL(RP), INTENT(OUT) :: f1, f3, f4, f2
  END SUBROUTINE transport
  END INTERFACE      
  
  interface
    subroutine SET_UVP_SURFACE_VOF (U, W, P, c, FLAG, flagib, imax, jmax, Fr, Re, GX, GZ, delt, dx, dz)
        use nrtype
        use defs
        implicit none
        integer, INTENT(IN) :: imax, jmax
        real(RP), INTENT(IN) :: Fr, Re, GX, GZ, delt
        integer(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG, flagib
        real(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, W, P
        real(DP), dimension(0:,0:), intent(INOUT) :: c
        real(RP), DIMENSION(0:), INTENT(IN) :: dx, dz
    end subroutine SET_UVP_SURFACE_VOF
  end interface
  
  interface
      subroutine COMP_NUT (FLAG, NUT, KA, EP, imax, jmax)
      use nrtype
      use defs
      implicit none
	  integer, intent(IN) :: imax, jmax
	  integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
	  real(RP), dimension(0:,0:), intent(IN) :: KA, EP
	  real(RP), dimension(0:,0:), intent(OUT) :: NUT
    end subroutine COMP_NUT
  end interface
  
  interface
    subroutine COMP_KAEP (FLAG, U, W, NUT, KANP1, EPNP1, KAN, EPN, imax, jmax, delt, dx, dz, x, z, &
                          gamma, Re, St)
   use nrtype
   use defs
   implicit none
   integer, intent(IN) :: imax, jmax
   integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
   real(RP), dimension(0:,0:), intent(IN) :: U, W, NUT, KAN, EPN
   real(RP), dimension(0:), intent(IN) :: dx, dz, x, z
   real(RP), dimension(0:,0:), intent(INOUT) :: KANP1, EPNP1
   real(RP), intent(IN) :: delt, gamma, Re, St
   end subroutine COMP_KAEP
  end interface
  
  !interface
  !  real FUNCTION A_X(i, j, P, dx)
	 ! use nrtype
	 ! implicit none
	 ! integer, intent(IN) :: i, j
	 ! real(RP), dimension(0:,0:), intent(IN) :: P
	 ! real(RP), dimension(0:), intent(IN) :: dx
  !  end function A_X
  !  end interface
  !  
  !  interface
  !      real function B_X(i, j, P, dx)
	 !       use nrtype
	 !       implicit none
	 !       integer, intent(IN) :: i, j
	 !       real(RP), dimension(0:,0:), intent(IN) :: P
	 !       real(RP), dimension(0:), intent(IN) :: dx
  !      end function B_X
  !  end interface
  !  
  !  interface
  !      real function A_Z(i, j, P, dz)
	 !       use nrtype
	 !       implicit none
	 !       integer, intent(IN) :: i, j
	 !       real(RP), dimension(0:,0:), intent(IN) :: P
	 !       real(RP), dimension(0:), intent(IN) :: dz
  !      end function A_Z
  !  end interface
  !  
  !  interface
  !      real function B_Z(i, j, P, dz)
	 !       use nrtype
	 !       implicit none
	 !       integer, intent(IN) :: i, j
	 !       real(RP), dimension(0:,0:), intent(IN) :: P
	 !       real(RP), dimension(0:), intent(IN) :: dz
  !      end function B_Z
  !  end interface
  !  
    interface
        real function CONV_FI(i, j, gamma, U, W, dx, dz, FI)
            use nrtype
            use defs
            IMPLICIT NONE
            integer, INTENT (IN):: i, j
            real(RP), INTENT(IN):: gamma
            real(RP), DIMENSION(0:), INTENT (IN):: dx, dz
            real(RP), DIMENSION(0:,0:), INTENT(IN):: U, W, FI
        end function CONV_FI
    end interface
    
    interface
        real function DIFF_FI(i, j, U, W, dx, dz, FI, NUT, sgm, Re)
            use nrtype
            use defs
            IMPLICIT NONE
            integer, INTENT (IN):: i, j
            real(RP):: sgm, Re
            real(RP), DIMENSION(0:), INTENT (IN):: dx, dz
            real(RP), DIMENSION(0:,0:), INTENT(IN):: U, W, FI, NUT
        end function DIFF_FI
    end interface
    
    interface
        real function GRADU(i, j, U, W, dx, dz)
            use nrtype
            use defs
            IMPLICIT NONE
            integer, INTENT (IN):: i, j
            real(RP), DIMENSION(0:), INTENT (IN):: dx, dz
            real(RP), DIMENSION(0:,0:), INTENT(IN):: U, W
        end function GRADU
    end interface
    
    interface
        real function TDiff_U(i, j, U, W, NUT, KA, dx, dz, Re)
            use nrtype
            IMPLICIT NONE
            integer, INTENT(IN):: i, j
            real(RP):: Re
            real(RP), DIMENSION(0:), INTENT(IN):: dx, dz
            real(RP), DIMENSION(0:, 0:), INTENT(IN):: U, W, NUT, KA
        end function TDiff_U
    end interface
    
    interface
        real function TDiff_W(i, j, U, W, NUT, KA, dx, dz, Re)
            use nrtype
            IMPLICIT NONE
            integer, INTENT(IN):: i, j
            real(RP):: Re
            real(RP), DIMENSION(0:), INTENT(IN):: dx, dz
            real(RP), DIMENSION(0:, 0:), INTENT(IN):: U, W, NUT, KA
        end function TDiff_W
    end interface
    
    interface
        real function Wall(y, u, Re, uto)
    	    use nrtype
	        implicit none
	        real(RP), intent(IN):: y, u, Re, uto
        end function Wall
    end interface
    
    interface
        real function dWall(y, u, Re, uto)
    	    use nrtype
	        implicit none
	        real(RP), intent(IN):: y, u, Re, uto
        end function dWall
    end interface
    
    interface
        subroutine WR_TURBULENT_QUANT(imax, jmax, U, KA, EP, x, z, dx, dz, xlen, h0, Re, ro)
            use nrtype
	        implicit none
	        integer, intent(IN):: imax, jmax
	        real(RP), intent(IN):: xlen, h0, Re, ro
	        real(RP), dimension(0:), intent(IN):: x, z, dx, dz
	        real(RP), dimension(0:,0:), intent(IN):: U, KA, EP
        end subroutine WR_TURBULENT_QUANT
    end interface
    
    interface
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
        end subroutine NEIGHBORINGPRESSURE
    end interface
    
    interface
        subroutine WR_DAM_BASE(imax, jmax, t, ro, GZ, P, h0, dz, H, alpha)
            use nrtype
            implicit none
            integer, intent(IN):: imax, jmax
            real(RP), intent(IN):: t, ro, GZ, h0, alpha
            real(RP), dimension(0:,0:), intent(IN):: P
            real(RP), dimension(0:), intent(IN):: dz, H
        end subroutine WR_DAM_BASE
    end interface
    
    interface
    
    subroutine DAM_FACE(problem, FLAG, imax, jmax, P, H, dz, t, h0, ro, GZ, alpha)
        use nrtype
        use defs
	    character (LEN=30), INTENT(IN) :: problem
        integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
	    integer, intent(IN):: imax, jmax
	    real(RP), intent(IN):: t, h0, ro, GZ, alpha
	    real(RP), dimension(0:), intent(IN):: dz, H
	    real(RP), dimension(0:,0:), intent(IN):: P
    end subroutine DAM_FACE
    end interface
    
    interface
        subroutine WR_SOLITARY_SHELF(imax, dx, H, t, h0)
            use nrtype
            implicit none
            integer, intent(IN):: imax
            real(RP), intent(IN):: t, h0
            real(RP), dimension(0:), intent(IN):: dx, H
        end subroutine WR_SOLITARY_SHELF
    end interface
    
    
    interface
        subroutine CHOPRA (jmax, h0, a, t, T0, alpha, ro, GZ, z)
            use nrtype
            implicit none
            integer, intent(IN) :: jmax
            real(RP), intent(IN) :: h0, a, t, T0, alpha, ro, GZ
            real(RP), dimension(0:), INTENT(IN) :: z
        end subroutine CHOPRA
    end interface
    
    !interface
    !integer function FINDFLUIDX(flag, i, j, const)
    !use nrtype
    !use defs
    !implicit none
    !integer, intent(IN) :: i, j, const
    !integer(I2B), dimension(0:,0:), intent(IN) :: flag
    !end function FINDFLUIDX
    !end interface
    !
    !interface
    !integer function FINDFLUIDY(flag, i, j, const)
    !use nrtype
    !use defs
    !implicit none
    !integer, intent(IN) :: i, j, const
    !integer(I2B), dimension(0:,0:), intent(IN) :: flag
    !end function FINDFLUIDY
    !end interface
    
    
    
    interface
	  subroutine READ_EARTHQUAKE_ACCEL(EC, delt)
	     use nrtype
	     implicit none
	     real(RP), dimension(0:), intent(OUT):: EC
	     real(RP), intent(OUT):: delt
	  end subroutine READ_EARTHQUAKE_ACCEL
    end interface
    
interface
	  real function SEISMIC(t, dt, EA)
	     use nrtype
	     implicit none
	     real(RP), intent(IN):: t, dt
	     real(RP), dimension(0:), intent(IN):: EA
	  end function SEISMIC
end interface

     interface
        subroutine HALLOW (problem, imax, jmax, FLAG, vos, C)
        use nrtype
        use defs
        implicit none
        character (LEN=30), INTENT(IN) :: problem
        integer, intent(IN):: imax, jmax
        integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
        real(DP), dimension(0:,0:), intent(INOUT) :: vos, C
        end subroutine HALLOW
     end interface



!interface
!    real(DP) function fxcir(x, y)
!        use nrtype
!        implicit none
!        real(DP), intent(IN):: x, y
!    end function
!end interface
    
    end module interfaces



