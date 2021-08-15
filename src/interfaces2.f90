MODULE interfaces2
    INTERFACE
     SUBROUTINE ASSIGNINTERNALBOUND(problem, imax, jmax, flag)
     USE nrtype
     USE defs
     IMPLICIT NONE
     CHARACTER (LEN=30), INTENT(IN) :: problem
     INTEGER, INTENT(IN) :: imax, jmax
     INTEGER, DIMENSION(0:,0:), INTENT(INOUT) :: flag
     END SUBROUTINE ASSIGNINTERNALBOUND
     END INTERFACE
     
     interface
     subroutine ASSIGNINTERNALVOS(problem, imax, jmax, vos)
        use nrtype
        use defs
        implicit none
        integer, intent(IN) :: imax, jmax
        real(DP), dimension(0:,0:), intent(INOUT) :: vos
        character (LEN=30), INTENT(IN) :: problem
     end subroutine ASSIGNINTERNALVOS
     end interface
     
    INTERFACE
    INTEGER FUNCTION FINDFLUIDX(imax, ac, i, j, const)
    USE nrtype
    USE defs
    IMPLICIT NONE
    INTEGER, intent(IN) :: imax, i, j, const
    INTEGER, dimension(:,:), intent(IN) :: ac
    END FUNCTION FINDFLUIDX
    END INTERFACE

    INTERFACE
    INTEGER FUNCTION FINDFLUIDY(jmax, ac, i, j, const)
    USE nrtype
    USE defs
    IMPLICIT NONE
    INTEGER, intent(IN) :: jmax, i, j, const
    INTEGER, dimension(:,:), intent(IN) :: ac
    END FUNCTION FINDFLUIDY
    END INTERFACE

    INTERFACE
    subroutine MATRIXINVERSA(pp, r, a)
     USE nrtype
     USE defs
     IMPLICIT NONE
    integer, intent(IN) :: pp
    real(DP), dimension(:,:), intent(IN) :: r
    real(DP), dimension(:,:,:), intent(OUT) :: a
    end subroutine MATRIXINVERSA
    END INTERFACE

    INTERFACE
    subroutine MATRIXINVERSB(pp, r, a)
    use nrtype
    implicit none
    integer, intent(IN) :: pp
    real(DP), dimension(:,:), intent(IN) :: r
    real(DP), dimension(:,:,:), intent(OUT) :: a
    end subroutine MATRIXINVERSB
    END INTERFACE

    INTERFACE
    INTEGER FUNCTION J_OF_Z(jmax, h, dz)
    USE nrtype
    IMPLICIT NONE
    INTEGER, INTENT(IN):: jmax
    REAL(RP), INTENT(IN):: h
    REAL(RP), DIMENSION(0:), INTENT(IN):: dz
    END FUNCTION  J_OF_Z
    END INTERFACE

    interface
    real FUNCTION CONV_UQ(i, j, imax, jmax, U, W, dx, dz)
    use nrtype
    implicit none
    integer, INTENT (IN):: i, j, imax, jmax
    real(RP), DIMENSION(0:), INTENT (IN):: dx, dz
    real(RP), DIMENSION(0:,0:), INTENT(IN):: U, W
    end FUNCTION CONV_UQ
    end interface

    interface
    real FUNCTION CONV_U(i, j, U, W, dx, dz, gamma)
    use nrtype
    implicit none
    integer, INTENT (IN):: i, j
    real(RP), DIMENSION(0:), INTENT (IN):: dx, dz
    real(RP), DIMENSION(0:,0:), INTENT(IN):: U, W
    real(RP), INTENT(IN)::gamma
    end FUNCTION CONV_U
    end interface

    interface
    real FUNCTION DIFF_U(i, j, U, dx, dz)
    use nrtype
    implicit none
    integer, INTENT(IN):: i, j
    real(RP), DIMENSION(0:), INTENT(IN):: dx, dz
    real(RP), DIMENSION(0:, 0:), INTENT(IN):: U
    end FUNCTION DIFF_U
    end interface

    interface
    real FUNCTION COMP_X(i, j, U, W, dx, dz)
    use nrtype
    implicit none
    integer, INTENT(IN):: i, j
    real(RP), DIMENSION(0:), INTENT(IN):: dx, dz
    real(RP), DIMENSION(0:, 0:), INTENT(IN):: U, W
    end FUNCTION COMP_X
    end interface
    !
    ! momentumw.f90 interfaces
    !
    interface
    real FUNCTION CONV_WQ(i, j, imax, jmax, U, W, dx, dz)
    use nrtype
    implicit none
    integer, INTENT (IN):: i, j, imax, jmax
    real(RP), DIMENSION(0:), INTENT (IN):: dx, dz
    real(RP), DIMENSION(0:, 0:), INTENT (IN)::U, W
    end FUNCTION CONV_WQ
    end interface

    interface
    real FUNCTION CONV_W(i, j, U, W, dx, dz, gamma)
    use nrtype
    implicit none
    integer, INTENT (IN):: i, j
    real(RP), DIMENSION(0:), INTENT (IN):: dx, dz
    real(RP), DIMENSION(0:, 0:), INTENT (IN)::U, W
    real(RP), INTENT(IN)::gamma
    end FUNCTION CONV_W
    end interface

    interface
    real FUNCTION DIFF_W(i, j, W, dx, dz)
    use nrtype
    implicit none
    integer, INTENT(IN):: i, j
    real(RP), DIMENSION(0:), INTENT(IN):: dx, dz
    real(RP), DIMENSION(0:,0:), INTENT(IN):: W
    end FUNCTION DIFF_W
    end interface

    interface
    real FUNCTION SRC_W(i, j, GZ, dx, dz)
    use nrtype
    implicit none
    integer, INTENT(IN)::i, j
    real(RP), INTENT(IN):: GZ
    real(RP), DIMENSION(0:), INTENT(IN):: dx, dz
    end FUNCTION SRC_W
    end interface

    interface
    real FUNCTION COMP_Z(i, j, U, W, dx, dz)
    use nrtype
    implicit none
    integer, INTENT(IN):: i, j
    real(RP), DIMENSION(0:), INTENT(IN):: dx, dz
    real(RP), DIMENSION(0:,0:), INTENT(IN):: U, W
    end FUNCTION COMP_Z
    end interface

    INTERFACE
     SUBROUTINE swpx (c,u,w,imax,kmax, wW, wE, delt,dx, dz)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imax, kmax, wW, wE
         REAL(RP), INTENT(IN) :: delt
         REAL(RP), DIMENSION(0:), INTENT (IN):: dx, dz
         REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: u,w
         REAL(DP), DIMENSION(0:,0:), INTENT(INOUT) :: c
         END SUBROUTINE swpx
         end interface

         INTERFACE
     SUBROUTINE swpz (c,u,w,imax,kmax, wE, delt,dx, dz)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imax, kmax, wE
         REAL(RP), INTENT(IN) :: delt
         REAL(RP), DIMENSION(0:), INTENT (IN):: dx, dz
         REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: u,w
         REAL(DP), DIMENSION(0:,0:), INTENT(INOUT) :: c
     END SUBROUTINE swpz

     END INTERFACE
      interface
        subroutine VOFBOUND (c, FLAG, dx, dz, imax, jmax)
        use nrtype
        use defs
        implicit none
        integer, intent(IN) :: imax, jmax
        real(RP), dimension(0:), intent (IN):: dx, dz
        real(DP), dimension(0:,0:), intent(INOUT) :: c
        integer(I2B), dimension(0:,0:), intent(IN) :: FLAG
        end subroutine VOFBOUND
      end interface
      
      interface
    subroutine RegulaFalsi(fnf, a, b, roots, iFind, iFound, sec)
        use nrtype
        implicit none
        real(DP), intent(IN):: a, b, sec
        integer, intent(IN):: iFind
        integer, intent(OUT):: iFound
        real(DP), dimension(:), allocatable, intent(OUT) :: roots
        real(DP), external:: fnf
    end subroutine RegulaFalsi
    end interface
    
    !interface
    !subroutine FCT(imax, jmax, H, cfl, dx)
    !use nrtype
    !implicit none
    !integer, INTENT(IN):: imax, jmax
    !real(RP), INTENT(IN):: cfl
    !real(RP), DIMENSION(0:), INTENT(INOUT):: H
    !real(RP), DIMENSION(0:), INTENT(IN):: dx
    !end subroutine FCT
    !end interface
    
     interface
    real FUNCTION SGN(x)
    use nrtype
    implicit none
    real(RP), INTENT(IN)::  x
    end FUNCTION SGN
    end interface


END MODULE interfaces2
