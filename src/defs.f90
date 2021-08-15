    module defs
    use nrtype
    INTEGER(I2B), parameter :: C_B = 0
    INTEGER(I2B), parameter :: B_N = 1   ! IBSET(C_B,0)
    INTEGER(I2B), parameter :: B_S = 2   ! IBSET(C_B,1)
    INTEGER(I2B), parameter :: B_W = 4   ! IBSET(C_B,2)
    INTEGER(I2B), parameter :: B_E = 8   ! IBSET(C_B,3)
    
    INTEGER(I2B), parameter :: B_NW = B_N + B_W
    INTEGER(I2B), parameter :: B_SW = B_S + B_W
    INTEGER(I2B), parameter :: B_NE = B_N + B_E
    INTEGER(I2B), parameter :: B_SE = B_S + B_E
    INTEGER(I2B), parameter :: B_WE = B_W + B_E
    INTEGER(I2B), parameter :: B_NS = B_N + B_S
    
    INTEGER(I2B), parameter :: B_SWE=B_S + B_W + B_E       
    INTEGER(I2B), parameter :: B_NSW=B_N + B_S + B_W        
    INTEGER(I2B), parameter :: B_NWE=B_N + B_W + B_E        
    INTEGER(I2B), parameter :: B_NSE=B_N + B_S + B_E     
    
    INTEGER(I2B), parameter :: B_NSWE=B_N + B_S + B_W + B_E

    INTEGER(I2B), parameter :: C_F = 16 ! IBSET(C_B,4) C_F=0010 (fluid)
    INTEGER(I2B), parameter :: C_X = C_F - 1
    INTEGER(I2B), parameter :: C_A = C_F + B_N + B_S + B_E + B_W
    !
    !  for free-boundary flows
    !
    INTEGER(I2B), parameter :: C_O = 256  ! IBSET(C_B,8)  C_O=0100
    INTEGER(I2B), parameter :: C_W = 512  ! IBSET(C_B,9)  C_W=0200
    INTEGER(I2B), parameter :: C_S = 1024 ! IBSET(C_B,10) C_S=0400
    INTEGER(I2B), parameter :: C_N = 2048 ! IBSET(C_B,11) C_N=0800
    !
    INTEGER(I2B), parameter :: C_E = 4096 ! IBSET(C_B,12) C_E=1000 (empty)
    !
    INTEGER(I2B), parameter :: C_WO= C_W + C_O              !C_WO=0300
    INTEGER(I2B), parameter :: C_NS= C_N + C_S              !C_NS=0c00
    INTEGER(I2B), parameter :: C_SW= C_S + C_W              !C_SW=0600
    INTEGER(I2B), parameter :: C_NW= C_N + C_W              !C_NW=0a00
    INTEGER(I2B), parameter :: C_NO= C_N + C_O              !C_NO=0900
    INTEGER(I2B), parameter :: C_SO= C_S + C_O              !C_SO=0500
    !
    INTEGER(I2B), parameter :: C_SWO=C_S + C_W + C_O        !C_SWO=0700
    INTEGER(I2B), parameter :: C_NSW=C_N + C_S + C_W        !C_NSW=0e00
    INTEGER(I2B), parameter :: C_NWO=C_N + C_W + C_O        !C_NWO=0b00
    INTEGER(I2B), parameter :: C_NSO=C_N + C_S + C_O        !C_NSO=0d00
    !
    INTEGER(I2B), parameter :: C_NSWO=C_N + C_S + C_W + C_O !C_NSWO=0f00
    
  
    logical VOLUM_TRACKING
    logical immersed
    logical movbound
    logical hf
    logical curvature
    logical absorption
    integer swp_dir
    real(DP) emf
    real(DP) emf1
    end module defs

