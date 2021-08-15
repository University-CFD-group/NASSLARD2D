    module nrtype
    integer, parameter :: I1B = SELECTED_INT_KIND(2)
    integer, parameter :: I2B = SELECTED_INT_KIND(4)
    integer, parameter :: I4B = SELECTED_INT_KIND(9)
    integer, parameter :: SP = KIND(1.0)
    integer, parameter :: DP = KIND(1.d0)
    integer, parameter :: RP = SP
    end module nrtype
