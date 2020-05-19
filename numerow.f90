    MODULE NUMEROW
    USE constants
    IMPLICIT NONE

    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Applying Numerow method to solve
    !    d^2 Y(x)
    !    _______  =  - G(x)*Y(X) + S(X)
    !      d x^2
    !
    ! by integration from inside ( X(1) ) to outside ( X(NR) )
    ! ASSUMING PROPER INITIAL VALUES AT THE BOUNDARY!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE  INTEGRATE_Y_USING_NUMEROW_FROM_INSIDE(Y,G,S,NR,DR)
    IMPLICIT NONE
    INTEGER :: NR ! number of components in functions Y,G,S
    REAL(p) :: Y(NR),G(NR),S(NR) ! functions
    REAL(p) :: DR ! step size
    ! local
    INTEGER :: I
    
    DO I = 2, NR - 1
        Y(I + 1) = (2d0 * Y(I) * (1d0 - 5d0/12d0*DR**2 * G(I)) - Y(I - 1) * &
            (1d0 + G(I-1) * DR**2 / 12d0) + DR**2/12d0 * &
            (S(I+1) + 10d0 * S(I) + S(I - 1))) / (1d0 + DR**2 * G(I+1) / 12d0)
    END DO

    END SUBROUTINE INTEGRATE_Y_USING_NUMEROW_FROM_INSIDE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Applying Numerow method to solve
    !    d^2 Y(x)
    !    _______  =  - G(x)*Y(X) + S(X)
    !      d x^2
    !
    ! by integration from outside ( X(NR) ) to inside ( X(1) )
    ! ASSUMING PROPER INITIAL VALUES AT THE BOUNDARY!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    SUBROUTINE  INTEGRATE_Y_USING_NUMEROW_FROM_OUTSIDE(Y,G,S,NR,DR)
    IMPLICIT NONE
    INTEGER :: NR ! number of components in functions Y,G,S
    REAL(p) :: Y(NR),G(NR),S(NR) ! functions
    REAL(p) :: DR ! step size
    ! local
    INTEGER :: I

    DO I = NR - 1, 2
        Y(I - 1) = (2d0 * Y(I) * (1d0 - 5d0/12d0*DR**2 * G(I)) - Y(I + 1) * &
            (1d0 + G(I-1) * DR**2 / 12d0) + DR**2/12d0 * &
            (S(I+1) + 10d0 * S(I) + S(I - 1))) / (1d0 + DR**2 * G(I-1) / 12d0)
    END DO

    END SUBROUTINE INTEGRATE_Y_USING_NUMEROW_FROM_OUTSIDE


    END MODULE
