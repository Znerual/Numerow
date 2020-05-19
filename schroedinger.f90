    MODULE SCHROEDINGER
    USE constants
    IMPLICIT NONE

    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Main module that solves the radial Schroedinger equation
    !
    !
    !    d^2 U(x)
    !    ___________  =  -2*(E-Z/r+POT) U(x)
    !      d x^2
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE SOLVE_SE(EMIN,EMAX,Z,POT,R,NR,DR,U)
    USE numerow
    IMPLICIT NONE
    INTEGER, intent(in) :: NR
    REAL(p), intent(in) :: DR
    REAL(p), intent(inout) :: EMIN,EMAX,Z
    REAL(p),  intent(in) :: POT(NR)
    REAL(p),  intent(inout) :: U(NR)
    REAL(p),  intent(in) :: R(NR)
    ! local
    REAL(p)  :: G(NR)
    REAL(p)  :: S(NR) ! dummy S=0 for numerow
    REAL(p) :: E ! Trial energy
    INTEGER :: NODES,I


    S=0.0_p
    call INIT_U_FOR_SE(U,R,Z,NR)
    
    DO WHILE(abs(EMAX - EMIN) > tol)   
        E=(EMAX+EMIN)/2.0_p       
        call CALC_G(E,Z,POT,R,NR,G)
        call INTEGRATE_Y_USING_NUMEROW_FROM_INSIDE(U,G,S,NR,DR)
        IF (CALC_NODES(U,NR) == 0) THEN
            EMIN = E
        ELSE
            EMAX = E
        END IF
    END DO
    
    
    !###########################################################
    !   TODO: CONSTRUCT A WHILE LOOP THAT FINDS THE SOLUTION TO THE SCHROEDINGER EQUATION:

    !         CALL THE REQUIRED ROUTINES IN THE CORRECT ORDER
    !         TO OBTAIN THE SOLUTION FOR THE SCHROEDINGER EQUATION FOR EACH TRIAL ENERGY
    !         (about 6 lines)
    !
    !         ONCE THE WAVEFUNCTIOM HAS BEEN OBTAINED COMPUTE ITS NUMBER OF NODES AND
    !         MODIFY EMIN and EMAX ACCORDINGLY. COMPUTE NEW TRIAL ENERGY
    !         (about 4 lines lines)
    !###########################################################


    CALL NORMALIZE_U(U,R,NR,DR)

    END SUBROUTINE SOLVE_SE


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !   Calculate function g(x) for radial Schroedinger equation in
    !   the way it is needed for the numerow method (depends on Z, E and \phi(r) ) )
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CALC_G(E,Z,POT,R,NR,G)
    IMPLICIT NONE
    REAL(p), intent(in) :: E,Z
    REAL(p), intent(in) :: R(NR),POT(NR)
    INTEGER, intent(in) :: NR
    REAL(p), intent(out) :: G(NR)
    ! local
    INTEGER :: I
    
    DO I = 1, NR
        G(I) = 2d0 * (E + Z/R(I) - POT(I))
    END DO

    END SUBROUTINE CALC_G

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  Set proper initial values of U at the boundary
    !  as needed by the Numerow method to solve the Schroedinger equation.
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE INIT_U_FOR_SE(U,R,Z,NR)
    IMPLICIT NONE
    REAL(p), DIMENSION(1), intent(inout) :: U(NR)
    REAL(p), DIMENSION(1), intent(in) :: R(NR)
    INTEGER  :: NR
    REAL(p)  :: Z

    U(1)=R(1)*(1-Z*R(1))
    U(2)=R(2)*(1-Z*R(2))

    END SUBROUTINE INIT_U_FOR_SE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  Normalize U properly
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE NORMALIZE_U(U,R,NR,DR)
    IMPLICIT NONE
    REAL(p), intent(inout) :: U(NR)
    REAL(p), intent(in) :: R(NR)
    INTEGER, intent(in) :: NR
    REAL(p), intent(in) :: DR
    ! local
    REAL(p) :: NORM
    INTEGER :: I

    NORM=0.0_p
    
    !TRAPEZREGEL
    NORM = NORM + U(1)**2
    DO I = 2, NR - 1
        NORM = NORM + 2d0 * U(I)**2
    END DO
    NORM = NORM + U(NR)**2
    
    NORM = sqrt(NORM * DR * 0.5d0)
    
    U = U / NORM
    !###########################################################
    !   TODO: Complete this routine to normalize U() (approximately 3 lines)
    !###########################################################


    END SUBROUTINE NORMALIZE_U

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  This function computes the number of nodes of the radial wavefunction U
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    FUNCTION CALC_NODES(U,NR)
    IMPLICIT NONE
    INTEGER :: CALC_NODES
    INTEGER, intent(in) :: NR
    REAL(p), intent(in) :: U(NR)
    ! local
    INTEGER :: I

    CALC_NODES=0
    
    DO I = 2, NR
        IF (U(I-1) * U(I) < 0d0) CALC_NODES = CALC_NODES + 1
    END DO
    !###########################################################
    !   TODO: Complete this function to compute the number of nodes of the function U() (approximately 3 lines)
    !###########################################################


    RETURN
    END FUNCTION CALC_NODES


    END MODULE
