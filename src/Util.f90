MODULE Util
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: SOUTH=1, WEST=2, NORTH=3, EAST=4
    REAL(8) :: pi=3.141592
    
    TYPE direction
        REAL(8) :: x, y, z
    END TYPE
    
    TYPE neutron
        REAL(8) :: xyz(3)
        INTEGER :: g
        INTEGER :: idxx, idxy, idxz
        REAL(8) :: wt
        TYPE(direction) :: dir
    END TYPE
    
    TYPE(fsource), POINTER, DIMENSION(:) :: que, quenxt, quetmp
    
    INTEGER :: nq, nqnxt
    INTEGER :: nx, ny
    REAL(8) :: gk=0
    REAL(8), PRIVATE :: keff
    REAL(8) :: ksum=0., ksqsum=0.
    INTEGER, PARAMETER :: CTABS=1, CTSCT=2
    INTEGER, PARAMETER :: nofis=1, yesfis=2
    INTEGER :: nacc=0
    
    
    
    
    FUNCTION chkboundary(ix, iy)
        INTEGER, INTENT(IN) :: ix, iy
        LOGICAL :: chkboundary
        chkboundary = .FALSE.
        IF (ix .lt. 1 .or. ix .gt. nx) chkboundary = .TRUE.
        IF (iy .lt. 1 .or. iy .gt. ny) chkboundary = .TRUE.
    END FUNCTION
    
    
        
    
END MODULE Util