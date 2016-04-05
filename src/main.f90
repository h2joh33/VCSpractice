PROGRAM VCSpractice
    USE Util
    
    IMPLICIT NONE
    
    TYPE(neutron) :: nst
    INTEGER :: i, icyc
    INTEGER :: nht=10000
    INTEGER :: ncyc=150, ninact=50
    INTEGER :: nact
    LOGICAL :: refl(4)
    
    nact=ncyc-ninact
    refl = .TRUE.
    
    
    
    
END PROGRAM VCSpractice