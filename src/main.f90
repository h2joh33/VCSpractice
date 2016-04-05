PROGRAM VCSpractice
<<<<<<< HEAD
    USE Util
    
=======
    USE Util    
>>>>>>> Implementation_of_Utility_Module
    IMPLICIT NONE
    
    TYPE(neutron) :: nst
    INTEGER :: i, icyc
    INTEGER :: nht=10000
    INTEGER :: ncyc=150, ninact=50
    INTEGER :: nact
<<<<<<< HEAD
    LOGICAL :: refl(4)
=======
    LOGICAL :: refl(4), term
    REAL(8) :: wtadj, dtc, dts, keff, kavg, kstd
    INTEGER :: surf, cotype
>>>>>>> Implementation_of_Utility_Module
    
    nact=ncyc-ninact
    refl = .TRUE.
    
<<<<<<< HEAD
=======
    CALL InitMCUtil(nht)
    
    DO icyc=1, ncyc
        wtadj=GetWtAdj(nht)
        DO i=1, GetNcurQue()
            CALL GenFS(nst, i)
            nst.wt=wtadj
            term=.FALSE.
            CALL SmpDir(nst.dir)
            DO WHILE(.not. term)
                dtc=SmpDTC(nst)
                CALL calDTS(nst, nst.dir, dts, surf)
                IF (dts .lt. dtc) THEN
                    CALL movetrck(nst.xy, dts, nst.dir, nst)
                    SELECT CASE(surf)
                    CASE(SOUTH)
                        nst.idxy=nst.idxy-1
                    CASE(EAST)
                        nst.idxx=nst.idxx+1
                    CASE(WEST)
                        nst.idxx=nst.idxx-1
                    CASE(NORTH)
                        nst.idxy=nst.idxy+1
                    END SELECT
                    IF (chkboundary(nst.idxx, nst.idxy)) THEN
                        IF (refl(surf)) THEN
                            SELECT CASE(surf)
                            CASE(SOUTH)
                                nst.idxy=nst.idxy+1
                                nst.dir.y=nst.dir.y*-1
                            case(EAST)
                                nst.idxx=nst.idxx-1
                                nst.dir.x=nst.dir.x*-1
                            case(WEST)
                                nst.idxx=nst.idxx+1
                                nst.dir.x=nst.dir.x*-1
                            case(NORTH)
                                nst.idxy=nst.idxy-1
                                nst.dir.y=nst.dir.y*-1
                            END SELECT
                        ELSE
                            term = .TRUE.
                        END IF
                    END IF
                ELSE
                    CALL movetrck(nst.xy, dtc, nst.dir, nst)
                    cotype=Collision(nst)
                    IF (cotype .ne. CTABS) THEN
                        term = .FALSE.
                        CALL SmpDir(nst.dir)
                    ELSE
                        term = .TRUE.
                    END IF
                END IF
                IF (term) CYCLE
            END DO
        END DO
    
        keff=EstimateKeff(nht)
        IF (icyc>ninact) THEN
            CALL AccTally()
        END IF
        CALL PrepareNxtCycle
        PRINT *, icyc, keff
    END DO
    CALL GetAvgTally(kavg, kstd)
    PRINT *, "kavg",kavg,"std=",kstd
>>>>>>> Implementation_of_Utility_Module
    
    
    
END PROGRAM VCSpractice