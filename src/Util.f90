MODULE Util
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: SOUTH=1, WEST=2, NORTH=3, EAST=4
    REAL(8) :: pi=3.141592
    
    TYPE direction
        REAL(8) :: x, y
    END TYPE
    
    TYPE neutron
        REAL(8) :: xy(2)
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
    REAL(8) :: Lx=51, Ly=51, h=1
    TYPE xsec
        REAL(8), POINTER, DIMENSION(:) :: tot, nuf, abs
        REAL(8), POINTER, DIMENSION(:,:) :: sct, sctcdf
    END TYPE
    TYPE(xsec), POINTER, DIMENSION(:) :: xseclist
    ALLOCATE(xseclist(1))
    ALLOCATE(xseclist(1).tot(2)); ALLOCATE(xseclist(1).nuf(2)); ALLOCATE(xseclist(1).abs(2)); 
    ALLOCATE(xseclist(1).sct(2,2)); ALLOCATE(xseclist(1).sctcdf(2,2)); 
    xseclist(1).tot(1)=1.19623E-01; xseclist(1).tot(2)=1.44581E-01
    xseclist(1).nuf(1)=7.24142E-02; xseclist(1).nuf(2)=3.95781E-02
    xseclist(1).abs(1)=1.25593E-02; xseclist(1).abs(2)=1.54202E-02
    xseclist(1).sct(1,1)=3.27683E-02; xseclist(1).sct(1,2)=1.87221E-02
    xseclist(1).sct(2,1)=0.00000E-00; xseclist(2).sct(2,2)=5.51891E-02
    xseclist(1).sctcdf(1,1)=xseclist(1)sct(1,1)/xseclist(1).sct(1,2); xseclist(1).sctcdf(1,2)=xseclist(1)sct(1,2)/xseclist(1).sct(1,2); 
    xseclist(1).sctcdf(2,1)=xseclist(1)sct(2,1)/xseclist(1).sct(2,2); xseclist(1).sctcdf(2,2)=xseclist(1)sct(2,2)/xseclist(1).sct(2,2); 
    
    FUNCTION chkboundary(ix, iy)
        INTEGER, INTENT(IN) :: ix, iy
        LOGICAL :: chkboundary
        chkboundary = .FALSE.
        IF (ix .lt. 1 .or. ix .gt. nx) chkboundary = .TRUE.
        IF (iy .lt. 1 .or. iy .gt. ny) chkboundary = .TRUE.
    END FUNCTION
    
    SUBROUTINE InitMCUtil(nht)
        INTEGER :: ng
        INTEGER, INTENT(IN) :: nht
        TYPE(neutron) :: ns
        INTEGER :: i, idx(2)
        REAL(8) :: pos(2)
        
        nx=Lx/h
        ny=Ly/h
        
        ng=getNG()
        ALLOCATE(que(nht*2), quenxt(nht*2))
        nq=0
        DO i=1, nht
            CALL SmpPos(pos, idx)
            que(i).idxx=idx(1)
            que(i).idxy=idx(2)
            que(i).xy=pos
            nq=nq+1
        END DO
        nqnxt=0
        keff=1.
        ALLOCATE(compmap(nx,ny))
        compmap=1
    END SUBROUTINE
    
    SUBROUTINE GenFS(ns, i)
        TYPE(neutron), INTENT(INOUT) :: ns
        INTEGER :: i
        ns.xy=que(i).xy
        ns.idxx=que(i).idxx
        ns.idxy=que(i).idxy
        ns.g=1
    END SUBROUTINE
    
    SUBROUTINE SmpFN(ns)
        TYPE(neutron), INTENT(IN) :: ns
        INTEGER :: nfn, i, tcomp
        REAL(8) :: rn, xstot, xsnuf
        tcomp=compmap(ns.idxx, ns.idxy)
        xstot=xseclist(tcomp).tot(ns.g)
        xsnuf=xseclist(tcomp).nuf(ns.g)
        
        CALL RANDOM_NUMBER(rn)
        nfn=INTEGER(xsnuf/xstot/keff*ns.wt+rn)
        IF (nfn .gt. 0 ) THEN
            DO i=1, nfn
                nqnxt=nqnxt+1
                quenxt(nqnxt).xy=ns.xy
                quenxt(nqnxt).idxx=ns.idxx
                quenxt(nqnxt).idxy=ns.idxy
            END DO
        END IF
        
    END SUBROUTINE
    
    SUBROUTINE SmpPos(pos, idx)
        REAL(8) :: rn1, rn2
        REAL(8), INTENT(OUT) :: pos(2)
        INTEGER, INTENT(OUT) :: idx(2)
        CALL RANDOM_NUMBER(rn1)
        pos(1)=rn1*Lx
        CALL RANDOM_NUMBER(rn2)
        pos(2)=rn2*Ly
        idx(1)=pos(1)/h+1
        idx(2)=pos(2)/h+1
    END SUBROUTINE
    
    SUBROUTINE SmpDir(dir)
        REAL(8) :: rn1, rn2
        REAL(8) :: mu, ms, alpha
        REAL(8) :: pi=3.141592
        TYPE(direction), INTENT(OUT) :: dir
        CALL RANDOM_NUMBER(rn1)
        CALL RANDOM_NUMBER(rn2)
        mu=2*rn1-1
        alpha=2*pi*rn2
        ms=sqrt(1._8-mu**2)
        dir.x=ms*cos(alpha)
        dir.y=ms*sin(alpha)
        dir.z=mu
    END SUBROUTINE
    
    FUNCTION SmpDTC(ns) RESULT(r)
        TYPE(neutron), INTENT(IN) :: ns
        INTEGER :: tcomp
        REAL(8) :: xstot
        REAL(8) :: rn
        REAL(8) :: r
        tcomp=compmap(nx,ny)
        xstot=xseclist(tcomp).tot(ns.g)
        CALL RANDOM_NUMBER(rn)
        r=-DLOG(rn)/xstot
    END FUNCTION
    
    SUBROUTINE calDTS(ns, dir, dts, surf)
        TYPE(neutron), INTENT(IN) :: ns
        TYPE(direction), INTENT(IN) :: dir
        REAL(8), INTENT(OUT) :: dts
        INTEGER, INTENT(OUT) :: surf
        REAL(8) :: dx, dy, lx, ly
        IF (dir.x .gt. 0) THEN
            dx=ns.idxx*h-ns.xy(1)
        ELSE
            dx=(ns.idxx-1)*h-ns.xy(1)
        END IF
        IF (dir.y .gt. 0 )THEN
            dy=ns.idxy*h-ns.xy(2)
        ELSE 
            dy=(ns.idxy-1)*h-ns.xy(2)
        ENDIF
        
        lx=dx/dir.x
        ly=dy/dir.y
        
        IF (lx .lt. ly) THEN
            dts=lx
            IF (dir.x .gt. 0) THEN
                surf=EAST
            ELSE
                surf=WEST
            END IF
        ELSE
            dts=ly
            IF (dir.y .gt. 0 ) THEN
                surf=NOTRH
            ELSE
                surf=SOUTH
            END IF
        END IF
    END SUBROUTINE
    
    FUNCTION cmpdt(ns, dtc, dts) RESULT(dt)
        TYPE(neutron) :: ns
        REAL(8), INTENT(IN) :: dtc, dts
        REAL(8) :: dt
        IF (dtc .gt. dts) THEN
            dt=dts
        ELSE 
            dt=dtc
        END IF
    END FUNCTION
    
    SUBROUTINE movetrck(pos, length, dir, ns)
        REAL(8), INTENT(INOUT) :: pos(2)
        TYPE(neutron), INTENT(IN) :: ns
        REAL(8), INTENT(IN) :: length
        TYPE(direction), INTENT(IN) :: dir      
        
        pos(1)=pos(1)+length*dir.x
        pos(2)=pox(2)+length*dir.y
        
    END SUBROUTINE
    
    FUNCTION collision(ns) RESULT(cotype)
        TYPE(neutron), INTENT(INOUT) :: ns
        REAL(8) :: rn
        LOGICAL :: cotype
        REAL(8) :: lock, xsabs, xstot, xsnuf
        INTEGER :: idx, tcomp, gnxt
        
        tcomp=compmap(ns.idxx, ns.idxy)
        xsnuf=xseclist(tcomp).nuf(ns.g)
        xstot=xseclist(tcomp).tot(ns.g)
        xsabs=xseclist(tcomp).abs(ns.g)
        lock=xsnuf/xstot*ns.wt
        gk=gk+lock
        
        CALL SmpFN(ns)
        
        CALL RANDOM_NUMBER(rn)
        IF (rn .lt. xsabs/xstot) THEN
            cotype=CTABS
        ELSE
            cotype=CTSCT
        END IF
        
        IF (cotype .eq. CTSCT) THEN
            CALL RANDOM_NUMBER(rn)
            DO gnxt=1, getNG()-1
                IF (rn<xseclist(tcomp).sctcdf(ns.g, gnxt)) EXIT
            END DO
            ns.g=gnxt
        END IF
    END FUNCTION
    
    FUNCTION EstimateKeff(nht) RESULT(retkeff)
        INTEGER, INTENT(IN) :: nht
        REAL(8) :: retkeff
        keff=gk/nht
        retkeff=keff
        gk=0.
    END FUNCTION
    
    SUBROUTINE PrepareNxtCycle
        quetmp=>que
        que=>quenxt
        quenxt=>quetmp
        nq=nqnxt
        nqnxt=0
    END SUBROUTINE
    
    SUBROUTINE AccTally
        nacc=nacc+1
        ksum=ksum+keff
        ksqsum=ksqsum+keff**2
    END SUBROUTINE
    
    SUBROUTINE GetAvgTally(kavg, kstd)
        REAL(8), INTENT(OUT) :: kavg, kstd
        kavg=0.
        kstd=0.
        kavg=ksum/nacc
        kstd=sqrt((ksqsum-ksum**2/nacc)/(nacc-1)/nacc)
    END SUBROUTINE
        
    
END MODULE Util