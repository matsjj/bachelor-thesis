!Diverse Mathematische Funktionen, umgeschrieben auf F90
!Quelle: http://jin.ece.illinois.edu/routines/routines.html

!
!       ==================================================
!       Purpose: Compute complete elliptic integrals K(k)
!                and E(k)
!       Input  : K  --- Modulus k ( 0 ó k ó 1 )
!       Output : CK --- K(k)
!                CE --- E(k)
!       ==================================================
!

SUBROUTINE COMELP(HK,CK)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PK=1.0D0-HK*HK
        IF (HK.EQ.1.0) THEN
           CK=1.0D+300
           CE=1.0D0
        ELSE
           AK=(((.01451196212D0*PK+.03742563713D0)*PK+.03590092383D0)*PK+.09666344259D0)*PK+1.38629436112D0
           BK=(((.00441787012D0*PK+.03328355346D0)*PK+.06880248576D0)*PK+.12498593597D0)*PK+.5D0
           CK=AK-BK*DLOG(PK)
        ENDIF
        RETURN
        END


SUBROUTINE JELP(U,HK,ESN,ECN,EDN,EPH)
!
!       ========================================================
!       Purpose: Compute Jacobian elliptic functions sn u, cn u
!                and dn u
!       Input  : u   --- Argument of Jacobian elliptic fuctions
!                Hk  --- Modulus k ( 0 ó k ó 1 )
!       Output : ESN --- sn u
!                ECN --- cn u
!                EDN --- dn u
!                EPH --- phi ( in degrees )
!       ========================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION R(40)
        PI=3.14159265358979D0
        A0=1.0D0
        B0=DSQRT(1.0D0-HK*HK)
        DO 10 N=1,40
           A=(A0+B0)/2.0D0
           B=DSQRT(A0*B0)
           C=(A0-B0)/2.0D0
           R(N)=C/A
           IF (C.LT.1.0D-7) GO TO 15
           A0=A
10         B0=B
15      DN=2.0D0**N*A*U
        DO 20 J=N,1,-1
           T=R(J)*DSIN(DN)
           SA=DATAN(T/DSQRT(DABS(1.0D0-T*T)))
           D=.5D0*(DN+SA)
20         DN=D
        EPH=D*180.0D0/PI
        ESN=DSIN(D)
        ECN=DCOS(D)
        EDN=DSQRT(1.0D0-HK*HK*ESN*ESN)
        RETURN
        END





!Dreiecksfunktion
!eps2 hat (bisher) keine weitere Funktion
subroutine BASIS_FUNKTION(eps1, eps2, value)
    implicit none
    real *8, intent(in)                ::              eps1, eps2
    real *8, intent(out)               ::              value

    if(eps1.LE.0.5) then
        value = 2*eps1
    else
        value = 1-((0.5-eps1)*2)
    end if
    return
 end
