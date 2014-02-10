!Formel (20) in einer Subroutine
!

module calculate_kernel
contains
SUBROUTINE CALC_KERNEL(wire_radius, wellenzahl, coord_z, coord_z_strich, integrating_nodes, integrating_weights, &
                                integrating_limit, value)
    implicit none
    real *8, intent(in)                                 ::          wire_radius, wellenzahl
    integer, intent(in)                                 ::          integrating_limit
    real *8, intent(in), dimension(integrating_limit)   ::          integrating_weights, integrating_nodes
    real *8, intent(in), dimension(3)                   ::          coord_z, coord_z_strich
    complex *8, intent(inout)                              ::          value


    real *8                                             ::          roh, K_von_beta, beta, rmax, summe, temp, dn_result, &
								dn_argument, am, sn, phi
    complex *8						::	    c_temp, c_summe
    integer                                             ::          j

    include 'constants.f90'

    roh = sqrt( (coord_z(1)-coord_z_strich(1))**2 + (coord_z(2)-coord_z_strich(2))**2 )
!    if (roh.EQ.0) then 
!	roh = wire_radius
!    end if
    rmax = sqrt( ( (coord_z(3) - coord_z_strich(3))**2 ) + ( (wire_radius + roh ) ** 2 ) )

    beta = 2 * sqrt(roh * wire_radius) / rmax
    call COMELP(beta, K_von_beta)

    !print *, coord_z, '    ', coord_z_strich
	c_summe = 0
    do j = 1, integrating_limit
        dn_argument = K_von_beta * integrating_nodes(j)
!dn_argument = 0.5
!beta = 0.5

        call JELP(dn_argument, beta, am, sn, dn_result, phi)
!print *, dn_argument, '  ', beta, dn_result
        c_temp = EULER ** ( (0.0, -1.0) * wellenzahl * rmax * dn_result )

	!print*, c_temp
        c_summe = c_summe + ( c_temp * integrating_weights(j))
!	print *, 'beta', beta, ' ctemp: ',c_temp, 'euler_factor: ',  wellenzahl * rmax * dn_result

    end do

    value = K_von_beta / ( 2 * (PI ** 2) * rmax) * c_summe
	!print *, ' c_summe: ',value, 'rmax: ', rmax, 'beta: ', beta ,'coord_z: ', &
	!	coord_z(3), 'coord_z_strich', coord_z_strich(3)

!	print *, ' K_von_beta: ',K_von_beta, 'rmax: ', rmax, 'beta: ', beta ,'coord_z: ', &
!		coord_z(3), 'coord_z_strich', coord_z_strich(3)

    
    end subroutine




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
        END subroutine


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
        
        END subroutine





!Dreiecksfunktion
!eps2 hat (bisher) keine weitere Funktion
!modus = 1 -> nur steigend
!modus = 2 -> normal
!modul = 3 -> nur fallend

subroutine BASIS_FUNKTION(eps1, eps2, modus, value)
    implicit none
    real *8, intent(in)                ::              eps1, eps2
    integer, intent(in)		       ::		modus
    real *8, intent(out)               ::              value
	if (modus.EQ.2) then
		if(eps1.LE.0.5) then
		    value = 2*eps1
		else
		    value = (1-(eps1))*2
		end if
		!print *, 'eps1: ', eps1, 'value: ', value
	else if (modus.EQ.1) then
		value = eps1
	else 
		value = 1-eps1
	end if
    end subroutine
end module



