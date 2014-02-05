!Formel (20) in einer Subroutine
!

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
    if (roh.EQ.0) then 
	roh = wire_radius
    end if
    rmax = sqrt( ( (coord_z(3) - coord_z_strich(3))**2 ) + ( (wire_radius + roh ) ** 2 ) )

    beta = 2 * sqrt(roh * wire_radius) / rmax
    call COMELP(beta, K_von_beta)

!    print *, coord_z, '    ', coord_z_strich
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
	!print *, 'beta', beta, ' ctemp: ',c_temp, 'euler_factor: ',  wellenzahl * rmax * dn_result

    end do

    value = K_von_beta / ( 2 * (PI ** 2) * rmax) * c_summe
	!print *, ' c_summe: ',value, 'rmax: ', rmax, 'beta: ', beta

    return
    end



