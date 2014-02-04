!Formel (20) in einer Subroutine
!

SUBROUTINE CALC_KERNEL(wire_radius, wellenzahl, coord_z, coord_z_strich, integrating_weights, integrating_nodes, &
                                integrating_limit, value)
    implicit none
    real *8, intent(in)                                 ::          wire_radius, wellenzahl
    integer, intent(in)                                 ::          integrating_limit
    real *8, intent(in), dimension(integrating_limit)   ::          integrating_weights, integrating_nodes
    real *8, intent(in), dimension(3)                   ::          coord_z, coord_z_strich
    real *8, intent(inout)                              ::          value


    real *8                                             ::          roh, K_von_beta, beta, rmax, summe, temp, dn_result, dn_argument
    integer                                             ::          j

    include 'constants.f90'

    summe = 0

    roh = sqrt( (coord_z(1)-coord_z_strich(1))**2 + (coord_z(2)-coord_z_strich(2))**2 )
    rmax = sqrt( ( (coord_z(3) - coord_z_strich(3))**2 ) + ( (wire_radius + roh ) ** 2 ) )

    beta = 2 * sqrt(roh * wire_radius) / rmax
    call COMELP(beta, K_von_beta)


    do j = 1, integrating_limit
        dn_argument = K_von_beta * integrating_nodes(j)
        call JELP(dn_argument, beta, temp, temp, dn_result, temp)

        temp = EULER ** ( (-1) * j * wellenzahl * rmax * dn_result )

        summe = summe + ( temp * integrating_weights(j))

    end do

    value = K_von_beta / ( 2 * (PI ** 2) * rmax) * summe

    return
    end




