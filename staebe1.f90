include 'math_functions.f90'
include 'calculate_kernel.f90'


program stab

!use mcomelp
implicit none

!Anzahl der Staebe
integer                                         ::                 staebe_anzahl

!Koordinaten von Anfang- und Endpunkten der Staebe
!staebe (Stabnummer; Art (1=Anfang; 2=Ende); Koordinate(1=x; 2=y; 3=z)
real *8, allocatable, dimension(:,:,:)          ::                staebe

!Radius der Staebe
real *8, allocatable, dimension(:)              ::                a

!Frequenz
real *8                                         ::                fq
real *8                                         ::                wellenzahl

!Segmente pro Stab
integer                                         ::                segmente_pro_stab

!Koordinaten der Segmentmittelpunkte
!Ab jetzt werden die Segmente über IDs angesprochen, die über das Programm hinweg eindeutig sind
!coord (Stabnummer; Segmentnummer; Ort(1=Anfang, 2=Mitte, 3=Ende); Koordinate (1=x; 2=y; 3=z))
real *8, allocatable, dimension(:,:,:,:)          ::                coord

!Vektor mit allen Strömen
real *8, allocatable, dimension(:)      ::        stroeme

!Matrix
real *8, allocatable, dimension(:,:)            ::                matrix



!Temporäre Variablen
real *8                                         ::                rmax, roh, beta, K_von_beta, temp, summe, delta_L, &
                                                                        basis, kernel, eps1, eps2, weight, z_u, z_l
integer                                         ::                i,j,k,l, id_obs, id_src, m, n
real *8                                         ::                b,c

real *8, dimension(3)                            ::                coord_z, coord_z_strich

!Integral Weights und Sample Points für Gauss-Legendre und MRW(Ma, Rouklin, Wandzura)
real *8, allocatable, dimension(:)                               ::                gauss_legendre_points, &
                                                                    gauss_legendre_weights, mrw_points, mrw_weights


include 'constants.f90'
include 'integrations_konstanten.f90'


staebe_anzahl = 2
segmente_pro_stab=5

fq = 100000
wellenzahl = 2*PI*fq*sqrt(MUE_NULL * EPSILON_NULL)

allocate(staebe(staebe_anzahl,2,3))
allocate(coord(staebe_anzahl, segmente_pro_stab, 3, 3))
allocate(matrix( segmente_pro_stab*staebe_anzahl, segmente_pro_stab*staebe_anzahl) )
allocate(a(staebe_anzahl))
allocate(stroeme(segmente_pro_stab * staebe_anzahl))

!Integration wird mit 10 Schritten momentan fest geschrieben
!TODO: Routine schreiben / finden, die die Nodes und Weights generisch berechnet



!Radius der Staebe wird manuell eingefügt
do i=1, staebe_anzahl
        a(i) = 5
end do

!Koordinaten der Staebe manuell zuweisen
staebe(1, 1, 1) = 0
staebe(1, 1, 2) = 0
staebe(1, 1, 3) = 0

staebe(1, 2, 1) = 0
staebe(1, 2, 2) = 0
staebe(1, 2, 3) = 100

staebe(2, 1, 1) = 10
staebe(2, 1, 2) = 10
staebe(2, 1, 3) = 0

staebe(2, 2, 1) = 10
staebe(2, 2, 2) = 10
staebe(2, 2, 3) = 100

open (10,file='staebe.out')
!Koordinaten der Segmentmittelpunkte ermitteln
do i=1, staebe_anzahl
        !Pro Koordinate einmal durchlaufen
        do k=1, 3
                !Entfernung zwischen Anfang und Ende
                b = (staebe(i, 2, k) - staebe(i, 1, k)) / segmente_pro_stab
                do j=1, segmente_pro_stab
                        coord(i, j, 1, k) = (j-1.0) * b + staebe(i, 1, k)
                        coord(i, j, 2, k) = (j-0.5) * b + staebe(i, 1, k)
                        coord(i, j, 3, k) = (j) * b + staebe(i, 1, k)
                end do
        end do
end do

!write(10,*) coord
!Beobachtungspunkt (z) Hier gucken wir uns das Feld an
do i=1, staebe_anzahl
        do j=1, segmente_pro_stab
                id_obs = (i-1)*segmente_pro_stab + j

                !Quellpunkt (z') Das iterieren wir über den Beobachtungspunkt
                do l=1, staebe_anzahl
                        do k=1, segmente_pro_stab
                                id_src = (l-1)*segmente_pro_stab + k
                                roh = sqrt( (coord(i,j,2,1)-coord(l,k,2,1))**2 + (coord(i,j,2,2)-coord(l,k,2,2))**2 )
                                !rmax = sqrt( ( (coord(i,j,3) - coord(l,k,3))**2 ) + ( (a(i) + roh ) ** 2 ) )
                                !write (10,*) rmax

                                !beta = 2 * sqrt(roh * a(i)) / rmax
                                !call COMELP(beta, K_von_beta)


                                !write (10, *) 'OBS:       Stab: ',i,' Segment: ',j,'      SRC:       Stab: ',l,' Segment: ',k
                                rmax = sqrt( ( (coord(i,j,2,3) - coord(l,k,2,3))**2 ) + ( (a(i) + roh ) ** 2 ) )

                                !write (10, *) 'SRC:       Stab: ',l,' Segment: ',k
                                !write (10, *) 'roh: ', roh, 'K: ', K_von_beta

            !SUBROUTINE CALC_KERNEL(wire_radius, wellenzahl, coord_z, coord_z_strich, integrating_weights, integrating_nodes, integrating_limit, return_value)

                                coord_z(1) = coord(i, j, 2, 1)
                                coord_z(2) = coord(i, j, 2, 2)
                                coord_z(3) = coord(i, j, 2, 3)

                                summe = 0
                                do m=1, 10

                                    eps1 = gauss_legendre_points(m)
                                    eps2 = 1 - eps1
                                    weight = gauss_legendre_weights(m)

                                    !Wir gehen davon aus, dass parallel zur z Achse
                                    !Daher ist delta_L nur von z abhängig
                                    delta_L = abs(coord(l, k, 3, 3) - coord(l, k, 1, 3))



                                    coord_z_strich(1) = coord(l, k, 1, 1)
                                    coord_z_strich(2) = coord(l, k, 1, 2)
                                    coord_z_strich(3) = coord(l, k, 1, 3) + delta_L*eps2

                                    call CALC_KERNEL(a(l), wellenzahl, coord_z, coord_z_strich, gauss_legendre_points, &
                                                                    gauss_legendre_weights, 10, kernel)
                                    call BASIS_FUNKTION(eps1, eps2, basis)

                                    summe = summe + (weight * basis * kernel)


                                end do
                                summe = delta_L * summe
                                matrix(id_obs, id_src) = summe
                                !write (10, *) 'Summe: ', summe



                         end do
                 end do
           end do
 end do

end



