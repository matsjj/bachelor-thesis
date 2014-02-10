!include 'math_functions.f90'
!include 'calculate_kernel.f90'


program stab
!use math_functions
use calculate_kernel
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
integer                                         ::                segmente_pro_stab, punkte_pro_stab

!Koordinaten der Segmentmittelpunkte
!Ab jetzt werden die Segmente über IDs angesprochen, die über das Programm hinweg eindeutig sind
!coord (Stabnummer; Segmentnummer; Ort(1=Anfang, 2=Mitte, 3=Ende); Koordinate (1=x; 2=y; 3=z))
real *8, allocatable, dimension(:,:,:,:)          ::                coord

!Koordinaten der Punkte (Am Anfang und am Ende von Segmenten)
!punkte (Stabnummer; Punktnummer; Koordinate (s.o.) )
real *8, allocatable, dimension(:,:,:)          ::                punkte


!Vektor mit allen Strömen
complex *8, allocatable, dimension(:)      ::        stroeme
complex *8, allocatable, dimension(:)      ::        anregung

complex *8, allocatable, dimension(:)		::		stromverteilung

!Matrix
complex *8, allocatable, dimension(:,:)            ::                matrix, a_matrix
complex *8, allocatable, dimension(:,:)		::		  psi_matrix, e_matrix



!Temporäre Variablen
real *8                                         ::                rmax, roh, beta, K_von_beta, temp, summe, delta_L, &
                                                                        basis, eps1, eps2, weight, z, z_u, z_l, z_k, div, &
																		eps_0_1, eps_0_2, u_1_u, u_1_l, u_1_k, delta_z, delta_eps, &
																		delta_eps_1, delta_eps_2, R_1_k, eps_1_k, eps_2_k, abstand
complex *8					::		  kernel, csumme_1, csumme_2, ctemp, int_1_1, int_1_2, int_2_1, int_2_2, int_3_1, int_3_2, int_4_1, &
													int_4_2
integer                                         ::                i,j,k,l, id_obs, id_src, m, n, info, &
									int_start, int_end, modus
integer, allocatable, dimension(:)		::		  pivot
real *8                                         ::                b,c

real *8, dimension(3)                            ::                coord_z, coord_z_strich

!Integral Weights und Sample Points für Gauss-Legendre und MRW(Ma, Rouklin, Wandzura)
real *8, allocatable, dimension(:)                               ::                gauss_legendre_points, &
                                                                    gauss_legendre_weights, mrw_points, mrw_weights


include 'constants.f90'



staebe_anzahl = 1
punkte_pro_stab = 1000
segmente_pro_stab= punkte_pro_stab - 1
!Ordnung der Matrizen, Anzahl der Unbekanten, etc..
n = punkte_pro_stab * staebe_anzahl

fq = 3000000
wellenzahl = 2*PI*fq*sqrt(MUE_NULL * EPSILON_NULL)

allocate(staebe(staebe_anzahl,2,3))
allocate(coord(staebe_anzahl, segmente_pro_stab, 3, 3))
allocate(punkte(staebe_anzahl, punkte_pro_stab, 3))
allocate(matrix( n, n) )
allocate(a_matrix( n, n) )
allocate(psi_matrix( n, n) )
allocate(e_matrix( n, n) )
allocate(a(staebe_anzahl))
allocate(stroeme(n))
allocate(anregung(n))
allocate(pivot(n))
allocate(stromverteilung(n*2))

do i=1, n
	anregung(i) = (0.0, 0.0)
	stroeme(i) = (0.0, 0.0)
end do
anregung(10) = (1.0, 0.0)

!Integration wird mit 10 Schritten momentan fest geschrieben
!TODO: Routine schreiben / finden, die die Nodes und Weights generisch berechnet
include 'integrations_konstanten.f90'


!Radius der Staebe wird manuell eingefügt
do i=1, staebe_anzahl
        a(i) = 0.01
end do

!Koordinaten der Staebe manuell zuweisen
staebe(1, 1, 1) = 0
staebe(1, 1, 2) = 0
staebe(1, 1, 3) = 0

staebe(1, 2, 1) = 0
staebe(1, 2, 2) = 0
staebe(1, 2, 3) = 1

!staebe(2, 1, 1) = 10
!staebe(2, 1, 2) = 10
!staebe(2, 1, 3) = 0

!staebe(2, 2, 1) = 10
!staebe(2, 2, 2) = 10
!staebe(2, 2, 3) = 100

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
				b = (staebe(i, 2, k) - staebe(i, 1, k)) / (punkte_pro_stab+1)
                do j=1, punkte_pro_stab
                        punkte(i, j, k) = (j-0) * b + staebe(i, 1, k)
                end do        
	end do
	
end do
!write(10,*) coord
!Beobachtungspunkt (z) Hier gucken wir uns das Feld an
do i=1, staebe_anzahl
        do j=1, punkte_pro_stab
                id_obs = (i-1)*punkte_pro_stab + j

				!Punkt, für den das Feld bestimmt werden soll
                coord_z(1) = punkte(i, j, 1)
                coord_z(2) = punkte(i, j, 2)
                coord_z(3) = punkte(i, j, 3)

                !Quellpunkt (z') Das iterieren wir über den Beobachtungspunkt
                do l=1, staebe_anzahl
                        do k=1, punkte_pro_stab
                                id_src = (l-1)*punkte_pro_stab + k

				!Prüfen, ob Punkt am Anfang / Ende sitzt
				!Falls, wird nur über die Halbe länge iteriert
				!Basisfunktion muss entsprechend angepasst werden

                                !Wir gehen davon aus, dass parallel zur z Achse
                                !Daher ist delta_L nur von z abhängig
								if (k.EQ.1) then
									delta_L = 2*abs(punkte(l, k, 3) - punkte(l, k+1, 3))
								else if (k.EQ.punkte_pro_stab) then
									delta_L = 2*abs(punkte(l, k-1, 3) - punkte(l, k, 3))
								else
									delta_L = abs(punkte(l, k-1, 3) - punkte(l, k+1, 3))
								end if
								modus = 2

                                csumme_1 = 0
                                csumme_2 = 0

								abstand = sqrt( (punkte(i, j, 1)-punkte(l, k, 1))**2 + &
									 (punkte(i, j, 2)-punkte(l, k, 2))**2 + (punkte(i, j, 3)-punkte(l, k, 3))**2 )


!print *, 'Abstand: ', abstand, 'Bedingung größer als ', 5*a(l)


								!Einfacher Fall, Beobachtungspunkt weit genug von der Quelle entfernt
								if(abstand.GE.(5*a(l))) then
		                            do m=1, 10

		                                eps1 = mrw_points(m)
		                                eps2 = 1 - eps1
		                                weight = mrw_weights(m)


										!Punkt, von dem aus ein Feld ausgestrahlt wird
		                                coord_z_strich(1) = punkte(l, k, 1)
		                                coord_z_strich(2) = punkte(l, k, 2)
										coord_z_strich(3) = punkte(l, k, 3) - delta_L / 2 + delta_L*eps2

		                                call CALC_KERNEL(a(l), wellenzahl, coord_z, coord_z_strich, mrw_points, &
		                                                                mrw_weights, 10, kernel)
		                                call BASIS_FUNKTION(eps1, eps2, modus, basis)

		                                csumme_1 = csumme_1 + (weight * basis * kernel)
										if(eps1.LE.0.5)	then
											div = 2
										else
											div = -2
										end if
		                                csumme_2 = csumme_2 + (weight * basis * kernel * div) 

		                            end do
		                            csumme_1 = delta_L * csumme_1
		                            csumme_2 = delta_L * csumme_2
								else

									!Schwieriger Fall mit Quelle nahe am Beobachtungspunkt
									!Intervalle int_1 bis int_4 aufstellen (24)
									z = punkte(i, j, 3)
									z_u = punkte(l, k, 3) + delta_L / 2
									z_l = punkte(l, k, 3) - delta_L / 2
print *, 'Quelle ist bei ', z, ' Segment geht von ', z_l, 'bis', z_u
									eps_0_1 = (z_u - z) / delta_L
									eps_0_2 = (z - z_l) / delta_L
									delta_eps = 5*a(l) / delta_L
									roh = sqrt( (punkte(i, j, 1)-punkte(l, k, 1))**2 + (punkte(i, j, 2)-punkte(l, k, 2))**2 )


print *, 'Abbildung auf den Stab: eps_0_1: ', eps_0_1, 'eps_0_2', eps_0_2
if (eps_0_1.LT.0) then
	print *, 'Intervall 1 über das ganze Segment!'
	eps_0_1 = 0
	eps_0_2 = 1
end if
	
									!Intervall 1:
									delta_eps_1 = min( eps_0_2, delta_eps)
									u_1_l = (-1) * asinh( (delta_L * eps_0_2) / (roh + a(l) ) ) 
									u_1_u = (-1) * asinh( (delta_L * max(delta_eps_1, (-1) * eps_0_1) ) / (roh - a(l) ) )
!Vermeintlicher Fehler: Integral muss bis zum Ende des Segments berechnet werden, nicht bis zu Eps_0_1									
!u_1_u = (-1) * asinh( (delta_L * max(delta_eps_1, (-1) * eps_0_1) ) / (roh - a(l) ) )
									!print *, 'eps_0_1', eps_0_1, 'eps_0_2', eps_0_2, 'delta_eps', delta_eps, 'delta_eps_1', delta_eps_1
									!print *, 'z_l', z_l, 'z_u', z_u,  'u_1_l', u_1_l, 'u_1_u', u_1_u

                           		    csumme_1 = 0
                                	csumme_2 = 0
									do m=1, 10
										eps1 = mrw_points(m)
		                                eps2 = 1 - eps1
		                                weight = mrw_weights(m)

										u_1_k = u_1_l * eps1 + u_1_u * eps2
										R_1_k = (roh + a(l)) * cosh(u_1_k)
										z_k = z + (roh + a(l)) * sinh(eps1)
										eps_2_k = eps_0_2 + ( (roh - a(l)) / (delta_L) ) * sinh(u_1_k)
										eps_1_k = 1-eps_2_k

										
										coord_z_strich(1) = punkte(l, k, 1)
		                                coord_z_strich(2) = punkte(l, k, 2)
										coord_z_strich(3) = z_k	

		                                call CALC_KERNEL(a(l), wellenzahl, coord_z, coord_z_strich, mrw_points, &
		                                                                mrw_weights, 10, kernel)
		                                call BASIS_FUNKTION(eps_1_k, eps_2_k, modus, basis)

										if(eps_1_k.LE.0.5)	then
											div = 2
										else
											div = -2
										end if

										temp = ( ( weight * (u_1_u - u_1_l) * R_1_k ) / delta_L )

		                                csumme_1 = csumme_1 + (basis * kernel * temp) 
		                                csumme_2 = csumme_2 + (basis * kernel * temp * div) 
									end do
									int_1_1 = csumme_1 * delta_L
									int_1_2 = csumme_2 * delta_L

									eps_0_1 = (z_u - z) / delta_L
									eps_0_2 = (z - z_l) / delta_L

									if (eps_0_1.LT.0) then
										int_2_1 = (0.0, 0.0)
										int_2_2 = (0.0, 0.0)
										int_3_1 = (0.0, 0.0)
										int_3_2 = (0.0, 0.0)
										int_4_1 = (0.0, 0.0)
										int_4_2 = (0.0, 0.0)
									else



										!Intervall 4:
										!Im Folgenden werden die Variablen wiederverwendet
										!Es werden KEINE neuen Variablen für jedes Intervall verwendet
										delta_eps_2 = min( eps_0_1, delta_eps)
										u_1_u = asinh( (delta_L * eps_0_1) / (roh + a(l) ) ) 
										u_1_l = asinh( (delta_L * max(delta_eps_2, (-1) * eps_0_2) ) / (roh - a(l) ) )


		                       		    csumme_1 = 0
		                            	csumme_2 = 0
										do m=1, 10
											eps1 = mrw_points(m)
				                            eps2 = 1 - eps1
				                            weight = mrw_weights(m)

											u_1_k = u_1_l * eps1 + u_1_u * eps2
											R_1_k = (roh + a(l)) * cosh(u_1_k)
											z_k = z + (roh + a(l)) * sinh(eps1)
											eps_2_k = eps_0_2 + ( (roh - a(l)) / (delta_L) ) * sinh(u_1_k)
											eps_1_k = 1-eps_2_k

										
											coord_z_strich(1) = punkte(l, k, 1)
				                            coord_z_strich(2) = punkte(l, k, 2)
											coord_z_strich(3) = z_k	

				                            call CALC_KERNEL(a(l), wellenzahl, coord_z, coord_z_strich, mrw_points, &
				                                                            mrw_weights, 10, kernel)
				                            call BASIS_FUNKTION(eps_1_k, eps_2_k, modus, basis)

											if(eps_1_k.LE.0.5)	then
												div = 2
											else
												div = -2
											end if

											temp = ( ( weight * (u_1_u - u_1_l) * R_1_k ) / delta_L )

				                            csumme_1 = csumme_1 + (basis * kernel * temp) 
				                            csumme_2 = csumme_2 + (basis * kernel * temp * div) 
										end do
										int_4_1 = csumme_1 * delta_L
										int_4_2 = csumme_2 * delta_L

										if (eps_0_2.LT.0) then
											int_1_1 = (0.0, 0.0)
											int_1_2 = (0.0, 0.0)
											int_2_1 = (0.0, 0.0)
											int_2_2 = (0.0, 0.0)
											int_3_1 = (0.0, 0.0)
											int_3_2 = (0.0, 0.0)
										else



											!Intervall 2
											delta_eps_1 = min( eps_0_2, delta_eps)
											delta_eps_2 = min( eps_0_1, delta_eps)

				                   		    csumme_1 = 0
				                        	csumme_2 = 0
											do m=1, 10
												eps1 = mrw_points(m)
						                        eps2 = 1 - eps1
						                        weight = mrw_weights(m)


												eps_2_k = eps_0_2 - (eps2 * delta_eps_1)
												eps_1_k = 1-eps_2_k


												z_k = z - ( eps2 * delta_L * delta_eps_1 )
										
												coord_z_strich(1) = punkte(l, k, 1)
						                        coord_z_strich(2) = punkte(l, k, 2)
												coord_z_strich(3) = z_k	

						                        call CALC_KERNEL(a(l), wellenzahl, coord_z, coord_z_strich, mrw_points, &
						                                                        mrw_weights, 10, kernel)
						                        call BASIS_FUNKTION(eps_1_k, eps_2_k, modus, basis)

												if(eps_1_k.LE.0.5)	then
													div = 2
												else
													div = -2
												end if

												temp = ( ( weight * delta_eps_1 ) )

						                        csumme_1 = csumme_1 + (basis * kernel * temp) 
						                        csumme_2 = csumme_2 + (basis * kernel * temp * div) 
											end do
											int_2_1 = csumme_1 * delta_L
											int_2_2 = csumme_2 * delta_L



											!Intervall 3
											delta_eps_1 = min( eps_0_2, delta_eps)
											delta_eps_2 = min( eps_0_1, delta_eps)

				                   		    csumme_1 = 0
				                        	csumme_2 = 0
											do m=1, 10
												eps1 = mrw_points(m)
						                        eps2 = 1 - eps1
						                        weight = mrw_weights(m)


												eps_2_k = eps_0_2 + (eps2 * delta_eps_2)
												eps_1_k = 1-eps_2_k


												z_k = z + ( eps2 * delta_L * delta_eps_2 )
										
												coord_z_strich(1) = punkte(l, k, 1)
						                        coord_z_strich(2) = punkte(l, k, 2)
												coord_z_strich(3) = z_k	

						                        call CALC_KERNEL(a(l), wellenzahl, coord_z, coord_z_strich, mrw_points, &
						                                                        mrw_weights, 10, kernel)
						                        call BASIS_FUNKTION(eps_1_k, eps_2_k, modus, basis)

												if(eps_1_k.LE.0.5)	then
													div = 2
												else
													div = -2
												end if

												temp = ( ( weight * delta_eps_2 ) )

						                        csumme_1 = csumme_1 + (basis * kernel * temp) 
						                        csumme_2 = csumme_2 + (basis * kernel * temp * div) 
											end do
											int_3_1 = csumme_1 * delta_L
											int_3_2 = csumme_2 * delta_L
										end if
									end if

									


									csumme_1 = int_1_1 + int_2_1 + int_3_1 + int_4_1
									csumme_2 = int_1_2 + int_2_2 + int_3_2 + int_4_2
print *, '1: ', int_1_1, '2: ', int_2_1, '3: ', int_3_1, '4: ', int_4_1, 'SUMME: ', csumme_1



								end if
				!print *,  'OBS:       Stab: ',i,' Punkt: ',j,'      SRC:       Stab: ',l,' Punkt: ',k, ' csumme: ', csumme
                                !a_matrix(id_src, id_obs) = csumme_1 * MUE_NULL * (0.0, -1.0) * (2 * PI * fq) 
                                !psi_matrix(id_src, id_obs) = csumme_2 * ((-1) / ((0.0, 1.0) * EPSILON_NULL * ( 2 * PI * fq ) ))
								a_matrix(id_src, id_obs) = csumme_1
								psi_matrix(id_src, id_obs) = csumme_2
				!print *, csumme
                                !write (10, *) 'Summe: ', csumme



                         end do
                 end do
           end do
 end do


a_matrix = a_matrix * MUE_NULL
psi_matrix = (-1) / ((0.0, 1.0) * EPSILON_NULL * ( 2 * PI * fq ) ) * psi_matrix
e_matrix =  a_matrix - psi_matrix

!write (10, *) matrix
!write(10, *) anregung
!write(10, *) e_matrix

!write(10, *) a_matrix
!write(10, *) psi_matrix


!call cgbsv(n, 0, 0, 1, e_matrix, n, pivot, anregung, n, info)
call cgesv(n, 1, e_matrix, n, pivot, anregung, n, info)
!write(10, *) info

!write(10, *) anregung
!do i=1, 100
!	write(10, *) real(anregung(i))
!end do

!write(10, *) e_matrix

do i=1, n
!	write(10, *) aimag(anregung(i))
	write(10, *) real(anregung(i))
!	write(10, *) real(e_matrix(i, 3))
end do
end


