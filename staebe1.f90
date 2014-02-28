!include 'math_functions.f90'
!include 'calculate_kernel.f90'


program stab
!use math_functions
use calculate_kernel
!use mcomelp
implicit none

!Anzahl der Staebe
integer                                         	::  	staebe_anzahl

!Koordinaten von Anfang- und Endpunkten der Staebe
!staebe (Stabnummer; Art (1=Anfang; 2=Ende); Koordinate(1=x; 2=y; 3=z)
real *8, allocatable, dimension(:,:,:)          	::   	staebe

!Radius der Staebe
real *8, allocatable, dimension(:)              	::   	a

!Frequenz
real *8                                         	::  	fq
real *8                                         	:: 		wellenzahl

!Segmente pro Stab
integer                                         	::  	segmente_pro_stab, punkte_pro_stab

!Koordinaten der Segmentmittelpunkte
!Ab jetzt werden die Segmente über IDs angesprochen, die über das Programm hinweg eindeutig sind
!coord (Stabnummer; Segmentnummer; Ort(1=Anfang, 2=Mitte, 3=Ende); Koordinate (1=x; 2=y; 3=z))
real *8, allocatable, dimension(:,:,:,:)          	::    	coord

!Koordinaten der Punkte (Am Anfang und am Ende von Segmenten)
!punkte (Stabnummer; Punktnummer; Koordinate (s.o.) )
!real *8, allocatable, dimension(:,:,:)          	::  	punkte

!Aufpunkte
!Die Aufpunkte liegen immer in der Mitte einer steigenden bzw. fallenden Flanke, bzw in der Mitte des Segments
!aufpunkte (Stabnummer; Punktnummer; Art (1=Steigend; 2=Mitte, 3=Fallend); Koordinate)
real *8, allocatable, dimension(:,:,:,:)			::		punkte


!Vektor mit allen Strömen
complex *8, allocatable, dimension(:)      			::  	stroeme
complex *8, allocatable, dimension(:)      			:: 		anregung

complex *8, allocatable, dimension(:)				::		stromverteilung

!Matrix
complex *8, allocatable, dimension(:,:)      		:: 		matrix, a_matrix
complex *8, allocatable, dimension(:,:)				::		psi_matrix, e_matrix

!Matrix mit Beziehung von Aufpunkten und Quellpunkten
!aufpunkt_matrix (Quellpunkt, Zielpunkt, Art (1=Steigend, 2=Fallend)
complex *8, allocatable, dimension(:,:,:)			::		aufpunkt_matrix



!Temporäre Variablen
real *8                                         	::		rmax, roh, beta, K_von_beta, temp, summe, delta_L, &
                                                                        basis, eps1, eps2, weight, z, z_u, z_l, z_k, div, &
																		eps_0_1, eps_0_2, u_1_u, u_1_l, u_1_k, delta_z, delta_eps, &
																		delta_eps_1, delta_eps_2, R_1_k, eps_1_k, eps_2_k, abstand, abst

complex *8											::		kernel, csumme_1, csumme_2, ctemp, int_1_1, int_1_2, int_2_1, int_2_2,&
																		int_3_1, int_3_2, int_4_1, int_4_2, c2

integer                                         	::  	i,j,k,l, id_obs, id_src, m, n, info, punktart, &
																		int_start, int_end, modus, mrw_n

integer, allocatable, dimension(:)					::		pivot
real *8                                         	:: 		b,c

real *8, dimension(3)                            	:: 		coord_z, coord_z_strich

!Integral Weights und Sample Points für Gauss-Legendre und MRW(Ma, Rouklin, Wandzura)
real *8, allocatable, dimension(:)   				:: 		gl_points, &
                                                                    gl_weights, mrw_points, mrw_weights


include 'constants.f90'

!Integration wird mit 20 Schritten momentan fest geschrieben
!TODO: Routine schreiben / finden, die die Nodes und Weights generisch berechnet
include 'integrations_konstanten.f90'

staebe_anzahl = 1
punkte_pro_stab = 151
segmente_pro_stab= punkte_pro_stab - 1

!Ordnung der Matrizen, Anzahl der Unbekanten, etc..
n = punkte_pro_stab * staebe_anzahl

fq = 2000000000
wellenzahl = 2*PI*fq*sqrt(MUE_NULL * EPSILON_NULL)

allocate(staebe(staebe_anzahl,2,3))
allocate(coord(staebe_anzahl, segmente_pro_stab, 3, 3))
!allocate(punkte(staebe_anzahl, punkte_pro_stab, 3))
allocate(punkte(staebe_anzahl, punkte_pro_stab, 3, 3))
allocate(matrix( n, n) )
allocate(a_matrix( n, n) )
allocate(psi_matrix( n, n) )
allocate(e_matrix( n, n) )
allocate(aufpunkt_matrix( n, n, 2))
allocate(a(staebe_anzahl))
allocate(stroeme(n))
allocate(anregung(n))
allocate(pivot(n))
allocate(stromverteilung(n*2))

do i=1, n
	anregung(i) = (0.0, 0.0)
	stroeme(i) = (0.0, 0.0)
end do
anregung(punkte_pro_stab/2+1) = (-1.0, 0.0)



print *, 'Wellenzahl: ', wellenzahl
print *, 'Wellenlänge: ', fq*sqrt(MUE_NULL * EPSILON_NULL)


!Radius der Staebe wird manuell eingefügt
do i=1, staebe_anzahl
        a(i) = 0.00001
end do

print *, 'ka << 1: ka = ',wellenzahl * a(1)

!Koordinaten der Staebe manuell zuweisen
staebe(1, 1, 1) = 0
staebe(1, 1, 2) = 0
staebe(1, 1, 3) = 0

staebe(1, 2, 1) = 0
staebe(1, 2, 2) = 0
staebe(1, 2, 3) = 1

!staebe(2, 1, 1) = 0
!staebe(2, 1, 2) = 0
!staebe(2, 1, 3) = 1

!staebe(2, 2, 1) = 0
!staebe(2, 2, 2) = 0
!staebe(2, 2, 3) = 2

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
                !do j=1, punkte_pro_stab
                !        punkte(i, j, k) = (j) * b + staebe(i, 1, k)
                !end do        
				!Aufpunkte zwischen den Punkten bestimmen
				do j=1, punkte_pro_stab
						punkte(i, j, 1, k) = (j-0.5) * b + staebe(i, 1, k)
						punkte(i, j, 2, k) = (j) * b + staebe(i, 1, k)
						punkte(i, j, 3, k) = (j+0.5) * b + staebe(i, 1, k)
						!print *, punkte(i, j, 1, k), punkte(i, j, 2, k), punkte(i, j, 3, k)
				end do
	end do
	
end do
!write(10,*) coord
!Beobachtungspunkt (z) Hier gucken wir uns das Feld an
do i=1, staebe_anzahl
        do j=1, punkte_pro_stab
			!Aufpunkt Anfang, Aufpunkt Mitte, Aufpunkt Ende?
			do punktart = 1, 3
                id_obs = (i-1)*punkte_pro_stab + j

				!Punkt, für den das Feld bestimmt werden soll
                coord_z(1) = punkte(i, j, punktart, 1)
                coord_z(2) = punkte(i, j, punktart, 2)
                coord_z(3) = punkte(i, j, punktart, 3)


                !Quellpunkt (z') Das iterieren wir über den Beobachtungspunkt
                do l=1, staebe_anzahl
                        do k=1, punkte_pro_stab
                                id_src = (l-1)*punkte_pro_stab + k

                                !Wir gehen davon aus, dass parallel zur z Achse
                                !Daher ist delta_L nur von z abhängig
								!Für Quellpunkte wird immer der Punkt in der Mitte des Segments gewählt = 2
								if (k.EQ.1) then
									delta_L = 2*abs(punkte(l, k, 2, 3) - punkte(l, k+1, 2, 3))
								else if (k.EQ.punkte_pro_stab) then
									delta_L = 2*abs(punkte(l, k-1, 2, 3) - punkte(l, k, 2, 3))
								else
									delta_L = abs(punkte(l, k-1, 2, 3) - punkte(l, k+1, 2, 3))
								end if

								modus = 2

                                csumme_1 = (0.0, 0.0)
                                csumme_2 = (0.0, 0.0)

								abstand = sqrt( (punkte(i, j, punktart, 1)-punkte(l, k, 2, 1))**2 + &
									 (punkte(i, j, punktart, 2)-punkte(l, k, 2, 2))**2 + (punkte(i, j, punktart, 3)-punkte(l, k, 2, 3))**2 )


!print *, 'Abstand: ', abstand, 'Bedingung größer als ', 10*a(l)
!print *, coord_z, delta_L
								!Einfacher Fall, Beobachtungspunkt weit genug von der Quelle entfernt
								if(abstand.GE.(0*a(l))) then
		                            do m=1, mrw_n

		                                eps1 = gl_points(m)
		                                eps2 = 1 - eps1
		                                weight = gl_weights(m)

										
								
										!Punkt, von dem aus ein Feld ausgestrahlt wird
		                                coord_z_strich(1) = punkte(l, k, 2, 1)
		                                coord_z_strich(2) = punkte(l, k, 2, 2)
										coord_z_strich(3) = punkte(l, k, 2, 3) - delta_L / 2 + delta_L*eps2
!print *, coord_z_strich

		                                call CALC_KERNEL(a(l), wellenzahl, coord_z, coord_z_strich, gl_points, &
		                                                                gl_weights, mrw_n, kernel)
		                                call BASIS_FUNKTION(eps1, eps2, modus, basis)

										call BASIS_CHARGE(eps1, eps2, div)

		                                csumme_1 = csumme_1 + (weight * kernel * basis)
										csumme_2 = csumme_2 + (weight * kernel * div / delta_L) 
!print *,  coord_z_strich, (weight * kernel * basis), (weight * kernel * div) 
!csumme_1 = csumme_1 + (weight * eps1)
										!csumme_2 = csumme_1
										!print *, 'csumme_1', csumme_1, 'csumme_2', csumme_2
!if(id_obs.eq.5 .and. punktart.eq.1) write(10, *) real(csumme_2)
											if ((id_obs.eq.24 .or. id_obs.eq.25) .and. id_src.eq.31 .and. punktart.eq.1) then
!write(10, *) real(kernel)
										!		print *, coord_z(3), coord_z_strich(3), (coord_z(3) - coord_z_strich(3)), div, kernel, csumme_1, csumme_2
											end if

		                            end do	


					!print *, ''
					!print *, csumme_1
					!print*, ' Summe2', csumme_2
					!print *, 
		                            csumme_1 = delta_L * csumme_1
		                            csumme_2 = delta_L * csumme_2
									c2 = csumme_2
									!if(id_src.eq.21 .and. punktart.eq.2) print *, punktart, abstand, csumme_1
if (id_src.eq.31 .and. punktart.eq.2) then
!	print *, csumme_1, csumme_2
!	write(10, *) real(csumme_2)
end if

								else

									!Schwieriger Fall mit Quelle nahe am Beobachtungspunkt
									!Intervalle int_1 bis int_4 aufstellen (24)
									int_1_1 = (0.0, 0.0)
									int_1_2 = (0.0, 0.0)
									int_2_1 = (0.0, 0.0)
									int_2_2 = (0.0, 0.0)
									int_3_1 = (0.0, 0.0)
									int_3_2 = (0.0, 0.0)
									int_4_1 = (0.0, 0.0)
									int_4_2 = (0.0, 0.0)
									z = punkte(i, j, punktart, 3)
									z_u = punkte(l, k, 2, 3) + delta_L / 2
									z_l = punkte(l, k, 2, 3) - delta_L / 2

									eps_0_1 = (z_u - z) / delta_L
									eps_0_2 = (z - z_l) / delta_L
									delta_eps = 5*a(l) / delta_L
									roh = sqrt( (punkte(i, j, punktart, 1)-punkte(l, k, 2, 1))**2 + (punkte(i, j, punktart, 2)-punkte(l, k, 2, 2))**2 )
									!if (punktart.ne.2) print *, punktart, z
if (id_src.eq.31 .and. punktart.eq.1) then
!print *,'eps01', eps_0_1, 'eps02', eps_0_2, 'delta_eps', delta_eps
!print *, 'Obs Koordinate: ',z, ' Segment Obs: ', id_obs
end if

								if( eps_0_2 .GT. delta_eps ) then
									!Intervall 1:
									delta_eps_1 = min( eps_0_2, delta_eps)
									u_1_l = (-1) * asinh( (delta_L * eps_0_2) / (roh + a(l) ) ) 
									u_1_u = (-1) * asinh( (delta_L * max(delta_eps_1, (-1) * eps_0_1))  / (roh + a(l) ) )

                           		    csumme_1 = 0
                                	csumme_2 = 0
									do m=1, mrw_n
										eps1 = gl_points(m)
		                                eps2 = 1 - eps1
		                                weight = gl_weights(m)/20

										u_1_k = u_1_l * eps2 + u_1_u * eps1
										R_1_k = (roh + a(l)) * cosh(u_1_k)
										z_k = z + (roh + a(l)) * sinh(u_1_k)
										eps_2_k = eps_0_2 + ( (roh + a(l)) / (delta_L) ) * sinh(u_1_k)
										eps_1_k = 1-eps_2_k

										
										coord_z_strich(1) = punkte(l, k, 2, 1)
		                                coord_z_strich(2) = punkte(l, k, 2, 2)
										coord_z_strich(3) = z_k	

		                                call CALC_KERNEL(a(l), wellenzahl, coord_z, coord_z_strich, gl_points, &
		                                                                gl_weights, mrw_n, kernel)
		                                call BASIS_FUNKTION(eps_1_k, eps_2_k, modus, basis)

										call BASIS_CHARGE(eps_2_k, eps2, div)

										temp = ( ( weight * (u_1_u - u_1_l) * R_1_k ) / delta_L )

		                                csumme_1 = csumme_1 + (basis * kernel * temp) 
		                                csumme_2 = csumme_2 + (kernel * temp * div) /10
									end do
									int_1_1 = csumme_1 * delta_L
									int_1_2 = csumme_2 * delta_L

									eps_0_1 = (z_u - z) / delta_L
									eps_0_2 = (z - z_l) / delta_L

								end if

										!Intervall 4:
										!Im Folgenden werden die Variablen wiederverwendet
										!Es werden KEINE neuen Variablen für jedes Intervall verwendet

								if (eps_0_1 .GT. (delta_eps) ) then
										delta_eps_2 = min( eps_0_1, delta_eps)
										u_1_u = asinh( (delta_L * eps_0_1) / (roh + a(l) ) ) 
										u_1_l = asinh( (delta_L * max(delta_eps_2, (-1) * eps_0_2) ) / (roh + a(l) ) )


		                       		    csumme_1 = 0
		                            	csumme_2 = 0
										do m=1, mrw_n
											eps1 = gl_points(m)
				                            eps2 = 1 - eps1
				                            weight = gl_weights(m)

											u_1_k = u_1_l * eps2 + u_1_u * eps1
											R_1_k = (roh + a(l)) * cosh(u_1_k)
											z_k = z + (roh + a(l)) * sinh(u_1_k)
											eps_2_k = eps_0_2 + ( (roh + a(l)) / (delta_L) ) * sinh(u_1_k)
											eps_1_k = 1-eps_2_k

										
											coord_z_strich(1) = punkte(l, k, 2, 1)
				                            coord_z_strich(2) = punkte(l, k, 2, 2)
											coord_z_strich(3) = z_k	

				                            call CALC_KERNEL(a(l), wellenzahl, coord_z, coord_z_strich, gl_points, &
				                                                            gl_weights, mrw_n, kernel)
				                            call BASIS_FUNKTION(eps_1_k, eps_2_k, modus, basis)

											call BASIS_CHARGE(eps_2_k, eps2, div)

											temp = ( ( weight * (u_1_u - u_1_l) * R_1_k ) / delta_L )

				                            csumme_1 = csumme_1 + (basis * kernel * temp) 
				                            csumme_2 = csumme_2 + (kernel * temp *  div)/10 

											if ((id_obs.eq.24 .or. id_obs.eq.25) .and. id_src.eq.31 .and. punktart.eq.1) then
!write(10, *) real(kernel)
												print *, coord_z(3), coord_z_strich(3), (coord_z(3) - coord_z_strich(3)) , div, kernel, temp, csumme_1, csumme_2
											end if
										end do
										int_4_1 = csumme_1 * delta_L
										int_4_2 = csumme_2 * delta_L
											if ( (id_obs.eq.24 .or. id_obs.eq.25) .and. id_src.eq.31 .and. punktart.eq.1) then
												print *, '1: ', int_4_1, '2: ', int_4_2
											end if

								end if

								if (eps_0_1 .GT. (-delta_eps) .AND. eps_0_2 .GT. (-delta_eps) ) then 


											!Intervall 2
											delta_eps_1 = min( eps_0_2, delta_eps)
											delta_eps_2 = min( eps_0_1, delta_eps)

				                   		    csumme_1 = 0
				                        	csumme_2 = 0
											do m=1, mrw_n
												eps1 = mrw_points(m)
						                        eps2 = 1 - eps1
						                        weight = mrw_weights(m)


												eps_2_k = eps_0_2 - (eps2 * delta_eps_1)
												eps_1_k = 1-eps_2_k


												z_k = z - ( eps2 * delta_L * delta_eps_1 )
										
												coord_z_strich(1) = punkte(l, k, 2, 1)
						                        coord_z_strich(2) = punkte(l, k, 2, 2)
												coord_z_strich(3) = z_k	

						                        call CALC_KERNEL(a(l), wellenzahl, coord_z, coord_z_strich, gl_points, &
						                                                        gl_weights, mrw_n, kernel)
												!print *, aimag(kernel)
						                        call BASIS_FUNKTION(eps_1_k, eps_2_k, modus, basis)

												call BASIS_CHARGE(eps_1_k, eps2, div)

												temp = ( ( weight * delta_eps_1 ) )

						                        csumme_1 = csumme_1 + (basis * kernel * temp) 
						                        csumme_2 = csumme_2 + (kernel * temp * div) 
											end do
											int_2_1 = csumme_1 * delta_L
											int_2_2 = csumme_2 * delta_L


											!Intervall 3
											delta_eps_1 = min( eps_0_2, delta_eps)
											delta_eps_2 = min( eps_0_1, delta_eps)

				                   		    csumme_1 = 0
				                        	csumme_2 = 0
											do m=1, mrw_n
												eps1 = mrw_points(m)
						                        eps2 = 1 - eps1
						                        weight = mrw_weights(m)


												eps_2_k = eps_0_2 + (eps2 * delta_eps_2)
												eps_1_k = 1-eps_2_k


												z_k = z + ( eps2 * delta_L * delta_eps_2 )
										
												coord_z_strich(1) = punkte(l, k, 2, 1)
						                        coord_z_strich(2) = punkte(l, k, 2, 2)
												coord_z_strich(3) = z_k	

						                        call CALC_KERNEL(a(l), wellenzahl, coord_z, coord_z_strich, gl_points, &
						                                                        gl_weights, mrw_n, kernel)
						                        call BASIS_FUNKTION(eps_1_k, eps_2_k, modus, basis)

												call BASIS_CHARGE(eps_1_k, eps2, div)

												temp = ( ( weight * delta_eps_2 ) )

						                        csumme_1 = csumme_1 + (basis * kernel * temp) 
						                        csumme_2 = csumme_2 + (kernel * temp * div) 
											end do
											int_3_1 = csumme_1 * delta_L
											int_3_2 = csumme_2 * delta_L

										
									end if

									


									csumme_1 = int_1_1 + int_2_1 + int_3_1 + int_4_1
									csumme_2 = int_1_2 + int_2_2 + int_3_2 + int_4_2
											if ( (id_obs.eq.23 .or. id_obs.eq.24) .and. id_src.eq.31 .and. punktart.eq.1) then
											!	print *, '1: ', int_4_1, '2: ', int_4_2, 'SUMME: ', csumme_2
											end if
!print *, '1: ', int_1_1, '2: ', int_2_1, '3: ', int_3_1, '4: ', int_4_1, 'SUMME: ', csumme_1
!print *, id_src
if(punktart.eq.1 .and. id_src.eq.31) then
!print *, punktart, id_src
print *, '1: ', int_1_1, '2: ', int_2_1, '3: ', int_3_1, '4: ', int_4_1, 'SUMME: ', csumme_1
print *, '1: ', int_1_2, '2: ', int_2_2, '3: ', int_3_2, '4: ', int_4_2, 'SUMME: ', csumme_2
!write (10, *) aimag(csumme_2)
!write (10, *) aimag(c2)
end if


								end if
									abst = sqrt( (punkte(i, j, punktart, 1)-punkte(l, k, 2, 1))**2 + &
											(punkte(i, j, punktart, 2)-punkte(l, k, 2, 2))**2 + (punkte(i, j, punktart, 3)-punkte(l, k, 2, 3))**2 ) 
									!abst = sqrt( (punkte(i, j, 1, 1)-punkte(i, j, 3, 1))**2 + &
									!		(punkte(i, j, 1, 2)-punkte(i, j, 3, 2))**2 + (punkte(i, j, 1, 3)-punkte(i, j, 3, 3))**2 ) 

!print *, punktart, id_src, csumme_1
								if (punktart.EQ.1) then
									aufpunkt_matrix(id_src, id_obs, 1) = csumme_2
								else if (punktart.EQ.3) then
									aufpunkt_matrix(id_src, id_obs, 2) = csumme_2
								else
!print *, abst
									a_matrix(id_src, id_obs) = csumme_1 * abst
								end if

								if(punktart.EQ.3) then

!print *, 'Koordinaten Aufpunkt 1: ' ,punkte(i, j, 1, 3), 'Koordinaten Aufpunkt 2: ' ,&
!										punkte(i, j, 3, 3), 'Koordinaten Quelle: ', punkte(l, k, 2, 3)
!print *, 'Ergebnis 1: ', aufpunkt_matrix(id_src, id_obs, 1), 'Ergebnis 2: ', aufpunkt_matrix(id_src, id_obs, 2) 


									abst = sqrt( (punkte(i, j, 1, 1)-punkte(i, j, 3, 1))**2 + &
											(punkte(i, j, 1, 2)-punkte(i, j, 3, 2))**2 + (punkte(i, j, 1, 3)-punkte(i, j, 3, 3))**2 ) 
									psi_matrix(id_src, id_obs) = (aufpunkt_matrix(id_src, id_obs, 2) &
											- aufpunkt_matrix(id_src, id_obs, 1) )!  / abst
!print *, 'Differenz: ', (aufpunkt_matrix(id_src, id_obs, 2) - aufpunkt_matrix(id_src, id_obs, 1)), &
!		 'Abstand: ', abst, 'gradient: ', psi_matrix(id_src, id_obs) 
									!print *,  'OBS:       Stab: ',i,' Punkt: ',j,'      SRC:       Stab: ',l,' Punkt: ',k, ' csumme: ', csumme_1
									!print *, aufpunkt_matrix(id_src, id_obs, 2), aufpunkt_matrix(id_src, id_obs, 1), 'psi:', psi_matrix(id_src, id_obs)
								end if

								!psi_matrix(id_src, id_obs) = csumme_2
				!print *, csumme
                                !write (10, *) 'Summe: ', csumme


                         end do
					end do
                 end do
			!print *, ''
           end do
 end do

!do i=1, n
!	write(10, *) aufpunkt_matrix(i, :, 1)
!end do

!do i=1, n
!	write(10, *) aufpunkt_matrix(i, :, 2)
!end do

!write(10, *) a_matrix
a_matrix = ( (0.0, -1.0) * (2 * PI * fq ) * a_matrix )  * MUE_NULL
psi_matrix = (-1) / ((0.0, 1.0) * EPSILON_NULL * ( 2 * PI * fq ) ) * psi_matrix
e_matrix =  a_matrix - psi_matrix
print *,  (0.0, -1.0) * (2 * PI * fq ) *  MUE_NULL
print *, (-1) / ((0.0, 1.0) * EPSILON_NULL * ( 2 * PI * fq ) )
!print *, a_matrix
!print *, ''
!print *, psi_matrix
!print *, ''
!print *, e_matrix

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
!print *, a_matrix(31, n)
	write(10, *) real(anregung(i))
!	write(10, *) real(e_matrix(punkte_pro_stab/2+1, i))
!	write(10, *) aimag(anregung(i))
!	write (10, *) sqrt(real(anregung(i))**2 + aimag(anregung(i))**2)
!	write(10, *) real(psi_matrix(76, i))
!	write(10, *) aimag(psi_matrix(punkte_pro_stab/2+1, i))
!write(10, *) aimag(e_matrix(51, i))
!	write(10, *) aimag(aufpunkt_matrix(31, i, 1))
!	write(10, *) real(aufpunkt_matrix(i, 8, 2))
!	write(10, *) aimag(aufpunkt_matrix(i, 8, 2))
!	write (10, *) real(e_matrix(11, i))
!	do j=1, n
!		print *, aufpunkt_matrix(i, j, 1)
!	end do
!	print *, ''
!write(10, *) aimag(aufpunkt_matrix(76, i, 1))
end do
do i=1, n
	write(10, *) aimag(anregung(i))
!	write(10, *) aimag(e_matrix(punkte_pro_stab/2+1, i))
!write(10, *) (aimag(aufpunkt_matrix(punkte_pro_stab/2, i, 2))-aimag(aufpunkt_matrix(punkte_pro_stab/2, i, 1)))/abst
!write(10, *) aimag(aufpunkt_matrix(50, i, 2))
!write(10, *) aimag(aufpunkt_matrix(50, i, 1))
!	write(10, *) aimag(anregung(i))
end do
end


