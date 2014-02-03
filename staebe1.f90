module global
implicit none
integer                                         ::                test
end module global

program stab
use global
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
!coord (Stabnummer; Segmentnummer; Koordinate (1=x; 2=y; 3=z))
real *8, allocatable, dimension(:,:,:)          ::                coord

!Vektor mit allen Strömen
!Ab jetzt werden die Segmente über IDs angesprochen, die über das Programm hinweg eindeutig sind
real *8, allocatable, dimension(:)		::		  stroeme

!Matrix
real *8, allocatable, dimension(:,:)            ::                matrix



!Temporäre Variablen
real *8                                         ::                rmax
integer                                         ::                i,j,k
real *8                                         ::                b,c

!Integral Weights und Sample Points für Gauss-Legendre und MRW(Ma, Rouklin, Wandzura)
!gauss_legendre(Points, Art(1=Nodes, 2=Weights)
real *8, allocatable, dimension(:, :)                               ::                gauss_legendre
real *8, allocatable, dimension(:, :)                               ::                mrw


staebe_anzahl = 2
segmente_pro_stab=5

allocate(staebe(staebe_anzahl,2,3))
allocate(coord(staebe_anzahl, segmente_pro_stab, 3))
allocate(matrix( segmente_pro_stab*staebe_anzahl, segmente_pro_stab*staebe_anzahl) )
allocate(a(staebe_anzahl))
allocate(stroeme(segmente_pro_stab * staebe_anzahl))

!Integration wird mit 10 Schritten momentan fest geschrieben 
!TODO: Routine schreiben / finden, die die Nodes und Weights generisch berechnet
allocate(gauss_legendre(10, 2))
allocate(mrw(10, 2))
!Quelle: http://keisan.casio.com/exec/system/1280624821
gauss_legendre(1, 1) =	0.9931285991850949247861	
gauss_legendre(2, 1) =	0.9639719272779137912677	
gauss_legendre(3, 1) =	0.9122344282513259058678	
gauss_legendre(4, 1) =	0.8391169718222188233945	
gauss_legendre(5, 1) =	0.7463319064601507926143	
gauss_legendre(6, 1) =	0.6360536807265150254528	
gauss_legendre(7, 1) =	0.5108670019508270980044	
gauss_legendre(8, 1) =	0.3737060887154195606725	
gauss_legendre(9, 1) =	0.2277858511416450780805	
gauss_legendre(10, 1)= 	0.0765265211334973337546	

gauss_legendre(1, 2) =	0.017614007139152118312
gauss_legendre(2, 2) =	0.040601429800386941331
gauss_legendre(3, 2) =	0.06267204833410906357
gauss_legendre(4, 2) =	0.0832767415767047487248
gauss_legendre(5, 2) =	0.1019301198172404350368
gauss_legendre(6, 2) =	0.1181945319615184173124
gauss_legendre(7, 2) =	0.131688638449176626898
gauss_legendre(8, 2) =	0.1420961093183820513293
gauss_legendre(9, 2) =	0.149172986472603746788
gauss_legendre(10, 2)= 	0.1527533871307258506981

mrw(1, 1) =	0.482961710689630E-03 
mrw(2, 1) =	0.698862921431577E-02 
mrw(3, 1) =	0.326113965946776E-01 
mrw(4, 1) =	0.928257573891660E-01 
mrw(5, 1) =	0.198327256895404E+00 
mrw(6, 1) =	0.348880142979353E+00 
mrw(7, 1) =	0.530440555787956E+00 
mrw(8, 1) =	0.716764648511655E+00  
mrw(9, 1) =	0.875234557506234E+00  
mrw(10, 1) =	0.975245698684393E+00 

mrw(1, 2) =	0.183340007378985E-02 
mrw(2, 2) =	0.134531223459918E-01 
mrw(3, 2) =	0.404971943169583E-01 
mrw(4, 2) =	0.818223696589036E-01 
mrw(5, 2) =	0.129192342770138E+00 
mrw(6, 2) =	0.169545319547259E+00 
mrw(7, 2) =	0.189100216532996E+00 
mrw(8, 2) =	0.177965753961471E+00 
mrw(9, 2) =	0.133724770615462E+00 
mrw(10, 2) =	0.628655101770325E-01



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
                        coord(i, j, k) = (j-0.5) * b + staebe(i, 1, k)
                        write (10, *) coord(i,j,k)
                end do
        end do
end do

write (10, *) mrw


end

                

                
                        

