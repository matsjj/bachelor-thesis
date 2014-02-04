!Integration wird mit 10 Schritten momentan fest geschrieben
!TODO: Routine schreiben / finden, die die Nodes und Weights generisch berechnet
allocate(gauss_legendre_points(10))
allocate(mrw_points(10))
allocate(gauss_legendre_weights(10))
allocate(mrw_weights(10))
!Quelle: http://keisan.casio.com/exec/system/1280624821
gauss_legendre_points(1) =  0.9931285991850949247861
gauss_legendre_points(2) =  0.9639719272779137912677
gauss_legendre_points(3) =  0.9122344282513259058678
gauss_legendre_points(4) =  0.8391169718222188233945
gauss_legendre_points(5) =  0.7463319064601507926143
gauss_legendre_points(6) =  0.6360536807265150254528
gauss_legendre_points(7) =  0.5108670019508270980044
gauss_legendre_points(8) =  0.3737060887154195606725
gauss_legendre_points(9) =  0.2277858511416450780805
gauss_legendre_points(10)=  0.0765265211334973337546

gauss_legendre_weights(1) =  0.017614007139152118312
gauss_legendre_weights(2) =  0.040601429800386941331
gauss_legendre_weights(3) =  0.06267204833410906357
gauss_legendre_weights(4) =  0.0832767415767047487248
gauss_legendre_weights(5) =  0.1019301198172404350368
gauss_legendre_weights(6) =  0.1181945319615184173124
gauss_legendre_weights(7) =  0.131688638449176626898
gauss_legendre_weights(8) =  0.1420961093183820513293
gauss_legendre_weights(9) =  0.149172986472603746788
gauss_legendre_weights(10)=  0.1527533871307258506981

mrw_points(1) = 0.482961710689630E-03
mrw_points(2) = 0.698862921431577E-02
mrw_points(3) = 0.326113965946776E-01
mrw_points(4) = 0.928257573891660E-01
mrw_points(5) = 0.198327256895404E+00
mrw_points(6) = 0.348880142979353E+00
mrw_points(7) = 0.530440555787956E+00
mrw_points(8) = 0.716764648511655E+00
mrw_points(9) = 0.875234557506234E+00
mrw_points(10) =    0.975245698684393E+00

mrw_weights(1) = 0.183340007378985E-02
mrw_weights(2) = 0.134531223459918E-01
mrw_weights(3) = 0.404971943169583E-01
mrw_weights(4) = 0.818223696589036E-01
mrw_weights(5) = 0.129192342770138E+00
mrw_weights(6) = 0.169545319547259E+00
mrw_weights(7) = 0.189100216532996E+00
mrw_weights(8) = 0.177965753961471E+00
mrw_weights(9) = 0.133724770615462E+00
mrw_weights(10) =    0.628655101770325E-01
