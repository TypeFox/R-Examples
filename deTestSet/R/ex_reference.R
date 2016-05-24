

### ============================================================================
### Finds exact, reference solution for certain problem
### ============================================================================
reference <- function(name=c("andrews","beam","caraxis","crank","E5","emep","fekete",
    "vdpol","hires","nand","orego","pleiades","pollution","ring","rober","transistor",
    "tube","twobit","wheelset")) {
  switch(name,
    "vdpol" = c(0.1706167732170469e1, -0.8928097010248125e-3),
    
    "orego" = c(0.1000814870318523e1, 0.1228178521549917e4, 
              0.1320554942846706e3),
              
    "hires" = c(y1=0.7371312573325668e-3,y2=0.1442485726316185e-3,
               y3=0.5888729740967575e-4, y4=0.1175651343283149e-2,
               y5=0.2386356198831331e-2, y6=0.6238968252742796e-2,
               y7=0.2849998395185769e-2, y8=0.2850001604814231e-2),
               
    "rober" = c(0.2083340149701255e-7, 0.8333360770334713e-13, 
              0.9999999791665050),
              
    "E5" = c(0.1152903278711829e-290,0.8867655517642120e-22,
          0.8854814626268838e-22,0.),
          
    "transistor" = c(-0.5562145012262709e-2,0.3006522471903042e1,
                      0.2849958788608128e1, 0.2926422536206241e1,
                      0.2704617865010554e1,0.2761837778393145e1,
                      0.4770927631616772e1,0.1236995868091548e1),
                      
    "caraxis" = c(0.493455784275402809122e-1,0.496989460230171153861,
                  0.104174252488542151681e1,0.373911027265361256927,
                  -0.770583684040972357970e-1, 0.744686658723778553466e-2,
                  0.175568157537232222276e-1,0.770341043779251976443,
                 -0.473688659084893324729e-2,-0.110468033125734368808e-2), 
    
    "wheelset" = c(0.86355386965811e-2,0.13038281022727e-4,
                   -0.93635784016818e-4,-0.13642299804033e-1,
                   0.15292895005422e-2,-0.76985374142666e-1,
                   -0.25151106429207e-3,0.20541188079539e-2,
                   -0.23904837703692,-0.13633468454173e-1,
                   -0.24421377661131,-0.33666751972196e-3,
                   -0.15949425684022, 0.37839614386969e-3,
                   0.14173214964613,-0.10124044903201e-1,
                   -0.56285630573753e-2) ,
    "pleiades" =  .Fortran("pleiasoln", neqn=as.integer(28),y = as.double(rep(0.,28))  )$y,
    "beam" =      .Fortran("beamsoln", neqn=as.integer(80),y = as.double(rep(0.,80))   )$y,
    "andrews" =   .Fortran("andsoln", neqn=as.integer(27),y = as.double(rep(0.,27))    )$y,
    "crank" =     .Fortran("cranksoln", neqn=as.integer(24),y = as.double(rep(0.,24))  )$y,
    "fekete" =    .Fortran("feksoln", neqn=as.integer(160),y = as.double(rep(0.,160))  )$y,
    "emep"    =   .Fortran("emepsoln", neqn=as.integer(66),y = as.double(rep(0.,66))   )$y,
    "nand"  =     .Fortran("nandsoln", neqn=as.integer(14),y = as.double(rep(0.,14))   )$y,
    "pollution" = .Fortran("polsoln", neqn=as.integer(20),y = as.double(rep(0.,20))    )$y,
    "ring" =      .Fortran("ringsoln", neqn=as.integer(15),y = as.double(rep(0.,15))   )$y,
    "tube" =      .Fortran("tubesoln", neqn=as.integer(49),y = as.double(rep(0.,49))   )$y,
    "twobit" =    .Fortran("twobsoln", neqn=as.integer(350),y = as.double(rep(0.,350)) )$y ,
    stop(paste("cannot find ",name))
    
              )

}

