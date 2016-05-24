t_vol <-
function(vol, m)
 {
 d0   = 0.01
 dmax = 130
 V0 <- pi*(d0^3)/6
 Vmax <- pi*(dmax^3)/6# = 1150347
 s     <- 31*m
 ttt <- (-1/m)*log( 1 - (m/s)*log(vol/V0) ) / 365# time to reach volume vol in years
 return(ttt)
 }
