A2P.lin = function(angle, z, az, ATM){
  ceff = as.vector(ATM$c0 + c(ATM$wx0, ATM$wy0) %*% c(sin(az*pi/180), cos(az*pi/180)) + (z - ATM$z0) * (ATM$gc + c(ATM$gwx, ATM$gwy) %*% c(sin(az*pi/180), cos(az*pi/180))))
  p = sin(angle * pi/180)/ceff
  return(p)
}


