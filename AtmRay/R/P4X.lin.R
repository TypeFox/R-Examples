P4X.lin = function(x, y, zs, zr, ATM = CheckAtm.lin(list()), maxerror = 3){

#  dyn.load('AtmRay.so')
  L = .C('p4xlin', x = as.double(x), y = as.double(y), zs = as.double(zs) - ATM$z0, zr = as.double(zr) - ATM$z0, c0 = as.double(ATM$c0), gc = as.double(ATM$gc), wx0 = as.double(ATM$wx0), gwx = as.double(ATM$gwx), wy0 = as.double(ATM$wy0), gwy = as.double(ATM$gwy), rho0 = as.double(ATM$rho0), grho = as.double(ATM$grho), maxerror = as.double(maxerror), p = as.double(0), az = as.double(0), error = as.double(0))
  return(list(p = L$p, az = L$az, error = L$error))

}
