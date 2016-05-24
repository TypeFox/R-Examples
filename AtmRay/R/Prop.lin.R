Prop.lin = function(p, az, zs, zr, ATM = CheckAtm.lin(list())){
  
  ATM = CheckAtm.lin(ATM)
  L = .C('proplin', p = as.double(p), az = as.double(az), zs = as.double(zs) - ATM$z0, zr = as.double(zr) - ATM$z0, c0 = as.double(ATM$c0), gc = as.double(ATM$gc), wx0 = as.double(ATM$wx0), gwx = as.double(ATM$gwx), wy0 = as.double(ATM$wy0), gwy = as.double(ATM$gwy), rho0 = as.double(ATM$rho0), grho = as.double(ATM$grho), x = as.double(0), y = as.double(0), t = as.double(0), A = as.double(0))
  return(list(x = L$x, y = L$y, t = L$t, A = L$A, p = L$p))
}
