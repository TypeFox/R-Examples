MakeArrivals.lin = function(xs, ys, zs, xr, yr, zr, dt, nt, timing, ATM= CheckAtm.lin(list())){


  # check input
  ATM = CheckAtm.lin(ATM)
  
  if(length(xs) != length(ys) | length(xs) != length(zs)| length(xs) != length(timing)){
    nsrc = max(length(xs), length(ys), length(zs), length(timing))
    if(length(xs) == 1){
      xs = rep(xs, nsrc)
    }
    if(length(ys) == 1){
      ys = rep(ys, nsrc)
    }
    if(length(zs) == 1){
      zs = rep(zs, nsrc)
    }
    if(length(timing)==1){
      timing = rep(timing, nsrc)
    }
    
    if(length(xs) != length(ys) | length(xs) != length(zs) | length(xs) != length(timing)){
      stop('xs, ys, zs, and timing must be the same length')
    }
  }
  
  if(length(xr) != length(yr) | length(xr) != length(zr)){
    nrec = max(length(xr), length(yr), length(zr))
    if(length(xr) == 1){
      xr = rep(xr, nrec)
    }
    if(length(yr) == 1){
      yr = rep(yr, nrec)
    }
    if(length(zr) == 1){
      zr = rep(zr, nrec)
    }
    if(length(xr) != length(yr) | length(xr) != length(zr)){
      stop('xr, yr, and zr must be the same length')
    }
  }

#  dyn.load('AtmRay.so')
  
  nsrc = length(xs)
  nrec = length(xr)
  P = rep(0, nt * nrec)
#  browser()

  L = .C('makearrivals', xs = as.double(xs), ys = as.double(ys), zs = as.double(zs) - ATM$z0, xr = as.double(xr), yr = as.double(yr), zr = as.double(zr) - ATM$z0, dt = as.double(dt), nt = as.integer(nt), timing = as.double(timing), c0 = as.double(ATM$c0), gc = as.double(ATM$gc), wx0 = as.double(ATM$wx0), gwx = as.double(ATM$gwx), wy0 = as.double(ATM$wy0), gwy = as.double(ATM$gwy), rho0 = as.double(ATM$rho0), grho = as.double(ATM$grho), nsrc = as.integer(nsrc), nrec = as.integer(nrec), P = as.double(P), NAOK = TRUE)
#  browser()
  return(matrix(L$P, nt, nrec))
}
