"rho.weight" <-
function(svals, ips = 1, xk = 1.06){
  n <- length(svals)
  fvals <- double(n)
  storage.mode(svals) <- "double"
  f.res <- .Fortran("srrhoamm",
    n = as.integer(n),
    svals = svals,
    fvals = fvals,
    as.integer(ips),
    as.double(xk))
  f.res$fvals}

