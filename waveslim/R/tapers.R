dpss.taper <- function(n, k, nw = 4, nmax = 2^(ceiling(log(n,2)))) {
  if(n > nmax)
    stop("length of taper is greater than nmax")
  w <- nw/n
  if(w > 0.5)
    stop("half-bandwidth parameter (w) is greater than 1/2")
  if(k <= 0)
    stop("positive dpss order (k) required")
  v <- matrix(0, nrow = nmax, ncol = (k + 1))
  storage.mode(v) <- "double"
  out <- .Fortran("dpss",
                  nmax = as.integer(nmax),
                  kmax = as.integer(k),
                  n = as.integer(n),
                  w = as.double(w),
                  v = v,
                  sig = double(k + 1),
                  totit = integer(1),
                  sines = double(n),
                  vold = double(n),
                  u = double(n),
                  scr1 = double(n),
                  ifault = integer(1),
                  PACKAGE="waveslim")
  ##list(v = out$v[1:n, 1:k], eigen = out$sig[-1] + 1, iter = 
  ##     out$totiTRUE, n = out$n, w = out$w, ifault = out$ifault)
  return(out$v[1:n, 1:k])
}

sine.taper <- function(n, k) {
  tapers <- NULL
  for(i in 1:k)
    tapers <- cbind(tapers, sqrt(2/(n+1)) * sin((pi*i*1:n)/(n+1)))
  return(tapers)
}
