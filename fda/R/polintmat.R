polintmat <- function(xa, ya, x) {
#  Polynomial extrapolion for a converging sequence
#  YA is an 3-D array with 1st D same as XA
  n     <- length(xa)
  dimya <- dim(as.array(ya))
  if (length(dimya) == 1) ya <- array(ya,c(dimya[1],1,1))
  if (length(dimya) == 2) ya <- array(ya,c(dimya[1],dimya[2],1))
  if (dimya[1] != n)      stop('First dimension of YA must match XA')
  difx <- xa - x
  absxmxa <- abs(difx)
  ns <- min((1:n)[absxmxa == min(absxmxa)])
  cs <- ds <- ya
  y  <- ya[ns,,]
  ns <- ns - 1
  for (m in 1:(n-1)) {
    for (i in 1:(n-m)) {
      ho      <- difx[i]
      hp      <- difx[i+m]
      w       <- (cs[i+1,,] - ds[i,,])/(ho - hp)
      ds[i,,] <- hp*w
      cs[i,,] <- ho*w
    }
    if (2*ns < n-m) {
      dy <- cs[ns+1,,]
    } else {
      dy <- ds[ns,,]
      ns <- ns - 1
    }
    y <- y + dy
  }
  return( list(y, dy) )
}
