"censNID" <-
  function(y, n, cy, ctyp=c("left", "right"))
  {
    ctyp <- match.arg(ctyp,choices=c("left", "right"))
    L <- rep(-Inf, n)
    U <- rep(Inf, n)
    z <- numeric(n)
    m <- length(y)
    if (m==n) warning("length(y)==n, means that there is no censoring!")
    z[1:m] <- y
    z[(m+1):n] <- cy
    if (ctyp=="left")
      U[(m+1):n] <- cy else L[(m+1):n] <- cy
    cenNID(z, L, U)
  }