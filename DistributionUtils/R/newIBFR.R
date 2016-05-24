newIBFR <- function(x, y, nu,
                               tol = (.Machine$double.eps)^(0.85),
                               nmax = 120) {
    cat("x = ", x, " y = ", y, "\n")
  Am <- matrix(rep(0, (nmax + 1)^2), ncol = nmax + 1)
  An <- matrix(rep(0, (nmax + 1)^2), ncol = nmax + 1)
  Cnp <- numeric((nmax + 1)*(nmax + 2)/2)
  G <- numeric(nmax)
  GN <- numeric(nmax + 1)
  GM <- numeric(nmax)
  Cnp <- combinatorial(nmax)
  if(x >= y) {
    Am <- SSFcoef(nmax, nu - 1)
    An <- SSFcoef(nmax, -nu -1)
    GN[1] <- GDENOM(0, x, y, nu, An, nmax, Cnp)
    GN[2] <- GDENOM(1, x, y, nu, An, nmax, Cnp)
    GM[1] <- GNUM(1, x, y, nu, Am, An, nmax, Cnp, GM ,GN)
    G[1] <- x^nu*GM[1]/GN[2]
    for (n in (2:nmax)){
      GN[n + 1] <- GDENOM(n, x, y, nu, An, nmax, Cnp)
      GM[n] <- GNUM(n, x, y, nu, Am, An, nmax, Cnp, GM, GN)
      G[n] <- x^(nu)*GM[n]/GN[n+1]
      cat("GM[n] = ", GM[n], "\n")
      cat("G[n] = ", G[n], "abs(G[n] - G[n-1]) = ",
          abs(G[n] - G[n-1]), "\tn = ", n, "\n")
      if (abs(G[n] - G[n-1])<= tol) {
        break
      }
    }
  }
  if(y > x) {
    BK <- besselK(2*sqrt(x*y), nu)
    Am <- SSFcoef(nmax, -nu - 1)
    An <- SSFcoef(nmax, nu - 1)
    GN[1] <- GDENOM(0, y, x, -nu, An, nmax, Cnp)
    GN[2] <- GDENOM(1, y, x, -nu, An, nmax, Cnp)
    GM[1] <- GNUM(1, y, x, -nu, Am, An, nmax, Cnp, GM, GN)
    G[1] <- (y^(-nu))*GM[1]/GN[2]
    for (n in (2:nmax)){
      GN[n + 1] <- GDENOM(n, y, x, -nu, An, nmax, Cnp)
      GM[n] <- GNUM(n, y, x, -nu, Am, An, nmax, Cnp, GM, GN)
      G[n] <- (y^(-nu))*GM[n]/GN[n + 1]
      cat("GM[n] = ", GM[n], "\n")
      cat("G[n] = ", G[n], "abs(G[n] - G[n-1]) = ",
          abs(G[n] - G[n-1]), "\tn = ", n, "\n")
      if (abs(G[n] - G[n-1])<= tol) {
        G[n] <- 2*((x/y)^(nu/2))*BK - G[n]
        break
      }
    }
  }
  if(n == nmax) warning("Maximum G transformation order reached")
  IBF <- G[n]
  return(IBF)
}

