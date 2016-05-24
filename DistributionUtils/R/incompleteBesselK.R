### This file contains two versions of the incomplete Bessel K function
###
### incompleteBesselK calls the fortran function incompleteBesselK
### incompleteBesselKR is a pure R version of the incomplete Bessel K
### incompleteBesselKR also requires a number of additional R functions:
###     SSFcoef, combinatorial, GDENOM, and GNUM

incompleteBesselK <- function(x, y, nu, tol = (.Machine$double.eps)^(0.85),
                              nmax = 120) {

  KNu <- besselK(2*sqrt(x*y), nu)
  IBFOut <- .Fortran("incompleteBesselK",
                     as.double(x),
                     as.double(y),
                     as.double(nu),
                     as.double(tol),
                     as.integer(nmax),
                     as.double(KNu),
                     IBF = double(1),
                     status = integer(1)
                     )
  ## for debugging
  ##str(IBFOut)
  status <- IBFOut$status
  IBF <- IBFOut$IBF
  if(status == 1) warning("Maximum order exceeded\nResult may be unreliable")
  return(IBF)
}



incompleteBesselKR <- function(x, y, nu,
                               tol = (.Machine$double.eps)^(0.85),
                               nmax = 120) {
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
      ##cat("G[n] = ", G[n], "\tn = ", n, "\n")
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
      ##cat("G[n] = ", G[n], "\tn = ", n, "\n")
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


SSFcoef <- function(nmax, nu) {
  A <- matrix(rep(0, (nmax + 1)^2), ncol = nmax + 1)
  A[1,1] <- 1
  for(k in 2:(nmax+1)) {
    for(i in 2:k) {
      A[k, i] <- (-nu + i + k - 3)*A[k-1, i] + A[k-1, i-1]
    }
    A[k, k] <- 1
    A[k, 1] <- (-nu + k - 2)*A[k-1,1]
  }
  return(A)
}

combinatorial <- function(nu){
  Cnp <- numeric(0.5*(nu + 1)*(nu + 2))
  Cnp[1:min(3, 0.5*(nu + 1)*(nu + 2))] <- 1
  if (nu > 1) {
    for(n in 2:nu) {
      Cnp[n*(n + 1)/2 + 1] <- 1
      Cnp[n*(n + 1)/2 + n + 1] <- 1
      for(np in 1:(n - 1) ){
        Cnp[n*(n+1)/2 + np + 1] <-
          Cnp[n*(n-1)/2 + np] + Cnp[n*(n-1)/2 + np + 1]
      }
    }
  }
  return(Cnp)
}

GDENOM <- function(n, x, y, nu, An, nmax, Cnp){
  GN <- 0
  for (j in 0:n){
    terme <- 0
    for (i in 0:j){
      terme <- terme + An[j + 1, i + 1]*x^i
    }
    GN <- GN + Cnp[n*(n + 1)/2 + j + 1]*(-1/y)^j*terme
  }
  GN <- GN *(-x*y)^n*x^(nu + 1)*exp(x + y)
  return(GN)
}


GNUM <- function(n, x, y, nu, Am, An, nmax, Cnp, GM, GN){
  GM <- 0
  for (j in 1:n){
    terme <- 0
    for (k in 0:(j - 1)){
      termepr <- 0
      for (i in 0:k){
        termepr <- termepr + Am[k + 1, i + 1]*(-x)^i
      }
      terme <- terme + termepr*Cnp[j*(j - 1)/2 + k + 1]*(1/y)^k
    }
    GM <- GM + Cnp[n*(n + 1)/2 + j + 1]*(x*y)^j*GN[n - j + 1]*terme
  }
  GM <- GM *exp(-(x + y))*x^(-nu)/y
  return(GM)
}
