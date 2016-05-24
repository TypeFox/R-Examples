#' Penalized negative log likelihood
#' 
#' \code{pnll} computes the penalized negative log likelihood
#' 
#' @param tms list of matrices whose rows are the triple (t,m,sigma) for each band.
#' @param beta vector of intercepts
#' @param a vector of amplitudes
#' @param atilde prior vector
#' @param rho vector of phases.
#' @param omega frequency
#' @param gamma1 nonnegative regularization parameter for shrinking amplitudes
#' @param gamma2 nonnegative regularization parameter for shrinking phases
#' @export
pnll <- function(tms,beta,a,atilde,rho,omega,gamma1,gamma2) {
  B <- length(tms)
  f <- 0
  for (b in 1:B) {
    band <- tms[[b]]
#    n <- nrow(band)
    w <- 1/(band[,3]**2)
    s <- sin(omega*band[,1] + rho[b])
    m <- band[,2]
    r <- m - beta[b] - a[b]*s
    RSS <- 0.5*t(r)%*%(w*r)
    f <- f + RSS
  }
#  penalty <- 0
#  for (b in 1:(B-1)) {
#    for (bb in (b+1):B) {
#      penalty <- penalty + gamma1*sum((a[b] - a[bb])**2)
#      penalty <- penalty + gamma2*sum((rho[b] - rho[bb])**2)
#    }
#  }
  penalty <- gamma2*(sum(rho**2) - (1/B)*(sum(rho)**2))
  penalty <- penalty + gamma1*(sum(a**2) - (sum(atilde*a))**2)
  return(f + 0.5*penalty)
}
