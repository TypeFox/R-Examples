#' Check Gradients
#' 
#' \code{gradient_check} computes the block coordinate gradients.
#' 
#' @param tms list of matrices whose rows are the triple (t,m,sigma) for each band
#' @param beta initial intercept estimates
#' @param A initial amplitude estimates
#' @param at prior vector
#' @param rho Initial phase estimates
#' @param omega frequency
#' @param gamma1 nonnegative regularization parameter for shrinking amplitudes
#' @param gamma2 nonnegative regularization parameter for shrinking phases
#' @export
#' 
gradient_check <- function(tms,beta,A,at,rho,omega,gamma1,gamma2) {
  B <- length(tms)

  gbeta <- double(B)
  grho <- double(B)
  
  e <- double(B)
  xi <- matrix(NA,B,1)
  for (b in 1:B) {
    band <- tms[[b]]
    n <- nrow(band)
    s <- sin(omega*band[,1] + rho[b])
    c <- cos(omega*band[,1] + rho[b])
    w <- 1/(band[,3]**2)
    xi[b] <- t(s)%*%(w*(band[,2] - beta[b]))
    e[b] <- t(s)%*%(w*s) + gamma1
    gbeta[b] <- 2*beta[b]*sum(w) - 2*sum(w*(band[,2]-A[b]*s))
    grho[b] <- A[b]*t(beta[b] + A[b]*s - band[,2])%*%(w*c) + gamma2*rho[b] - (gamma2/B)*sum(rho)
  }
  gA <- c(e*A - gamma1*sum(A*at)*at - xi)
  
  return(list(gbeta=gbeta,gA=gA,grho=grho))
}