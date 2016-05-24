#' Update Beta0 parameter
#' 
#' \code{update_amplitude} solves a simple linear system of equations (rank 1 perturbation on a diagonal matrix) to update the beta0 estimates.
#' 
#' @param tms list of matrices whose rows are the triple (t,mu,sigma) for each band
#' @param a vector of the current amplitude estimates
#' @param rho vector of the current phase estimates
#' @param omega frequency
#' @param gamma nonnegative regularization parameter
#' @export
#' @examples
#' test_data <- synthetic_multiband()
#' B <- test_data$B
#' tms <- test_data$tms
#' beta <- test_data$beta
#' rho <- test_data$rho
#' omega <- test_data$omega
#' gamma <- 1
#' at <- rnorm(B)
#' at <- as.matrix(at/sqrt(sum(at**2)),ncol=1)
#' 
#' ## Check answer
#' a_next <- update_amplitude(tms,beta,rho,omega,at,gamma)
#' 
#' e <- double(B)
#' xi <- double(B)
#' for (b in 1:B) {
#'   nb <- length(tms[[b]][,1])
#'   s <- sin(omega*tms[[b]][,1] + rho[b])
#'   w <- 1/(tms[[b]][,3]**2)
#'   e[b] <- t(s)%*%(w*s) + gamma
#'   xi[b] <- t(s)%*%(w*(tms[[b]][,2]-beta[b]))
#' }
#' a_direct <- solve(diag(e)-gamma*at%*%t(at),xi)
#' norm(as.matrix(a_direct-a_next),'f')
update_beta_gamma <- function(tms,a,rho,omega,gamma) {
  B <- length(tms)
  e <- double(B)
  xi <- matrix(NA,B,1)
  for (b in 1:B) {
    band <- tms[[b]]
    nb <- length(band[,1])
    s <- sin(omega*band[,1] + rho[b])
    w <- 1/(band[,3]**2)
    iw <- band[,3]**2
    i <- rep(1,nb)
    xi[b] <- t(i)%*%((band[,2] - a[b]*s)/iw)
    e[b] <- t(i)%*%(i/iw) + gamma
  }
  at <- rep(1,B)/sqrt(B)
  eixi <- xi/e
  eiat <- at/e
  a <- eixi + gamma*(sum(eixi*at)/(1 - gamma*sum(at*eiat)))*eiat
  return(a)
}
