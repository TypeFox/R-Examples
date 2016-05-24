#' Update Amplitude parameter
#' 
#' \code{update_amplitude} solves a simple linear system of equations (rank 1 perturbation on a diagonal matrix) to update the amplitude estimates.
#' 
#' @param tms list of matrices whose rows are the triple (t,mu,sigma) for each band
#' @param beta vector of the current intercept estimates
#' @param rho vector of the current phase estimates
#' @param omega frequency
#' @param at prior vector
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
update_amplitude <- function(tms,beta,rho,omega,at,gamma) {
  B <- length(tms)
  e <- double(B)
  xi <- matrix(NA,B,1)
  for (b in 1:B) {
    band <- tms[[b]]
    nb <- length(band[,1])
    s <- sin(omega*band[,1] + rho[b])
    w <- 1/(band[,3]**2)
    iw <- band[,3]**2
    xi[b] <- t(s)%*%((band[,2] - beta[b])/iw)
    e[b] <- t(s)%*%(s/iw) + gamma
  }
  eixi <- xi/e
  eiat <- at/e
  a <- eixi + gamma*(sum(eixi*at)/(1 - gamma*sum(at*eiat)))*eiat
  return(a)
}
