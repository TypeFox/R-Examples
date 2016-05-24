#' Update Intercept beta
#' 
#' \code{update_beta} updates the vector of band specific intercepts.
#' 
#' @param tms list of matrices whose rows are the triple (t,m,sigma) for each band.
#' @param a Amplitude estimates
#' @param rho vector of the current estimates of the phase
#' @param omega frequency
#' @export
#' @examples
#' test_data <- synthetic_multiband()
#' tms <- test_data$tms
#' a <- test_data$a
#' rho <- test_data$rho
#' omega <- test_data$omega
#' 
#' update_beta(tms,a,rho,omega)
update_beta <- function(tms,a,rho,omega) {
  B <- length(tms)
  beta <- double(B)
  for (b in 1:B) {
    band <- tms[[b]]
    m <- band[,2]
    s <- sin(omega*band[,1] + rho[b])
#    w <- 1/band[,3]
#    beta[b] <- sum(w*(m - A[b]*s))/sum(w)
    beta[b] <- sum( (m-a[b]*s)/(band[,3]**2))/sum(1/band[,3]**2)
  }
  return(beta)
}
