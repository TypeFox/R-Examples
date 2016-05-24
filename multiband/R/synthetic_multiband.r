#' Generate Synthetic Multiband Data
#' 
#' \code{synethetic_multiband} generates multiband data of a given fequency.
#' 
#' @param B number of bands
#' @param omega frequency common to all bands
#' @param beta mean intercept
#' @param beta_sd standard deviation around mean intercept
#' @param a mean amplitude
#' @param a_sd standard deviation around mean amplitude
#' @param rho mean phase
#' @param rho_sd standard deviation around mean phase
#' @param n_max maximum number of time points per band
#' @param n_min minimum number of time points per band
#' @param seed initialization parameter for random number generator
#' @export
#' @examples
#' tms <- synthetic_multiband()
synthetic_multiband <- function(B=5,omega=0.2,beta=0,beta_sd=0.1,a=2,a_sd=0.25,rho=0,rho_sd=0.1,n_max=1e3,n_min=1e2,seed=12345) {
  set.seed(seed)
  tms <- vector(mode="list",length=B)
  nb <- sample(x=n_min:n_max,size=B,replace=TRUE)
  betab <- beta + sqrt(3)*beta_sd*runif(n=B,min=-1,max=1)
  Ab <- a + sqrt(3)*a_sd*runif(n=B,min=-1,max=1)
  rhob <- rho + sqrt(3)*rho_sd*runif(n=B,min=-1,max=1)
  for (b in 1:B) {
    t <- cumsum(runif(nb[b]))
    sigma <- runif(nb[b])
    m <- betab[b] + Ab[b]*sin(omega*t + rhob[b]) + sigma*rnorm(nb[b])
    tms[[b]] <- cbind(t,cbind(m,sigma))
  }
  return(list(tms=tms,B=B,omega=omega,beta=betab,a=Ab,rho=rhob))
}
