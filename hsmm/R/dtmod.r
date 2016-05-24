######################################################
# Density functions of the conditional distributions #
######################################################

# t distribution
# --------------
dtmod <- function(x, mu = 0, sigma = 1, nu = 30, log = FALSE){
  if (log == TRUE){
    x <- log(x)
    }
  nenner1  <- try(sigma*beta(1/2, nu/2))
  zaehler1 <- try(nu^(-1/2))
  nenner2  <- try(nu*sigma^2)
  zaehler2 <- try((x-mu)^2)
  dtmod    <- try(zaehler1/nenner1*(1+zaehler2/nenner2)^(-1*(nu+1)/2))
  return (dtmod)
  }
