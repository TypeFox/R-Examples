
# Generic synthetic NB1 model
# In Hilbe, J.M., Negative Binomial Regression, 2nd ed, Cambridge Univ Press
require(MASS)                     # nb1_syn.r
nb1_syn  <- function(nobs = 50000,
                    delta = 1,
                    xv = c(1, 0.75, -1.25))  {
  p <- length(xv) - 1
  X <- cbind(1, matrix(rnorm(nobs * p), ncol = p))
  xb <- X %*% xv
  d <- delta                  
  exb <- exp(xb)   
  idelta <- (1/delta)*exb            
  xg <- rgamma(n = nobs, shape = idelta, rate = idelta)  
  xbg <- exb*xg
  nb1y <- rpois(nobs, xbg)     
  out <- data.frame(cbind(nb1y, X[,-1]))
  names(out) <- c("nb1y", paste("x", 1:p, sep=""))
  return(out)
}
