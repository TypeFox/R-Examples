# Generic synthetic negative binomial (NB2) data and model
#  Hilbe, J.M Negative Binomial Regression, 2nd ed, Cambridge University Press   
nb2_syn  <- function(nobs = 50000, off = 0,
                    alpha = 1,
                    xv = c(1, 0.75, -1.5))  {
  p <- length(xv) - 1
  X <- cbind(1, matrix(rnorm(nobs * p), ncol = p))
  xb <- X %*% xv
  a <- alpha                  
  ia <- 1/a                   
  exb <- exp(xb + off)        
  xg <- rgamma(n = nobs, shape = a, rate = a)
  xbg <-exb*xg
  nby <- rpois(nobs, xbg) 
  out <- data.frame(cbind(nby, X[,-1]))
  names(out) <- c("nby", paste("x", 1:p, sep=""))
  return(out)
}

