# Generic synthetic NB-C model
#  Table 10.9, Hilbe, J.M., Negative Binomial Regression, 2nd ed.
#    Cambridge Univ Press       nbc_syn.r
require(MASS)
    nbc_syn  <- function(nobs = 50000,
                        alpha = 1.15,
                        xv = c(-1.5, -1.25, -.1))  {
      q <- length(xv) - 1
      X <- cbind(1, matrix(runif(nobs * q), ncol = q))
      xb <- X %*% xv
      a <- alpha
      mu <- 1/((exp(-xb)-1)*a)
      p <- 1/(1+a*mu)
      r <- 1/a
      nbcy <- rnbinom(nobs, size=r, prob = p)
      out <- data.frame(cbind(nbcy, X[,-1]))
      names(out) <- c("nbcy", paste("x", 1:q, sep=""))
      return(out)
    }


