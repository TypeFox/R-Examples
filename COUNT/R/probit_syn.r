# Generic synthetic probit regression  24Sept 2010
# Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
# Hilbe, Logistic Regression Models, Chapman & Hall/CRC
require(MASS)                  # probit_syn.r
probit_syn  <- function(nobs=50000, d = 1, xv = c(1, 0.5, -1.5))  {
   p <- length(xv) - 1
   X <- cbind(1, matrix(rnorm(nobs * p), ncol = p))
   xb <- X %*% xv
   pxb <- pnorm(xb)                      
   py  <- rbinom(nobs, size = d, prob =pxb)   
   dpy <- d - py
   out <- data.frame(cbind(cbind(py,dpy), X[,-1]))
   names(out) <- c("py","dpy", paste("x", 1:p, sep=""))
   return(out)
}


