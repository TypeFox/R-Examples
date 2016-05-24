# Generic synthetic binomial logistic regression  24Sept 2010
# Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
# Hilbe, Logistic Regression Models, Chapman & Hall/CRC
require(MASS)                                     # logit_syn.r
logit_syn  <- function(nobs=50000, d = 1, xv = c(1, 0.5, -1.5))  {
   p <- length(xv) - 1
   X <- cbind(1, matrix(rnorm(nobs * p), ncol = p))
   xb <- X %*% xv
   exb <- 1/(1+exp(-xb))                      
   by  <- rbinom(nobs, size = d, prob =exb)   
   dby <- d - by
   out <- data.frame(cbind(cbind(by,dby), X[,-1]))
   names(out) <- c("by","dby", paste("x", 1:p, sep=""))
   return(out)
}



