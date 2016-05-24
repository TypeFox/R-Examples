ldafitfun <- function(x, covt, mt1,mt2, prior)
  ## 14/10/10 IH
  ##
  ## Internal function of sprmsDA, prmsDA and predict functions
{
  d1 <- matrix(x, nrow=1)%*%solve.qr(qr(covt))%*%matrix(mt1, ncol=1) -1/2*matrix(mt1, nrow=1)%*%solve.qr(qr(covt))%*%matrix(mt1, ncol=1) + log(prior[1])
  d2 <- matrix(x, nrow=1)%*%solve.qr(qr(covt))%*%matrix(mt2, ncol=1) -1/2*matrix(mt2, nrow=1)%*%solve.qr(qr(covt))%*%matrix(mt2, ncol=1) + log(prior[2])
  return(c(d1,d2))
}
