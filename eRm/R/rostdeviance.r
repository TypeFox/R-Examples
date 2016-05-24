rostdeviance <- function(object)
{
# Analysis of Deviance Table (Test against a saturated model)
# object... object of class ppar

#---------------saturated model---------------------
  X <- object$X
  N <- dim(X)[1]                     #number of subjects
  K <- dim(X)[2]                     #number of items
  x.ch <- apply(X,1,toString)        #response patters as string vectors
  nx <- as.vector(table(x.ch))       #pattern frequencies
  lsat <- sum(nx*(log(nx/N)))        #log-likelihood of saturated model (Rost, p.334)
  #npar.sat <- length(nx)
  npar.sat <- prod(apply(X, 2, max) + 1) - 1  #number of possible response patterns - 1
#------------end saturated model--------------------

  rv <- rowSums(X, na.rm = TRUE)                          #person raw scores
  lmml <- sum(table(rv)*log(table(rv)/N))+object$loglik.cml   #MML likelihood
  npar.mml <- dim(object$W)[2]        #+ length(table(rv)) ... not sure about that
  
  dev <- -2*(lmml - lsat)             #deviance
  df.chi <- npar.sat - npar.mml
  p.value <- 1-pchisq(dev,df.chi)
  result <- list(value = dev, df = df.chi, p.value = p.value)
  return(result)
}
