#############################################################
#                                                           #
#	residualsAnscombe function                          #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: October, 15, 2012                             #
#	Version: 0.1-4                                      #
#                                                           #
#	Copyright (C) 2012 Claudio Agostinelli              #
#                                                           #
#############################################################

residualsAnscombe <- function(y, mu, family, ...) {
  family <- family$family
  ny <- length(y)
  ## Binomial
  if (family=='binomial') {
    ## nel caso di binomuale mu deve essere la probabilita' di successo
    if (any(mu > 1 | mu < 0))
      stop("in the binomial family 'mu' must be the probability of success")
###    resid <- (3*y^(1/3)*((y^(1/3))/2 - 1) - 3*mu^(1/3)*((mu^(1/3))/2 - 1))/(mu^(1/6)*(1-mu^(-1/3))*sqrt(1-mu))
    ## usata formula in COX and SNELL 1968 invece che di quella in gill+1998.pdf
    ## che appare commentata qui sopra
    if (NCOL(y)) {
      resid <- beta(2/3,2/3)*(pbeta(y, 2/3, 2/3) - pbeta(mu-(1-2*mu)/6, 2/3, 2/3))/(mu^(1/6)*(1-mu)^(1/6))
    } else {
      m <- apply(y, 1, sum)
      y <- y[,1]
      resid <- beta(2/3,2/3)*(pbeta(y/m, 2/3, 2/3) - pbeta(mu-(1-2*mu)/(6*m), 2/3, 2/3))/((mu^(1/6)*(1-mu)^(1/6))/sqrt(m))
    }
  ## Poisson    
  } else if (family=='poisson') {
    ## usata formula in COX and SNELL 1968
    mu1 <- ifelse(mu>1/6, mu-1/6, 0)
    resid <- (3/2*(y^(2/3) - (mu1)^(2/3)))/(mu^(1/6))
  ## Gamma
  } else if (family=='Gamma') {
    ## usata formula in McCullagh & Nelder 1989      
    resid <- 3*(y^(1/3) - mu^(1/3))/(mu^(1/3))
  } else if (family=='gaussian') {
  ## Normal I have to check it!!!!!!!!!!!
    resid <- y - mu
  } else if (family=='inverse.gaussian') {
    ## usata formula in McCullagh & Nelder 1989    
    resid <- (log(y)-log(mu))/sqrt(mu)
  } else {
     resid <- rep(NA, ny)
     warning('the function is not implemented for the family', family$family, '\n')
  }
  return(resid)
}
