################################################
#### AUTHOR:     Arnost Komarek             ####
####             (2005)                     ####
####                                        ####
#### FILE:       bayessurvreg2.priorBeta.R  ####
####                                        ####
#### FUNCTIONS:  bayessurvreg2.priorBeta    ####
####
#### 24/11/2008: bug appearing when design$nX == 0 & design$randomInt == TRUE fixed
####
################################################

### ======================================
### bayessurvreg2.priorBeta
### ======================================
## Manipulation with the prior specification for the regression parameters
## (and means of random effects)
## - version for bayessurvreg2 (AFT model with G-spline error)
##
## 24/01/2005
## =========================================================================
bayessurvreg2.priorBeta <- function(prior.beta, init, design)
{  
  if (design$nX){
    ## Initials
    if(length(init) == 0) ininit <- "arnost"
    else                  ininit <- names(init)    
    tmp <- match("beta", ininit, nomatch=NA)
    if(is.na(tmp)){     
      init.beta <- rep(0, design$nX)
    }
    else{
      if (length(init$beta) < design$nX) stop("Incorrect init$beta parameter supplied.")
      init.beta <- init$beta[1:design$nX]
    }
    if (sum(is.na(init.beta))) stop("Incorrect init$beta parameter supplied.")

    
    ## Specification of priors    
    if(length(prior.beta) == 0) inprior <- "arnost"
    else                        inprior <- names(prior.beta)
    tmp <- match("mean.prior", inprior, nomatch=NA)
    if(is.na(tmp)) stop("Prior means for betas must be specified.")
    tmp <- match("var.prior", inprior, nomatch=NA)
    if(is.na(tmp)) stop("Prior vars for betas must be specified.")
    if (length(prior.beta$mean.prior) != design$nX) stop("Incorrect length of a vector of prior means for betas.")
    if (length(prior.beta$var.prior) != design$nX) stop("Incorrect length of a vector of prior vars for betas.")  
    if (sum(is.na(prior.beta$mean.prior))) stop("Prior means for betas must not be missing.")
    if (sum(is.na(prior.beta$var.prior))) stop("Prior vars for betas must not be missing.")  
    if (sum(prior.beta$var.prior <= 0)) stop("Prior vars for betas must be all positive.")
    mean.prior <- prior.beta$mean.prior
    var.prior <- prior.beta$var.prior

    
    ## Input parameters for classBetaGamma constructor
    ngamma <- design$nrandom - 1*(design$randomInt)
    nFixed <- design$nX - ngamma
    parmI <- c(design$nX, nFixed, ngamma,  1*design$randomInt, design$indb)
    parmD <- c(init.beta, mean.prior, var.prior)
    names(parmI) <- c("nbeta", "nFixed", "ngamma", "randomIntcpt", paste("indbA", 1:design$nX, sep=""))
    names(parmD) <- c(paste("beta", 1:design$nX, sep=""), paste("mean.beta", 1:design$nX, sep=""), paste("var.beta", 1:design$nX, sep=""))
  }
  else{
    if (!is.null(design$randomInt)){
      init.beta <- numeric(0)    
      mean.prior <- numeric(0)
      var.prior <- numeric(0)
      if (is.null(design$indb)) parmI <- c(0, 0, 0, 1*design$randomInt, 0)
      else                      parmI <- c(0, 0, 0, 1*design$randomInt, design$indb)
      parmD <- c(0, 0, 0)
      names(parmI) <- c("nbeta", "nFixed", "ngamma", "randomIntcpt", "indbA")
      names(parmD) <- c("beta", "mean.beta", "var.beta")    
    }
    else{
      init.beta <- numeric(0)    
      mean.prior <- numeric(0)
      var.prior <- numeric(0)
      parmI <- c(0, 0, 0, 0, 0)
      parmD <- c(0, 0, 0)
      names(parmI) <- c("nbeta", "nFixed", "ngamma", "randomIntcpt", "indbA")
      names(parmD) <- c("beta", "mean.beta", "var.beta")    
    }      
  }    

  toreturn <- list(parmI=parmI, parmD=parmD)
  attr(toreturn, "init") <- init.beta
  attr(toreturn, "prior.beta") <- list(mean.prior=mean.prior, var.prior=var.prior)

  return(toreturn)
} 


