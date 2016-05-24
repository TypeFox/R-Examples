##
##  PURPOSE:   Handle init.alpha or init2.alpha argument of GLMM_MCMC function
##
##             THIS IS A HELP FUNCTION, NOT TO BE CALLED BY ORDINARY USERS
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    02/11/2011    
##
##  FUNCTIONS:  GLMM_MCMCinit.alpha
##
## ================================================================================================

## *************************************************************
## GLMM_MCMCinit.alpha
## *************************************************************
##
GLMM_MCMCinit.alpha <- function(init.alpha, lalpha, ialpha, number="")
{
  if (lalpha){
    if (missing(init.alpha)) init.alpha <- ialpha
    if (!is.numeric(init.alpha)) stop(paste("init", number, ".alpha must be numeric", sep=""))
    if (length(init.alpha) == 1) init.alpha <- rep(init.alpha, lalpha)
    if (length(init.alpha) != lalpha) stop(paste("init", number, ".alpha must be of length ", lalpha, sep=""))
    if (is.null(names(init.alpha))) names(init.alpha) <- paste("alpha", 1:lalpha, sep="")
  }else{
    init.alpha <- 0
  }  

  return(init.alpha)
}  
