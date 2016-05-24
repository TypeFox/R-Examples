##
##  PURPOSE:   Handle init.eps or init2.eps argument of GLMM_MCMC function
##
##             THIS IS A HELP FUNCTION, NOT TO BE CALLED BY ORDINARY USERS
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    02/11/2011    
##
##  FUNCTIONS:  GLMM_MCMCinit.eps
##
## ================================================================================================

## *************************************************************
## GLMM_MCMCinit.eps
## *************************************************************
##
GLMM_MCMCinit.eps <- function(init.eps, prior.eps, Rc, isigma, number="")
{
  if (Rc){
    if (missing(init.eps)) init.eps <- list()
    if (!is.list(init.eps)) stop(paste("init", number, ".eps must be a list", sep=""))    
    
    ininit.eps <- names(init.eps)
    ieps.sigma    <- match("sigma", ininit.eps, nomatch=NA)
    ieps.gammaInv <- match("gammaInv", ininit.eps, nomatch=NA)    

    ##### init.eps:  sigma
    ##### -----------------------------------------------
    if (is.na(ieps.sigma)) init.eps$sigma <- isigma
    if (length(init.eps$sigma) == 1) init.eps$sigma <- rep(init.eps$sigma, Rc)
    if (length(init.eps$sigma) != Rc) stop(paste("init", number, ".eps$sigma must be of length ", Rc, sep=""))
    if (any(is.na(init.eps$sigma))) stop(paste("NA in init", number, ".eps$sigma", sep=""))        
    if (any(init.eps$sigma <= 0)) stop(paste("init", number, ".eps$sigma must be higher than ", 0, sep=""))
    names(init.eps$sigma) <- paste("sigma", 1:Rc, sep="")

    ##### init.eps:  gammaInv
    ##### -----------------------------------------------
    if (is.na(ieps.gammaInv)) init.eps$gammaInv <- prior.eps$zeta * isigma^2
    if (length(init.eps$gammaInv) == 1) init.eps$gammaInv <- rep(init.eps$gammaInv, Rc)
    if (length(init.eps$gammaInv) != Rc) stop(paste("init", number, ".eps$gammaInv must be of length ", Rc, sep=""))
    if (any(is.na(init.eps$gammaInv))) stop(paste("NA in init", number, ".eps$gammaInv", sep=""))        
    if (any(init.eps$gammaInv <= 0)) stop(paste("init", number, ".eps$gammaInv must be higher than ", 0, sep=""))
    names(init.eps$gammaInv) <- paste("gammaInv", 1:Rc, sep="")        
  }else{
    init.eps <- list(sigma=0, gammaInv=0)
  }  
  
  return(init.eps)
}
