##
##  PURPOSE:   Handle prior.eps argument of GLMM_MCMC function
##
##             THIS IS A HELP FUNCTION, NOT TO BE CALLED BY ORDINARY USERS
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    05/08/2009
##
##  FUNCTIONS:  GLMM_MCMCprior.eps
##
## ================================================================================================

## *************************************************************
## GLMM_MCMCprior.eps
## *************************************************************
##
GLMM_MCMCprior.eps <- function(prior.eps, Rc, isigma, is.sigma)
{
##### Variables in the resulting object:
#####          prior.eps
#####          CpriorDouble_eps
##### -----------------------------------------------------------------------------------------------------------    
  if (Rc){
    if (missing(isigma)){
      isigma <- rep(1, Rc)
      is.sigma <- rep(TRUE, Rc)
    }else{
      if (missing(is.sigma)) stop("argument is.sigma may not be missing when isigma is given")
    }  
    Reps <- 6*isigma[is.sigma]
    
    if (missing(prior.eps)) prior.eps <- list()
    if (!is.list(prior.eps)) stop("prior.eps must be a list")    
    
    inprior.eps <- names(prior.eps)
    ieps.zeta   <- match("zeta", inprior.eps, nomatch=NA)
    ieps.g      <- match("g", inprior.eps, nomatch=NA)
    ieps.h      <- match("h", inprior.eps, nomatch=NA)    

    ##### prior.eps:  zeta
    ##### -----------------------------------------------
    if (is.na(ieps.zeta)) prior.eps$zeta <- rep(2, Rc)
    if (length(prior.eps$zeta) == 1) prior.eps$zeta <- rep(prior.eps$zeta, Rc)
    if (length(prior.eps$zeta) != Rc) stop(paste("prior.eps$zeta must be of length", Rc))
    if (any(is.na(prior.eps$zeta))) stop("NA in prior.eps$zeta")        
    if (any(prior.eps$zeta <= 0)) stop(paste("prior.eps$zeta must be higher than ", 0, sep=""))
    Cepszeta <- as.numeric(prior.eps$zeta)
    names(Cepszeta) <- names(prior.eps$zeta) <- paste("zeta", 1:Rc, sep="")

    ##### prior.eps:  g
    ##### -----------------------------------------------
    if (is.na(ieps.g)) prior.eps$g <- rep(0.2, Rc)
    if (length(prior.eps$g) == 1) prior.eps$g <- rep(prior.eps$g, Rc)
    if (length(prior.eps$g) != Rc) stop(paste("prior.eps$g must be of length", Rc))
    if (any(is.na(prior.eps$g))) stop("NA in prior.eps$g")    
    if (any(prior.eps$g <= 0)) stop(paste("prior.eps$g must be higher than ", 0, sep=""))
    Cepsg <- as.numeric(prior.eps$g)
    names(Cepsg) <- names(prior.eps$g) <- paste("g", 1:Rc, sep="")
    
    ##### prior.eps:  h
    ##### -----------------------------------------------
    if (is.na(ieps.h)) prior.eps$h <- 10/(Reps^2)
    if (length(prior.eps$h) == 1) prior.eps$h <- rep(prior.eps$h, Rc)
    if (length(prior.eps$h) != Rc) stop(paste("prior.eps$h must be of length", Rc))
    if (any(is.na(prior.eps$h))) stop("NA in prior.eps$h")    
    if (any(prior.eps$h <= 0)) stop(paste("prior.eps$h must be higher than ", 0, sep=""))
    Cepsh <- as.numeric(prior.eps$h)
    names(Cepsh) <- names(prior.eps$h) <- paste("h", 1:Rc, sep="")

    ##### put all together
    ##### -----------------------------------------------
    CpriorDouble_eps <- c(Cepszeta, Cepsg, Cepsh)    
  }else{
    prior.eps <- list(zeta=0, g=0, h=0)
    CpriorDouble_eps <- rep(0, 3)
  }  

  RET <- list(prior.eps        = prior.eps,
              CpriorDouble_eps = CpriorDouble_eps)
  return(RET)      
}  
