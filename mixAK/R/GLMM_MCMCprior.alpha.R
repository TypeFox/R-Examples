##
##  PURPOSE:   Handle prior.alpha argument of GLMM_MCMC function
##
##             THIS IS A HELP FUNCTION, NOT TO BE CALLED BY ORDINARY USERS
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    05/08/2009    as GLMM_MCMCprior.beta.R
##              20/01/2011    renamed to GLMM_MCMCprior.alpha
##
##  FUNCTIONS:  GLMM_MCMCprior.beta (originally)
##              GLMM_MCMCprior.alpha (from 20/01/2011)
##
## ================================================================================================

## *************************************************************
## GLMM_MCMCprior.alpha
## *************************************************************
##
GLMM_MCMCprior.alpha <- function(prior.alpha, lalpha)
{
##### Variables in the resulting object:
#####          prior.alpha  
#####          CpriorDouble_alpha
##### -----------------------------------------------------------------------------------------------------------    
  if (lalpha){
    if (missing(prior.alpha)) prior.alpha <- list()
    if (!is.list(prior.alpha)) stop("prior.alpha must be a list")

    inprior.alpha <- names(prior.alpha)
    ialpha.mean   <- match("mean", inprior.alpha, nomatch=NA)
    ialpha.var    <- match("var", inprior.alpha, nomatch=NA)    

    ##### prior.alpha:  mean
    ##### -----------------------------------------------
    if (is.na(ialpha.mean)) prior.alpha$mean <- rep(0, lalpha)
    if (length(prior.alpha$mean) == 1) prior.alpha$mean <- rep(prior.alpha$mean, lalpha)
    if (length(prior.alpha$mean) != lalpha) stop(paste("prior.alpha$mean must be of length", lalpha))
    if (any(is.na(prior.alpha$mean))) stop("NA in prior.alpha$mean")        
    Calphamean <- as.numeric(prior.alpha$mean)
    names(Calphamean) <- names(prior.alpha$mean) <- paste("alpha", 1:lalpha, ".mean", sep="")
    
    ##### prior.alpha:  var
    ##### -----------------------------------------------
    if (is.na(ialpha.var)) prior.alpha$var <- rep(10000, lalpha)
    if (length(prior.alpha$var) == 1) prior.alpha$var <- rep(prior.alpha$var, lalpha)
    if (length(prior.alpha$var) != lalpha) stop(paste("prior.alpha$var must be of length", lalpha))
    if (any(is.na(prior.alpha$var))) stop("NA in prior.alpha$var")
    if (any(prior.alpha$var <= 0)) stop(paste("prior.alpha$var must be higher than ", 0, sep=""))    
    Calphavar <- as.numeric(prior.alpha$var)
    names(Calphavar) <- names(prior.alpha$var) <- paste("alpha", 1:lalpha, ".var", sep="")

    ##### put all together
    ##### -----------------------------------------------
    CpriorDouble_alpha <- c(Calphamean, Calphavar)        
  }else{
    prior.alpha <- list(mean=0, var=0)
    CpriorDouble_alpha <- rep(0, 2)
  }  

  RET <- list(prior.alpha        = prior.alpha,
              CpriorDouble_alpha = CpriorDouble_alpha)
  return(RET)  
}
