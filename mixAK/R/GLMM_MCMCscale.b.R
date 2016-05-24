##
##  PURPOSE:   Handle scale.b argument of GLMM_MCMC function
##
##             THIS IS A HELP FUNCTION, NOT TO BE CALLED BY ORDINARY USERS
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    05/08/2009
##
##  FUNCTIONS:  GLMM_MCMCscale.b
##
## ================================================================================================

## *************************************************************
## GLMM_MCMCscale.b
## *************************************************************
##
GLMM_MCMCscale.b <- function(scale.b, dimb, iEranefVec, iSDranefVec)
{
##### Variables in the resulting object:
#####          scale.b
#####          CshiftScale_b   vector with shift and scale for random effects to be passed to C++  
##### -----------------------------------------------------------------------------------------------------------      
  if (dimb){
    if (missing(iEranefVec))  iEranefVec  <- rep(0, dimb)
    if (missing(iSDranefVec)) iSDranefVec <- rep(1, dimb)          
    if (missing(scale.b))     scale.b <- list(shift=iEranefVec, scale=iSDranefVec)
    if (!is.list(scale.b)) stop("scale.b must be a list")
    if (length(scale.b) != 2) stop("scale.b must have 2 components")
    inscale.b <- names(scale.b)
    iscale.b.shift <- match("shift", inscale.b, nomatch=NA)
    iscale.b.scale <- match("scale", inscale.b, nomatch=NA)
    if (is.na(iscale.b.shift)) stop("scale.b$shift is missing")
    if (length(scale.b$shift) == 1)    scale.b$shift <- rep(scale.b$shift, dimb)
    if (length(scale.b$shift) != dimb) stop(paste("scale.b$shift must be a vector of length ", dimb, sep=""))    
    if (is.na(iscale.b.scale)) stop("scale.b$scale is missing")
    if (length(scale.b$scale) == 1) scale.b$scale <- rep(scale.b$scale, dimb)
    if (length(scale.b$scale) != dimb) stop(paste("scale.b$scale must be a vector of length ", dimb, sep=""))
    if (any(scale.b$scale <= 0)) stop("all elements of scale.b$scale must be positive")
  }else{
    scale.b <- list(shift=0, scale=1)
  }  

  return(scale.b)  
}
