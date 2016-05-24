#'@title Regularization criterion types
#'@name Methods
#'@export
#'@rdname method
#'@description
#'\itemize{
#'  \item{\code{LPLS} for PLS-type SCGLR}
#'  \item{\code{SR} Method iterative normed gradient (ING) for Structural Relevance}
#'}
methodLPLS <- function() {
  structure(list(
      method="lpls"
    ),
    class="method.SCGLR",
    description="PLS-type SCGLR"
    )
}

#'@rdname method
#'@export
#'@param phi character string describing structura relevance used in the regularization process.
#' Allowed values are "vpi" for Variable Powered Inertia and "cv" for Component Variance. Default to "vpi".
#'@param l is a numeric argument (>1) tuning the importance of variable bundle locality.
#'@param s is a numeric argument (in [0,1]) tuning the strength of structural relevance with respect to goodness of fit.
#'@param maxiter integer for maximum number of iterations of \code{SR} function
#'@param epsilon positive convergence threshold
#'@param bailout integer argument 
methodSR <- function(phi="vpi",l=1,s=1/2,maxiter=1000,epsilon=1e-6,bailout=10) {
  # check arguments
  if(!(phi %in% c("vpi","cv"))) 
    stop("phi should be \"vpi\" or \"cv\"")
  if(!is.numeric(l) || l<1) 
    stop("l must be greater than 1")
  if(!is.numeric(s) || s<0 || s>1) 
    stop("s must be between 0 and 1")
  if(!is.numeric(maxiter) || maxiter<1) 
    stop("maxiter must be an integer greater than 1")
  if(!is.numeric(epsilon) || epsilon<=0) 
    stop("epsilon must be a positive numeric")
  if(!is.numeric(bailout) || bailout<1) 
    stop("bailout must be an integer greater than 1")
  
  structure(list(
    method="sr",
    phi=phi,
    l=l,
    s=s,
    maxiter=maxiter,
    epsilon=epsilon,
    bailout=bailout#1000
  ),
  class="method.SCGLR",
  description="Method iterative normed gradient (ING) for Structural Relevance"
  )
}
