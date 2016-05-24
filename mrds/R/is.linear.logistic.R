#' Collection of functions for logistic detection functions
#'
#' These functions are used to test whether a logistic detection function is a
#' linear function of distance (\code{is.linear.logistic}) or is constant
#' (varies by distance but no other covariates) \code{is.logistic.constant}).
#' Based on these tests, the most appropriate manner for integrating the
#' detection function with respect to distance is chosen.  The integrals are
#' needed to estimate the average detection probability for a given set of
#' covariates.
#'
#' If the logit is linear in distance then the integral can be computed
#' analytically. If the logit is constant or only varies by distance then only
#' one integral needs to be computed rather than an integral for each
#' observation.
#'
#' @param xmat data matrix
#' @param g0model logit model
#' @param zdim number of columns in design matrix
#' @param width transect width
#' @return Logical TRUE if condition holds and FALSE otherwise
#' @author Jeff Laake
#' @keywords utility
is.linear.logistic <- function(xmat,g0model,zdim,width){

  xmat$distance <- rep(width/2, nrow(xmat))
  beta <- rep(1,zdim)
  logit1 <- mean(beta %*% t(setcov(xmat, g0model)))
  xmat$distance <- rep(width, nrow(xmat))
  logit2 <- mean(beta %*% t(setcov(xmat, g0model)))
  xmat$distance <- rep(0, nrow(xmat))
  logit0 <- mean( beta %*% t(setcov(xmat, g0model)))

  if(is.nan(logit1) || is.nan(logit0)){
    integral.numeric <- TRUE
  }else if(logit1-logit0==0){
    integral.numeric <- FALSE
  }else if((logit2-logit0)/(logit1-logit0) <= 2.00001 &
          (logit2-logit0)/(logit1-logit0) >= 1.99999){
    integral.numeric <- FALSE
  }else{
    integral.numeric <- TRUE
  }

  return(integral.numeric)
}
