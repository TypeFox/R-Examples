#' Numeric Delta Method approximation for the variance-covariance matrix
#'
#' Computes delta method variance-covariance matrix of results of any generic function \code{fct} that computes a vector of estimates as a function of a set of estimated parameters \code{par}.
#'
#' The delta method (aka propogation of errors is based on Taylor series approximation - see Seber's book on Estimation of Animal Abundance). It uses the first derivative of \code{fct} with respect to \code{par} which is computed in this function numerically using the central-difference formula. It also uses the variance-covariance matrix of the estimated parameters which is derived in estimating the parameters and is an input argument.
#'
#' The first argument of \code{fct} should be \code{par} which is a vector of parameter estimates. It should return a single value (or vector) of estimate(s).  The remaining arguments of \code{fct} if any can be passed to \code{fct} by including them at the end of the call to \code{DeltaMethod} as \code{name=value} pairs.
#'
#' @param par vector of parameter values at which estimates should be constructed
#' @param fct function that constructs estimates from parameters \code{par}
#' @param vcov variance-covariance matrix of the parameters
#' @param delta proportional change in parameters used to numerically estimate first derivative with central-difference formula
#' @param \dots any additional arguments needed by \code{fct}
#' @export
#' @return a list with values \item{variance}{estimated variance-covariance matrix of estimates derived by \code{fct}} \item{partial}{ matrix (or vector) of partial derivatives of \code{fct} with respect to the parameters \code{par}}
#' @note This is a generic function that can be used in any setting beyond the \code{mrds} package. However this is an internal function for \code{mrds} and the user does not need to call it explicitly.
#' @author Jeff Laake
#' @keywords utility
DeltaMethod <- function(par, fct, vcov, delta, ...){

  # Construct theta call to fct
  theta <- function(par) fct(par,...)

  # Numerically compute the first derivative of the estimator with
  # respect to each parameter using a central difference formula
  savepar <- par
  value1 <- theta(par)
  partial <- matrix(0,nrow=length(par),ncol=length(value1))
  for(i in 1:length(par)){
    # Store the original parameters into par and then adjust the
    # ith parameter by adding the proportion based on delta
    if(savepar[i]!=0){
      deltap <- delta*savepar[i]
    }else{
      deltap <- delta
    }

    par <- savepar
    par[i] <- savepar[i]+deltap

    # With this new value call the function to compute the estimate
    value1 <- theta(par)

    # Next do the same thing as above except substract off the proportion
    # delta from the original parameter.

    par <- savepar
    par[i] <- savepar[i]-deltap
    value2 <- theta(par)

    # Compute the central difference formula for the first derivative
    # with respect to the ith parameter

    partial[i,] <- (value1-value2)/(2*deltap)
  }
  variance <- t(partial)%*%vcov%*%partial

  # return the v-c matrix and the first partial vector(matrix)
  return(list(variance = variance,
              partial  = partial))
}
