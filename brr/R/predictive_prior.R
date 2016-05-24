#' @name Prior_y 
#' @rdname Prior_y
#' @title Prior predictive distribution of the count in the control group
#' @description Density, distribution function, quantile function and random 
#' generation for the prior predictive distribution of the count in the control group.
#' @details The prior predictive distribution of the count \eqn{y} is the Poisson-Gamma 
#'  distribution 
#' with shape parameter \eqn{a} and hyperrate parameter \eqn{b/T}, which is also 
#' the negative binomial distribution with shape \eqn{a} and probability of success 
#' \eqn{b/(b+T)}.
#' 
#' @param y,q vector of non-negative \strong{integer} quantiles 
#' @param p vector of probabilities
#' @param a,b non-negative shape parameter and rate parameter of the Gamma prior distribution on the rate \eqn{\mu}
#' @param T sample size of the control group
#' @param n number of observations to be simulated
#' @param ... other arguments passed to \code{\link[stats]{NegBinomial}} 
#' or \code{\link{summary_nbinom}} (for \code{sprior_y})
#' 
#' @return \code{dprior_y} gives the density, \code{pprior_y} the distribution function, 
#' \code{qprior_y} the quantile function, \code{rprior_y} samples from the distribution, 
#' and \code{sprior_y} gives a summary of the distribution.
#' 
#' @note \code{Prior_y} is a generic name for the functions documented. 
#' 
#' @examples 
#' barplot(dprior_y(0:10, 2, 2, 1))
#' sprior_y(2, 2, 1, output="pandoc")
#' @importFrom stats setNames dnbinom pnbinom qnbinom rnbinom
NULL
#'
#' @rdname Prior_y
#' @export 
dprior_y<-function(y, a, b, T, ...){
  return( setNames(dnbinom(y, a, b/(b+T), ...), paste("y=",y,sep="")) )
}
#'
#' @rdname Prior_y
#' @export 
pprior_y<-function(q, a, b, T, ...){
  return( setNames(pnbinom(q, a, b/(b+T), ...), paste("y\u2264",q,sep="")) )
}
#'
#' @rdname Prior_y
#' @export 
qprior_y<-function(p, a, b, T, ...){
  return( qnbinom(p, a, b/(b+T), ...) )
}
#'
#' @rdname Prior_y
#' @export 
rprior_y <- function(n, a, b, T){
  return( rnbinom(n, a, b/(b+T)) )
}
#'
#' @rdname Prior_y
#' @export 
sprior_y <- function(a, b, T, ...){
  return( summary_nbinom(a, b/(b+T), ...) )
}


#' @name Prior_x 
#' @rdname Prior_x
#' @title Prior predictive distribution of the count in the treated group
#' @description Density, distribution function, quantile function and random 
#' generation for the prior predictive distribution of the count in the treated group.
#' @details The prior predictive distribution of the count \eqn{x} is the 
#' \code{\link[=PGB2Dist]{Poisson-Gamma-Beta2 distribution}} 
#' with shape parameters \eqn{a}, \eqn{d}, \eqn{c},  
#' and hyperrate parameter \eqn{b/(b+T)}. 
#' 
#' @param x,q vector of non-negative \strong{integer} quantiles 
#' @param p vector of probabilities
#' @param a,b non-negative shape parameter and rate parameter of the Gamma prior distribution on the rate \eqn{\mu}
#' @param c,d non-negative shape parameters of the prior distribution on \eqn{\phi} 
#' @param T sample size of the control group
#' @param n number of observations to be simulated
#' @param ... passed to \code{\link{summary_PGB2}}
#' 
#' @return \code{dprior_x} gives the density, \code{pprior_x} the distribution function, 
#' \code{qprior_x} the quantile function, \code{rprior_x} samples from the distribution, 
#' and \code{sprior_x} gives a summary of the distribution.
#' 
#' @note \code{Prior_x} is a generic name for the functions documented. 
#' 
#' @examples 
#' barplot(dprior_x(0:30, 2, 3, 4, 5, 10))
#' sprior_x(2, 3, 4, 5, 10, output="pandoc")
#' @importFrom stats setNames
NULL
#'
#' @rdname Prior_x
#' @export 
dprior_x <- function(x, a, b, c, d, T){
    return( setNames(dPGB2(x,a,d,c,b/(b+T)), paste("x=",x,sep="")) )
}
#'
#' @rdname Prior_x
#' @export 
pprior_x <- function(q, a, b, c, d, T){
  return( setNames(pPGB2(q,a,d,c,b/(b+T)), paste("x\u2264",q,sep="")) )
}
#'
#' @rdname Prior_x
#' @export 
qprior_x <- function(p, a, b, c, d, T){
  return( qPGB2(p,a,d,c,b/(b+T)) ) 
}
#'
#' @rdname Prior_x
#' @export 
rprior_x <- function(n, a, b, c, d, T){
  return( rPGB2(n,a,d,c,b/(b+T)) )
}
#'
#' @rdname Prior_x
#' @export 
sprior_x <- function(a, b, c, d, T, ...){
  return( summary_PGB2(a,d,c,b/(b+T), ...) )
}


#' @name Prior_x_given_y 
#' @rdname Prior_x_given_y
#' @title Prior predictive distribution of the count \eqn{x} in the treated group 
#' conditionally to the count \eqn{y} in the treated group
#' @description Density, distribution function, quantile function and random 
#' generation for the conditional prior predictive distribution of \eqn{x} 
#' given \eqn{y}. 
#' @details The prior predictive distribution of the count \eqn{x} is the 
#' \code{\link[=BNBDist]{Beta-negative binomial distribution}} 
#' with shape parameters \eqn{a+y}, \eqn{d}, \eqn{c}.
#' 
#' @param x,q vector of non-negative \strong{integer} quantiles 
#' @param y count (integer) in the control group
#' @param p vector of probabilities
#' @param a non-negative shape parameter of the Gamma prior distribution on the rate \eqn{\mu}
#' @param c,d non-negative shape parameters of the prior distribution on \eqn{\phi} 
#' @param n number of observations to be simulated
#' @param ... arguments passed to \code{\link{summary_beta_nbinom}}
#' 
#' @return \code{dprior_x_given_y} gives the density, 
#' \code{pprior_x_given_y} the distribution function, 
#' \code{qprior_x_given_y} the quantile function, 
#' \code{rprior_x_given_y} samples from the distribution, and
#' \code{sprior_x_given_y} gives a summary of the distribution.
#' 
#' @note \code{Prior_x_given_y} is a generic name for the functions documented. 
#' 
#' @examples 
#' barplot(dprior_x_given_y(0:10, 5, 3, 10, 20))
#' sprior_x_given_y(5, 3, 10, 20, output="pandoc") 
#' @importFrom stats setNames
NULL
#'
#' @rdname Prior_x_given_y
#' @export 
dprior_x_given_y <- function(x, y, a, c, d){
  if(length(y)==1){
    return( setNames(dbeta_nbinom(x,a+y,d,c), paste("x=",x,sep="")) )
  }else{
    return( t(vapply(setNames(y, paste("y=",y,sep="")), function(y) dprior_x_given_y(x, y, a, c, d), numeric(length(x)))) )
  }
}
#'
#' @rdname Prior_x_given_y
#' @export 
pprior_x_given_y <- function(q, y, a, c, d){
  if(length(y)==1){
    return( setNames(pbeta_nbinom(q,a+y,d,c),  paste("x\u2264",q,sep="")) ) 
  }else{
    return( t(vapply(setNames(y, paste("y=",y,sep="")), function(y) pprior_x_given_y(q, y, a, c, d), numeric(length(q)))) )
  }
}
#'
#' @rdname Prior_x_given_y
#' @export 
qprior_x_given_y <- function(p, y, a, c, d){
  if(length(y)==1){
    return( qbeta_nbinom(p,a+y,d,c) )
  }else{
    return( t(vapply(setNames(y, paste("y=",y,sep="")), function(y) qprior_x_given_y(p, y, a, c, d), numeric(length(p)))) )
  }
}
#'
#' @rdname Prior_x_given_y
#' @export 
rprior_x_given_y <- function(n, y, a, c, d){
  return( rbeta_nbinom(n,a+y,d,c) )
}
#'
#' @rdname Prior_x_given_y
#' @export 
sprior_x_given_y <- function(y, a, c, d, ...){
  return( summary_beta_nbinom(a+y,d,c,...) )
}
