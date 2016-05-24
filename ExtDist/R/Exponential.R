#' @title The Exponential Distribution.
#' @description Density, distribution, quantile, random number 
#' generation and parameter estimation functions for the exponential distribution. 
#' Parameter estimation can be based on a weighted or unweighted i.i.d sample and is carried out 
#' analytically. 
#' @rdname Exponential 
#' @name Exponential
#'
#' @aliases dExp
#' @aliases pExp
#' @aliases qExp 
#' @aliases rExp
#' @aliases eExp
#' @aliases lExp
#' @aliases sExp
#' @aliases iExp
#'
#' @param params A list that includes all named parameters
#' @param x,q A vector of sample values or quantiles.
#' @param w An optional vector of sample weights.
#' @param p A vector of probabilities.
#' @param n Number of observations.
#' @param scale scale parameter, called rate in other packages.
#' @param method Parameter estimation method.
#' @param logL logical; if TRUE, lExp gives the log-likelihood, otherwise the likelihood is given.
#' @param ... Additional parameters. 
#'      
#' @return dExp gives the density, pExp the distribution function,
#' qExp the quantile function, rExp generates random deviates, and 
#' eExp estimates the distribution parameters. lExp provides the log-likelihood function.
#' 
#' @details If \code{scale} is omitted, it assumes the default value 1 giving the 
#' standard exponential distribution.\cr
#' \cr 
#' The exponential distribution is a special case of the gamma distribution where the shape parameter 
#' \eqn{\alpha = 1}. The \code{dExp()}, \code{pExp()},
#' \code{qExp()},and \code{rExp()} functions serve as wrappers of the standard \code{\link[stats]{dexp}}, 
#' \code{\link[stats]{pexp}}, \code{\link[stats]{qexp}} and \code{\link[stats]{rexp}} functions  
#' in the \pkg{\link{stats}} package. They allow for the parameters to be declared not only as 
#' individual numerical values, but also as a list so parameter estimation can be carried out. \cr 
#' \cr
#' The probability density function for the exponential distribution with \code{scale}=\eqn{\beta} is
#' \deqn{f(x) = (1/\beta) * exp(-x/\beta)} 
#' for \eqn{\beta > 0 }, Johnson et.al (Chapter 19, p.494). Parameter estimation for the exponential distribution is 
#' carried out analytically using maximum likelihood estimation (p.506 Johnson et.al).\cr
#' \cr
#' The likelihood function of the exponential distribution is given by 
#' \deqn{l(\lambda|x) = n log \lambda - \lambda \sum xi.}    
#' It follows that the score function is given by
#' \deqn{dl(\lambda|x)/d\lambda = n/\lambda - \sum xi}
#' and Fisher's information given by
#' \deqn{E[-d^2l(\lambda|x)/d\lambda^2] = n/\lambda^2.}
#' 
#' @references Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions,
#'  volume 1, chapter 19, Wiley, New York.\cr
#' \cr
#' Kapadia. A.S., Chan, W. and Moye, L. (2005) Mathematical Statistics with Applications, Chapter 8, 
#' Chapman& Hall/CRC.
#' 
#' @author Jonathan R. Godfrey and Sarah Pirikahu.
#' 
#' @examples
#' # Parameter estimation for a distribution with known shape parameters
#' x <- rExp(n=500, scale=2)
#' est.par <- eExp(x); est.par
#' plot(est.par)
#' 
#' #  Fitted density curve and histogram
#' den.x <- seq(min(x),max(x),length=100)
#' den.y <- dExp(den.x,scale=est.par$scale)
#' hist(x, breaks=10, probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue")
#' lines(density(x), lty=2)
#'   
#' # Extracting the scale parameter
#' est.par[attributes(est.par)$par.type=="scale"]
#' 
#' # Parameter estimation for a distribution with unknown shape parameters
#' # Example from Kapadia et.al(2005), pp.380-381. 
#' # Parameter estimate as given by Kapadia et.al is scale=0.00277
#' cardio <- c(525, 719, 2880, 150, 30, 251, 45, 858, 15, 
#'            47, 90, 56, 68, 6, 139, 180, 60, 60, 294, 747)
#' est.par <- eExp(cardio, method="analytical.MLE"); est.par
#' plot(est.par)
#' 
#' # log-likelihood, score function and Fisher's information
#' lExp(cardio,param = est.par)
#' sExp(cardio,param = est.par)
#' iExp(cardio,param = est.par)

#' @rdname Exponential
#' @export dExp
dExp <-function(x, scale = 1, params = list(scale = 1),...){
  if(!missing(params)){
    scale <- params$scale
  }
  out = stats::dexp(x, rate=scale)
  return(out)
}

#' @rdname Exponential
#' @export pExp
pExp <- function(q, scale = 1, params = list(scale = 1),...){
  if(!missing(params)){
    scale <- params$scale
  }
  out = stats::pexp(q,rate=scale)
  return(out)
}

#' @rdname Exponential
#' @export qExp
qExp <- function(p, scale = 1, params = list(scale = 1),...){
  if(!missing(params)){
    scale <- params$scale
  }
  out = stats::qexp(p, rate=scale)
  return(out)
}

#' @rdname Exponential
#' @export rExp
rExp <- function(n, scale = 1, params = list(scale = 1),...){
  if(!missing(params)){
    scale <- params$scale
  }
  out = stats::rexp(n, rate=scale)
  return(out)
}

#' @rdname Exponential
#' @export eExp
eExp<- function(x,w, method ="analytical.MLE",...){
  method <- match.arg(method)
  # Accounting for weights
  n <- length(x)
  if(missing(w)){
    w <- rep(1,n)
  } else {
    w <- n*w/sum(w)
  }
  method <- "analytical.MLE"
  # Maximum likelihood estimator p. 506 Johnson & Kotz.
  lambda <- 1/mean(x*w)
  SE.lambda <- lambda/sqrt(length(x))  
  
  est.par <- list(scale = lambda)
  est.par.se <- c(SE.lambda)
  
  attributes(est.par)$ob <- x
  attributes(est.par)$weights <- w
  attributes(est.par)$distname <- "Exp"
  attributes(est.par)$method <- method
  attributes(est.par)$par.name <- c("rate")
  attributes(est.par)$par.type <- c("scale")
  attributes(est.par)$par.vals <- c(est.par$scale)
  attributes(est.par)$par.s.e <-  est.par.se  
  
  class(est.par) <- "eDist"
  
  return(est.par)
}

#' @rdname Exponential
#' @export lExp
## (weighted) (log) likelihood function
lExp <- function(x, w, scale = 1, params = list(scale = 1), logL = TRUE,...){
    if(!missing(params)){
      scale <- params$scale
    }
    
    n <- length(x)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dExp(x=x,params = params)))
    # ll <- (shape-1)*sum(log(x)) - sum(x/scale) - n*shape*log(scale) + n*log(gamma(shape))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }

#' @rdname Exponential
#' @export sExp
sExp <- function(x, w, scale=1, params = list(scale=1),...){
  if(!missing(params)){
    scale <- params$scale
  }
  n <- length(x)
  if(missing(w)){
    w <- rep(1,n)
  } else {
    w <- n*w/sum(w)
  }
  score <- (n/scale) - sum(w*x)
  return(score)
}

#' @rdname Exponential
#' @export iExp
iExp <- function(x, w, scale = 1, params = list(scale=1), ...){
  if(!missing(params)){
    scale <- params$scale
  }
  n <- length(x)
  if(missing(w)){
    w <- rep(1,n)
  } else {
    w <- n*w/sum(w)
  }
  info <- n/scale^2
  return(info) 
}

  

