#' @title The Logistic Distribution.
#' @description Density, distribution, and quantile, random number generation, 
#' and parameter estimation functions for the logistic distribution with parameters \code{location} and \code{scale}. 
#' Parameter estimation can be based on a weighted or unweighted i.i.d. sample and can be carried out numerically.
#' @rdname Logistic
#' @name Logistic
#'
#' @aliases dLogistic
#' @aliases pLogistic
#' @aliases qLogistic
#' @aliases rLogistic
#' @aliases eLogistic
#' @aliases lLogistic
#' @details If \code{location} or \code{scale} are omitted, they assume the default values of 0 or 1 
#' respectively. \cr
#' \cr
#' The \code{dLogistic()}, \code{pLogistic()}, \code{qLogistic()},and \code{rLogistic()} functions serve as wrappers of the
#' standard \code{\link[stats]{dlogis}}, \code{\link[stats]{plogis}}, \code{\link[stats]{qlogis}}, and 
#' \code{\link[stats]{rlogis}} functions in the \pkg{\link{stats}} package. They allow for the parameters to be declared not only as 
#' individual numerical values, but also as a list so parameter estimation can be carried out. \cr  
#' \cr
#' The logistic distribution with \code{location} = \eqn{\alpha} and \code{scale} = \eqn{\beta} is most simply 
#' defined in terms of its cumulative distribution function (Johnson et.al pp.115-116)  
#' \deqn{F(x) =  1- [1 + exp((x-\alpha)/\beta)]^{-1}.}
#' The corresponding probability density function is given by
#' \deqn{f(x) = 1/\beta [exp(x-\alpha/\beta][1 + exp(x-\alpha/\beta)]^{-2}}
#' Parameter estimation is only implemented numerically. \cr
#' \cr
#' The score function and Fishers information are as given by Shi (1995) (See also Kotz & Nadarajah (2000)).
#' @param params A list that includes all named parameters.
#' @param x,q A vector of quantiles.
#' @param w An optional vector of sample weights.
#' @param p A vector of probabilities.
#' @param n Number of observations.
#' @param X Sample observations.
#' @param location Location parameter.
#' @param scale Scale parameter.
#' @param method Parameter estimation method.
#' @param logL logical; if TRUE, lLogistic gives the log-likelihood, otherwise the likelihood is given.
#' @param ... Additional parameters.
#'
#' @return dLogistic gives the density, pLogistic the distribution function,
#' qLogistic the quantile function, rLogistic generates random deviates, and  
#' eLogistic estimates the parameters. lLogistic provides the log-likelihood function.
#' @seealso \pkg{\link{ExtDist}} for other standard distributions.
#' @references Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions, volume 2,
#'  chapter 23. Wiley, New York.\cr
#'  \cr
#'  Shi, D. (1995) Fisher information for a multivariate extreme value distribution, Biometrika, vol 82, pp.644-649.\cr
#'  \cr
#'  Kotz, S. and Nadarajah (2000) Extreme Value Distributions Theory and Applications, chapter 3, Imperial Collage Press,
#'  Singapore.
#' @author Haizhen Wu and A. Jonathan R. Godfrey.\cr
#' Updates and bug fixes by Sarah Pirikahu.
#'
#' @examples
#' # Parameter estimation for a distribution with known shape parameters
#' X <- rLogistic(n=500, location=1.5, scale=0.5)
#' est.par <- eLogistic(X); est.par
#' plot(est.par)
#' #  Fitted density curve and histogram
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dLogistic(den.x,location=est.par$location,scale=est.par$scale)
#' hist(X, breaks=10, probability=TRUE, ylim = c(0,1.2*max(den.y)))
#' lines(den.x, den.y, col="blue")
#' lines(density(X), lty=2)
#' 
#' # Extracting location or scale parameters
#' est.par[attributes(est.par)$par.type=="location"]
#' est.par[attributes(est.par)$par.type=="scale"]
#' 
#' # log-likelihood function
#' lLogistic(X,param = est.par)
#' 
#' # Evaluation of the precision of the parameter estimates by the Hessian matrix
#' H <- attributes(est.par)$nll.hessian
#' fisher_info <- solve(H)
#' var <- sqrt(diag(fisher_info));var
#' 
#' # Example of parameter estimation for a distribution with 
#' # unknown parameters currently been sought after.

#' @rdname Logistic
#' @export dLogistic

dLogistic <-function(x, location = 0, scale = 1, params = list(location = 0, scale = 1),...){
    if(!missing(params)){
      location <- params$location
      scale <- params$scale
    }
    out = dlogis(x, location, scale)
    return(out)
  }

#' @rdname Logistic
#' @export pLogistic

pLogistic <- function(q, location = 0, scale = 1, params = list(location = 0, scale = 1),...){
    if(!missing(params)){
      location <- params$location
      scale <- params$scale
    }
    out = plogis(q,location,scale)
    return(out)
}

#' @rdname Logistic
#' @export qLogistic

qLogistic <- function(p, location = 0, scale = 1, params = list(location = 0, scale = 1),...){
    if(!missing(params)){
      location <- params$location
      scale <- params$scale
    }
    out = qlogis(p,location,scale)
    return(out)
}

#' @rdname Logistic
#' @export rLogistic
rLogistic <- function(n, location = 0, scale = 1, params = list(location = 0, scale = 1),...){
    if(!missing(params)){
      location <- params$location
      scale <- params$scale
    }
    out = rlogis(n,location,scale)
    return(out)
  }

#' @rdname Logistic
#' @export eLogistic
eLogistic <- function(X,w, method ="numerical.MLE",...){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
	{if(method != "numerical.MLE") warning(paste("method ", method, " is not avaliable, use numerial.MLE instead."))  
      method = "numerical.MLE"  
      est.par <- wmle(X=X, w=w, distname = "Logistic",
                      initial=list(location = 0, scale = 1),
                      lower=list(location = -Inf, scale = 0),
                      upper=list(location = Inf, scale = Inf))

      est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
      if(class(est.par.se) == "try-error") {
        est.par.se <- rep(NA, length(est.par))
      } 
    }
    
    attributes(est.par)$ob <- X
    attributes(est.par)$weights <- w
    attributes(est.par)$distname <- "Logistic"
    attributes(est.par)$method <- method
    attributes(est.par)$par.name <- c("location","scale")
    attributes(est.par)$par.type <- c("location","scale")
    attributes(est.par)$par.vals <- c(est.par$location, est.par$scale)
    attributes(est.par)$par.s.e <-  est.par.se  
    
    class(est.par) <- "eDist"
    
    return(est.par)
  }

#' @rdname Logistic
#' @export lLogistic
## (weighted) (log) likelihood function
lLogistic <- function(X, w, location = 0, scale = 1, params = list(location = 0, scale = 1), logL = TRUE,...){
    if(!missing(params)){
      location <- params$location
      scale <- params$scale
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dLogistic(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }
