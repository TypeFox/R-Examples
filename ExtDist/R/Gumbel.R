#' @title The Gumbel distribution
#' @description Density, distribution, quantile, random number 
#' generation, and parameter estimation functions for the Gumbel distribution with parameters
#' \code{location} and \code{scale}. 
#' Parameter estimation can be based on a weighted or unweighted i.i.d sample and can be performed
#' analytically or numerically.
#' @rdname Gumbel
#' @name Gumbel
#'
#' @aliases dGumbel 
#' @aliases pGumbel
#' @aliases qGumbel
#' @aliases rGumbel
#' @aliases eGumbel
#' @aliases lGumbel
#' @details The \code{dGumbel()}, \code{pGumbel()}, \code{qGumbel()},and \code{rGumbel()} functions serve as wrappers of the 
#' \code{\link[VGAM]{dgumbel}}, \code{\link[VGAM]{pgumbel}}, \code{\link[VGAM]{qgumbel}}, and \code{\link[VGAM]{rgumbel}} functions 
#' in the \pkg{{VGAM}} package.They allow for the parameters to be declared not only as 
#' individual numerical values, but also as a list so parameter estimation can be carried out. \cr
#' \cr
#' The Gumbel distribution is a special case of the generalised extreme value (GEV) distribution and
#' has probability density function,
#' \deqn{f(x) = exp{(-exp{-(x-\mu)/\sigma)}}}
#' where \eqn{\mu} = \code{location} and \eqn{\sigma} = \code{scale} which has the constraint \eqn{\sigma > 0}.
#' The analytical parameter estimations are as given by the \href{http://www.itl.nist.gov/div898/handbook/eda/section3/eda366g.htm}{Engineering Statistics Handbook} 
#' with corresponding standard errors given by Bury (p.273).\cr
#' \cr  
#' The log-likelihood function of the Gumbel distribution is given by 
#' \deqn{l(\mu, \sigma| x) = \sigma^{-n} exp(-\sum (x_{i}-\mu/\sigma) - \sum exp(-(x_{i}-\mu/\sigma))).} 
#' Shi (1995) provides the score function and Fishers information matrix.  
#' @param params A list that includes all named parameters
#' @param x,q A vector of quantiles.
#' @param w An optional vector of sample weights.
#' @param p A vector of probabilities.
#' @param n Number of observations.
#' @param X Sample observations.
#' @param location Location parameter.
#' @param scale Scale parameter.
#' @param method Parameter estimation method.
#' @param logL logical if TRUE, lGumbel gives the log-likelihood, otherwise the likelihood is given.
#' @param ... Additional parameters.

#' @return dGumbel gives the density, pGumbel the distribution function,
#' qGumbel the quantile function, rGumbel generates random deviates, and 
#' eGumbel estimate the distribution parameters. lGumbel provides the log-likelihood function.
#' @author Haizhen Wu and A. Jonathan R. Godfrey. \cr Updates and bug fixes by Sarah Pirikahu. 
#' @references Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions,
#'  volume 2, chapter 22, Wiley, New York.\cr
#'  \cr
#'  \href{http://www.itl.nist.gov/div898/handbook/eda/section3/eda366g.htm}{Engineering Statistics Handbook}.\cr
#'  \cr
#'  Bury, K. (1999) Statistical Distributions in Engineering, Chapter 15, pp.283-284, 
#'  Cambridge University Press.\cr
#'  \cr
#'  Shi, D. (1995). Multivariate extreme value distribution and its Fisher information matrix. Acta Mathematicae
#   Applicatae Sinica
#' @seealso \pkg{\link{ExtDist}} for other standard distributions.
#' @examples
#' # Parameter estimation for a distribution with known shape parameters
#' X <- rGumbel(n = 500, location = 1.5, scale = 0.5)
#' est.par <- eGumbel(X, method="moments"); est.par
#' plot(est.par)
#' 
#' # Extracting location and scale parameters
#' est.par[attributes(est.par)$par.type=="location"]
#' est.par[attributes(est.par)$par.type=="scale"]
#'
#' #  Fitted density curve and histogram
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dGumbel(den.x, location = est.par$location, scale= est.par$scale)
#' hist(X, breaks=10, probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue")
#' lines(density(X))
#' 
#' # Parameter Estimation for a distribution with unknown shape parameters
#' # Example from; Bury(1999) pp.283-284, parameter estimates as given by Bury are location = 33.5 
#' # and scale = 2.241 
#' data <- c(32.7, 30.4, 31.8, 33.2, 33.8, 35.3, 34.6, 33, 32, 35.7, 35.5, 36.8, 40.8, 38.7, 36.7)
#' est.par <- eGumbel(X=data, method="numerical.MLE"); est.par
#' plot(est.par)
#' 
#' # log-likelihood
#' lGumbel(data, param = est.par)
#'
#' # Evaluating the precision of the parameter estimates by the Hessian matrix
#' H <- attributes(est.par)$nll.hessian
#' var <- solve(H)
#' se <- sqrt(diag(var)); se

#' @rdname Gumbel
#' @export dGumbel

dGumbel <-
  function(x, location = 0, scale = 1, params = list(location = 0, scale = 1),...){
    if(!missing(params)){
      location <- params$location
      scale <- params$scale
    }
    out = VGAM::dgumbel(x, location, scale)
    return(out)
  }

#' @rdname Gumbel
#' @export pGumbel

pGumbel <- 
  function(q, location = 0, scale = 1, params = list(location = 0, scale = 1),...){
    if(!missing(params)){
      location <- params$location
      scale <- params$scale
    }
    out = VGAM::pgumbel(q,location,scale)
    return(out)
}

#' @rdname Gumbel
#' @export qGumbel

qGumbel <- 
  function(p, location = 0, scale = 1, params = list(location = 0, scale = 1),...){
    if(!missing(params)){
      location <- params$location
      scale <- params$scale
    }
    out = VGAM::qgumbel(p,location,scale)
    return(out)
}

#' @rdname Gumbel
#' @export rGumbel
rGumbel <- 
  function(n, location = 0, scale = 1, params = list(location = 0, scale = 1),...){
    if(!missing(params)){
      location <- params$location
      scale <- params$scale
    }
    out = VGAM::rgumbel(n,location,scale)
    return(out)
  }

#' @rdname Gumbel
#' @export eGumbel
eGumbel <- function(X,w, method =c("moments","numerical.MLE"),...){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    # Analytical method of moments parameter estimation
  if(method=="moments"){ 
    sd <- sd(X*w)       # Weighted sd
    mean <- mean(X*w)   # Weighted mean
    n <- length(X)
    
    # Parameter estimates (http://www.itl.nist.gov/div898/handbook/eda/section3/eda366g.htm)
    scale <- (sd*sqrt(6))/pi
    location <- mean-0.5772*scale
    
    # Analtical estimates of standard errors (Bury(1999) p.273)
    bn <- 1 + (2.2/(n^1.13))
    SE.location <- 0.77970*((sd*bn)/sqrt(n))
    SE.scale <- 0.50697*sd*sqrt(bn/n)
    
    est.par <- list(location = location, scale = scale)
    est.par.se <- c(SE.location, SE.scale)
  }  
	else{(method == "numerical.MLE")    
      est.par <- wmle(X=X, w=w, distname = "Gumbel",
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
  attributes(est.par)$distname <- "Gumbel"
  attributes(est.par)$method <- method
  attributes(est.par)$par.name <- c("location","scale")
  attributes(est.par)$par.type <- c("location","scale")
  attributes(est.par)$par.vals <- c(est.par$location, est.par$scale)
  attributes(est.par)$par.s.e <-  est.par.se
    
    class(est.par) <- "eDist"
    
    return(est.par)
  }

#' @rdname Gumbel
#' @export lGumbel
## (weighted) (log) likelihood function
lGumbel <- 
  function(X, w, location = 0, scale = 1, params = list(location = 0, scale = 1), logL = TRUE,...){
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
    
    ll <- sum(w*log(dGumbel(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }
