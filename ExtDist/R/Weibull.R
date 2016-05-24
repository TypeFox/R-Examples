#' @title The Weibull Distribution.
#' @description Density, distribution, quantile, random number 
#' generation, and parameter estimation functions for the Weibull distribution with parameters \code{shape} and 
#' \code{scale}. Parameter estimation can be based on a weighted or unweighted i.i.d sample and can be carried out analytically 
#' or numerically. 
#' @rdname Weibull
#' @name Weibull
#'
#' @aliases dWeibull 
#' @aliases pWeibull
#' @aliases qWeibull
#' @aliases rWeibull
#' @aliases eWeibull
#' @aliases lWeibull
#' @details The Weibull distribution is a special case of the generalised gamma distribution. The \code{dWeibull()}, \code{pWeibull()},
#' \code{qWeibull()},and \code{rWeibull()} functions serve as wrappers of the standard \code{\link[stats]{dgamma}}, 
#' \code{\link[stats]{pgamma}}, \code{\link[stats]{qgamma}}, and \code{\link[stats]{rgamma}} functions with  
#' in the \pkg{\link{stats}} package. They allow for the parameters to be declared not only as 
#' individual numerical values, but also as a list so parameter estimation can be carried out. \cr
#' \cr
#' The Weibull distribution with parameters \code{shape}=\eqn{a} and \code{scale}=\eqn{b} has probability density function,
#' \deqn{f(x)= (a/b)(x/b)^{a-1}exp(-(x/b)^a)}
#' for \eqn{x >0}. Parameter estimation can be carried out using the method of moments as done by Winston (2003) or by numerical
#' maximum likelihood estimation.\cr
#' \cr
#' The log-likelihood function of the Weibull distribution is given by
#' \deqn{l(a,b|x) = n(ln a - ln b) + (a-1)\sum ln(xi/b) - \sum(xi/b)^{a}}  
#' The score function and information matrix are as given by Rinne (p.412).
#' 
#' @param params A list that includes all named parameters
#' @param x,q A vector of quantiles.
#' @param w An optional vector of sample weights.
#' @param p A vector of probabilities.
#' @param n Number of observations.
#' @param X Sample observations.
#' @param shape Shape parameter.
#' @param scale Scale parameter.
#' @param method Parameter estimation method.
#' @param logL logical; if TRUE, lWeibull gives the log-likelihood, otherwise the likelihood is given.
#' @param ... Additional parameters. 
#'
#' @return dWeibull gives the density, pWeibull the distribution function,
#' qWeibull the quantile function, rWeibull generates random deviates, and 
#' eWeibull estimates the distribution parameters. lWeibull provides the log-likelihood function.
#' @seealso \pkg{\link{ExtDist}} for other standard distributions.
#' @references Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions,
#'  volume 1, chapter 21, Wiley, New York.\cr
#'  \cr
#'  Rinne, H. (2009) The Weibull Distribution A Handbook, chapter 11, Chapman & Hall/CRC.\cr
#'  \cr
#'  Winston, W.L (2003) Operations Research: Applications and algorithms, 4th Ed, Duxbury.
#' 
#' @author Haizhen Wu and A. Jonathan R. Godfrey. \cr
#' Updates and bug fixes by Sarah Pirikahu.
#'
#' @examples
#' # Parameter estimation for a distribution with known shape parameters
#' X <- rWeibull(n=1000, params=list(shape=1.5, scale=0.5))
#' est.par <- eWeibull(X=X, method="numerical.MLE"); est.par
#' plot(est.par)
#' 
#' #  Fitted density curve and histogram
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dWeibull(den.x,shape=est.par$shape,scale=est.par$scale)
#' hist(X, breaks=10, col="red", probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue", lwd=2)   # Original data
#' lines(density(X), lty=2)                 # Fitted curve
#' 
#' # Extracting shape and scale parameters
#' est.par[attributes(est.par)$par.type=="shape"]
#' est.par[attributes(est.par)$par.type=="scale"]
#'
#' # Parameter Estimation for a distribution with unknown shape parameters
#' # Example from: Rinne (2009) Dataset p.338 and example pp.418-419
#' # Parameter estimates are given as shape = 99.2079 and scale = 2.5957. The log-likelihood
#' # for this data and Rinne's parameter estimates is -1163.278.
#' data <- c(35,38,42,56,58,61,63,76,81,83,86,90,99,104,113,114,117,119,141,183)
#' est.par <- eWeibull(X=data, method="numerical.MLE"); est.par
#' plot(est.par)
#'
#' # Estimates calculated by eWeibull differ from those given by Rinne(2009).
#' # However, eWeibull's parameter estimates appear to be an improvement, due to a larger 
#' # log-likelihood of -99.09037 (as given by lWeibull below). 
#'
#'  # log-likelihood function
#' lWeibull(data, param = est.par)
#'
#' # evaluate the precision of estimation by Hessian matrix
#' H <- attributes(est.par)$nll.hessian
#' var <- solve(H)
#' se <- sqrt(diag(var));se

#' @rdname Weibull
#' @export dWeibull

dWeibull <-function(x, shape = 2, scale = 2, params = list(shape = 2, scale = 2)){
    if(!missing(params)){
      shape <- params$shape
      scale <- params$scale
    }
    out = stats::dgamma(x, shape, scale)
    return(out)
  }

#' @rdname Weibull
#' @export pWeibull

pWeibull <- function(q, shape = 2, scale = 2, params = list(shape = 2, scale = 2)){
    if(!missing(params)){
      shape <- params$shape
      scale <- params$scale
    }
    out = stats::pgamma(q,shape,scale)
    return(out)
  }

#' @rdname Weibull
#' @export qWeibull

qWeibull <- function(p, shape = 2, scale = 2, params = list(shape = 2, scale = 2)){
    if(!missing(params)){
      shape <- params$shape
      scale <- params$scale
    }
    out = stats::qgamma(p,shape,scale)
    return(out)
  }

#' @rdname Weibull
#' @export rWeibull
rWeibull <- function(n, shape = 2, scale = 2, params = list(shape = 2, scale = 2)){
    if(!missing(params)){
      shape <- params$shape
      scale <- params$scale
    }
    out = stats::rgamma(n,shape,scale)
    return(out)
  }

#' @rdname Weibull
#' @export eWeibull
eWeibull <- function(X,w, method =c("numerical.MLE","moments"),...){
  method <- match.arg(method)
  n <- length(X)
  if(missing(w)){
    w <- rep(1,n)
  } else {
    w <- n*w/sum(w)
  }
  
  # Method of moments parameter estimation
  if(method == "moments") {
    # Estimation of Weibull parameters based on first two moments only
    # built up from spreadsheet provided by Winston, W.L. Operations Research, 4th edition, Duxbury.
    
    WeibullParameters = function(x=NULL, Mean=mean(x), Var=var(x), Second=Var+Mean^2, Output=FALSE){
      # allows three options for data entry 
      # 1 - a vector, 
      # 2 - the mean and variance of the data
      # 3 - the mean and second moment of the data
      Ratio = Second/Mean^2
      Func = function (Estimate, Working){
        W1 =  exp(lgamma(1+(2/Estimate))) 
        W2 =  exp(lgamma(1+(1/Estimate)))^2
        return(Working - W1/W2)
      }
      Alpha = uniroot(Func, lower=1, upper=10, Working=Ratio)$root
      Beta = Mean/ exp(lgamma(1+(1/Alpha)))
      Param = c(Alpha, Beta)
      names(Param) = c("Alpha", "Beta")
      if(Output) cat(paste(names(Param),"parameter =",round(Param,4),"\n", sep=" "))
      return(Param)
    }
    # Parameter estimation via the WeibullParameters function
    params <- WeibullParameters(X) 
    
    # Extraction of shape and scale parameters
    shape <- as.numeric(params[1])
    scale <- as.numeric(params[2])
    
    est.par <- est.par <- list(shape = shape, scale = scale)
    est.par.se <- rep(NA, length(est.par))
  }     
  
{if( method == "numerical.MLE") 
  est.par <- wmle(X=X, w=w, distname = "Weibull",
                  initial=list(shape = 1, scale = 1),
                  lower=list(shape = 0, scale = 0),
                  upper=list(shape = Inf, scale = Inf))
 
 est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
 if(class(est.par.se) == "try-error") {
   est.par.se <- rep(NA, length(est.par))
 } 
}

attributes(est.par)$ob <- X
attributes(est.par)$weights <- w
attributes(est.par)$distname <- "Weibull"
attributes(est.par)$method <- method
attributes(est.par)$par.name <- c("shape","scale")
attributes(est.par)$par.type <- c("shape","scale")
attributes(est.par)$par.vals <- c(est.par$shape, est.par$scale)
attributes(est.par)$par.s.e <-  est.par.se  

class(est.par) <- "eDist"

return(est.par)
}

#' @rdname Weibull
#' @export lWeibull
## (weighted) (log) likelihood function
lWeibull <- 
  function(X, w, shape = 2, scale = 2, params = list(shape = 2, scale = 2), logL = TRUE){
    if(!missing(params)){
      shape <- params$shape
      scale <- params$scale
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dWeibull(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }
