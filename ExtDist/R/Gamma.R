#' @title The Gamma Distribution.
#' @description Density, distribution, quantile, random number 
#' generation, and parameter estimation functions for the gamma distribution with parameters \code{shape} and 
#' \code{scale}. Parameter estimation can be based on a weighted or unweighted i.i.d sample and can be carried out 
#' numerically. 
#' @rdname Gamma
#' @name Gamma
#'
#' @aliases dGamma
#' @aliases pGamma 
#' @aliases qGamma 
#' @aliases rGamma
#' @aliases eGamma
#' @aliases lGamma
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
#' @param logL logical; if TRUE, lBeta_ab gives the log-likelihood, otherwise the likelihood is given.
#' @param ... Additional parameters. 
#'
#' @details The \code{dGamma()}, \code{pGamma()}, \code{qGamma()},and \code{rGamma()} functions serve as wrappers of the standard
#' \code{\link[stats]{dgamma}}, \code{\link[stats]{pgamma}}, \code{\link[stats]{qgamma}}, and \code{\link[stats]{rgamma}} functions 
#' in the \pkg{\link{stats}} package. They allow for the parameters to be declared not only as 
#' individual numerical values, but also as a list so parameter estimation can be carried out. \cr
#' \cr
#' The gamma distribution with parameter \code{shape}=\eqn{\alpha} and \code{scale}=\eqn{\beta} has probability density function,
#' \deqn{f(x)= (1/\beta^\alpha \Gamma(\alpha))x^{\alpha-1}e^{-x/\beta}}
#'  where \eqn{\alpha > 0} and \eqn{\beta > 0}. Parameter estimation can be performed using the method of moments
#'  as given by Johnson et.al (pp.356-357). \cr
#'  \cr
#'  The log-likelihood function of the gamma distribution is given by,
#'  \deqn{l(\alpha, \beta |x) = (\alpha -1) \sum_i ln(x_i) - \sum_i(x_i/\beta) -n\alpha ln(\beta) + n ln \Gamma(\alpha)}
#'  where \eqn{\Gamma} is the \link[base]{gamma} function. The score function is provided by Rice (2007), p.270. 
#' @return dGamma gives the density, pGamma the distribution function,
#' qGamma the quantile function, rGamma generates random deviates, and 
#' eGamma estimates the distribution parameters.lgamma provides the log-likelihood function.
#' 
#' @references Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions,
#'  volume 1, chapter 17, Wiley, New York.\cr
#' \cr
#'  Bury, K. (1999) Statistical Distributions in Engineering, Chapter 13, pp.225-226, 
#'  Cambridge University Press.\cr
#'  \cr
#'  Rice, J.A. (2007) Mathematical Statistics and Data Analysis, 3rd Ed, Brookes/Cole.
#' @seealso \pkg{\link{ExtDist}} for other standard distributions.
#' @author Haizhen Wu and A. Jonathan R. Godfrey. \cr
#' Updates and bug fixes by Sarah Pirikahu.
#' @examples
#' # Parameter estimation for a distribution with known shape parameters
#' X <- rGamma(n=500, shape=1.5, scale=0.5)
#' est.par <- eGamma(X); est.par
#' plot(est.par)
#'  
#' #  Fitted density curve and histogram
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dGamma(den.x,shape=est.par$shape,scale=est.par$scale)
#' hist(X, breaks=10, probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue")
#' lines(density(X), lty=2)
#'  
#' # Extracting shape or scale parameters
#' est.par[attributes(est.par)$par.type=="shape"]
#' est.par[attributes(est.par)$par.type=="scale"]
#' 
#' # Parameter estimation for a distribution with unknown shape parameters
#' # Example from:  Bury(1999) pp.225-226, parameter estimates as given by Bury are 
#' # shape = 6.40 and scale=2.54. The log-likelihood for this data given 
#' # Bury's parameter estimates is -656.7921.
#' data <- c(16, 11.6, 19.9, 18.6, 18, 13.1, 29.1, 10.3, 12.2, 15.6, 12.7, 13.1, 
#'          19.2, 19.5, 23, 6.7, 7.1, 14.3, 20.6, 25.6, 8.2, 34.4, 16.1, 10.2, 12.3)
#' est.par <- eGamma(data, method="numerical.MLE"); est.par
#' plot(est.par)
#' 
#' # Estimates calculated by eGamma differ from those given by Bury(1999).
#' # However, eGamma's parameter estimates appear to be an improvement, due to a larger 
#' # log-likelihood of -80.68186 (as given by lGamma below). 
#' 
#' # log-likelihood
#' lGamma(data,param = est.par)
#' 
#' # Evaluating the precision of the parameter estimates by the Hessian matrix
#' H <- attributes(est.par)$nll.hessian
#' var <- solve(H)
#' se <- sqrt(diag(var));se


#' @rdname Gamma
#' @export dGamma

dGamma <-function(x, shape = 2, scale = 2, params = list(shape = 2, scale = 2),...){
    if(!missing(params)){
      shape <- params$shape
      scale <- params$scale
    }
    out = stats::dgamma(x, shape, scale)
    return(out)
  }

#' @rdname Gamma
#' @export pGamma

pGamma <- function(q, shape = 2, scale = 2, params = list(shape = 2, scale = 2),...){
    if(!missing(params)){
      shape <- params$shape
      scale <- params$scale
    }
    out = stats::pgamma(q,shape,scale)
    return(out)
}

#' @rdname Gamma
#' @export qGamma

qGamma <- function(p, shape = 2, scale = 2, params = list(shape = 2, scale = 2),...){
    if(!missing(params)){
      shape <- params$shape
      scale <- params$scale
    }
    out = stats::qgamma(p,shape,scale)
    return(out)
}

#' @rdname Gamma
#' @export rGamma
rGamma <- function(n, shape = 2, scale = 2, params = list(shape = 2, scale = 2),...){
    if(!missing(params)){
      shape <- params$shape
      scale <- params$scale
    }
    out = stats::rgamma(n,shape,scale)
    return(out)
  }

#' @rdname Gamma
#' @export eGamma
eGamma <- function(X,w, method =c("moments","numerical.MLE"),...){
    method <- match.arg(method)
    # Accounting for weights
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    if(method == "moments") {
      # Weighted parameter estimation 
      Smean <- mean(X)
      Ssd <- sd(X)
      
      shape <- (Smean/Ssd)^2    # Engineering statistics handbook (http://www.itl.nist.gov/div898/handbook/eda/section3/eda366b.htm)
      scale <- Smean/(Ssd^2)
      
      # Standard errors of estimates 
      SE.shape <- Ssd*sqrt(trigamma(shape)/(n*(shape*trigamma(shape)-1))) # Bury pg 216
      SE.scale <- sqrt(shape/(n*(shape*trigamma(shape)-1)))
     
      est.par <- list(shape = shape, scale = scale)
      est.par.se <- c(SE.shape, SE.scale)
    } 
   else 
	{method <- "numerical.MLE"   
      est.par <- wmle(X=X, w=w, distname = "Gamma",
                      initial=list(shape = 1, scale = 1/mean(X) ),
                      lower=list(shape = 0, scale = 0),
                      upper=list(shape = Inf, scale = Inf))

      est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
      if(class(est.par.se) == "try-error") {
        est.par.se <- rep(NA, length(est.par))
      } 
    }
    
    attributes(est.par)$ob <- X
    attributes(est.par)$weights <- w
    attributes(est.par)$distname <- "Gamma"
    attributes(est.par)$method <- method
    attributes(est.par)$par.name <- c("shape","scale")
    attributes(est.par)$par.type <- c("shape","scale")
    attributes(est.par)$par.vals <- c(est.par$shape, est.par$scale)
    attributes(est.par)$par.s.e <-  est.par.se  
    
    class(est.par) <- "eDist"
    
    return(est.par)
  }

#' @rdname Gamma
#' @export lGamma
## (weighted) (log) likelihood function
lGamma <- 
  function(X, w, shape = 2, scale = 2, params = list(shape = 2, scale = 2), logL = TRUE,...){
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
    
    ll <- sum(w*log(dGamma(x=X,params = params)))
    # ll <- (shape-1)*sum(log(X)) - sum(X/scale) - n*shape*log(scale) + n*log(gamma(shape))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }
