#' @title The Triangular Distribution.
#' 
#' @description  Density, distribution, quantile, random number 
#' generation and parameter estimation functions for the triangular distribution with support \eqn{[a,b]} 
#' and \code{shape} parameter \eqn{\theta}. Parameter estimation can be based on a weighted or unweighted i.i.d. sample
#' and can be performed numerically. 
#' 
#' @rdname Triangular
#' @name Triangular
#' 
#' @aliases dTriangular
#' @aliases pTriangular
#' @aliases qTriangular
#' @aliases rTriangular
#' @aliases eTriangular
#' @aliases lTriangular
#'
#' @details If \code{a}, \code{b} or \code{theta} are not specified they assume the default values of 0, 1 and 0.5 respectively.\cr
#' \cr
#' The \code{dTriangle()}, \code{pTriangle()}, \code{qTriangle()},and \code{rTriangle()} functions serve as wrappers of the
#' \code{\link[VGAM]{dtriangle}}, \code{\link[VGAM]{ptriangle}}, \code{\link[VGAM]{qtriangle}}, and 
#' \code{\link[VGAM]{rtriangle}} functions in the \pkg{{VGAM}} package. They allow for the parameters to be declared not only as 
#' individual numerical values, but also as a list so parameter estimation can be carried out. \cr  
#' \cr
#' The triangular distribution has a probability density function, defined in Forbes et.al (2010), that consists of two lines joined at \eqn{theta}, 
#' where \eqn{theta} is the location of the mode.
#' 
#' @param params A list that includes all named parameters.
#' @param x,q A vector of quantiles.
#' @param w An optional vector of sample weights.
#' @param p A vector of probabilities.
#' @param n Number of observations.
#' @param X Sample observations.
#' @param theta Shape parameters.
#' @param a,b Boundary parameters.
#' @param method Parameter estimation method.
#' @param logL logical, it is assumed that the log-likelihood is desired. Set to FALSE if the likelihood is wanted.
#' @param ... Additional parameters.
#'
#' @return dTriangular gives the density, pTriangular the distribution function,
#' qTriangular the quantile function, rTriangular generates random variables, and 
#' eTriangular estimates the parameters. lTriangular provides the log-likelihood function.
#' @seealso \pkg{\link{ExtDist}} for other standard distributions.
#' @author Haizhen Wu and A. Jonathan R. Godfrey. \cr
#' Updates and bug fixes by Sarah Pirikahu.
#'
#' @references Kotz, S. and van Dorp, J. R. (2004). Beyond Beta: Other Continuous
#' Families of Distributions with Bounded Support and Applications. Chapter 1.
#' World Scientific: Singapore.\cr
#' \cr
#' Forbes, C., Evans, M., Hastings, N. and Peacock, B. (2010) Triangular Distribution, 
#' in Statistical Distributions, Fourth Edition, John Wiley & Sons, Inc., Hoboken, NJ, USA. 
#' 

#' @rdname Triangular
#' @export dTriangular
dTriangular <-function(x, a=0, b=1, theta=0.5, params = list(a, b, theta),...){
    if(!missing(params)){
      a <- params$a; b <- params$b; theta = params$theta  
    }
    out <- VGAM::dtriangle(x, theta = theta*(b-a)+a, lower = a, upper = b)
    return(out)
  }

#' @rdname Triangular
#' @export pTriangular
pTriangular <-function(q, a=0, b=1, theta=0.5, params = list(a, b, theta),...){
    if(!missing(params)){
      a <- params$a; b <- params$b; theta = params$theta  
    }
    out <- VGAM::ptriangle(q, theta = theta*(b-a)+a, lower = a, upper = b)
    return(out)
  }

#' @rdname Triangular
#' @export qTriangular
qTriangular <-function(p, a=0, b=1, theta=0.5, params = list(a, b, theta),...){
    if(!missing(params)){
      a <- params$a; b <- params$b; theta = params$theta  
    }
    out <- p
    out[p==0] <- a
    out[p!=0] <- VGAM::qtriangle(p[p!=0], theta = theta*(b-a)+a, lower = a, upper = b)
    return(out)
  }

#' @rdname Triangular
#' @export rTriangular
rTriangular <-function(n, a=0, b=1, theta=0.5, params = list(a, b, theta),...){
    if(!missing(params)){
      a <- params$a; b <- params$b; theta = params$theta  
    }
    out <- VGAM::rtriangle(n, theta = theta*(b-a)+a, lower = a, upper = b)
    return(out)
  }

#' @rdname Triangular
#' @export eTriangular
eTriangular <-function(X,w, method ="numerical.MLE",...){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
{
  if(method != "numerical.MLE") warning(paste("method ", method, " is not avaliable, use numerial.MLE instead."))  
  method = "numerical.MLE"  
  
  d <- max(X)-min(X)
  est.par <- wmle(X=X, w=w, distname = "Triangular",
            initial=list(a=min(X)-0.1*d,b=max(X)+0.1*d, theta= 0.5),
            lower=list( a=-Inf,b=max(X), theta= 1e-10),
            upper=list(a=min(X),b=Inf, theta= 1- 1e-10))
  
  est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
  if(class(est.par.se) == "try-error") {
    est.par.se <- rep(NA, length(est.par))
  } 
} 

attributes(est.par)$ob <- X
attributes(est.par)$weights <- w
attributes(est.par)$distname <- "Triangular"
attributes(est.par)$method <- method
attributes(est.par)$par.name <- c("a","b","theta")
attributes(est.par)$par.type <- c("boundary","boundary","shape")
attributes(est.par)$par.vals <- c(est.par$a, est.par$b, est.par$theta)
attributes(est.par)$par.s.e <-  est.par.se  

class(est.par) <- "eDist"

return(est.par)
  }

#' @rdname Triangular
#' @export lTriangular
## (weighted) (log) likelihood function
lTriangular <-function(X, w, a=0, b=1, theta=0.5,  params = list(a, b, theta), logL = TRUE,...){
    if(!missing(params)){
      a <- params$a; b <- params$b; theta = params$theta  
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dTriangular(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }