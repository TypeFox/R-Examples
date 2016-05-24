#' @title The Burr Distribution.
#' @description Density, distribution, quantile, random number 
#' generation, and parameter estimation functions for the Burr distribution with parameters \code{location}, 
#' \code{scale} and \code{inequality}. Parameter estimation can be based on a weighted or unweighted i.i.d sample
#'  and can be performed numerically. 
#' @rdname Burr
#' @name Burr
#' @aliases dBurr
#' @aliases pBurr
#' @aliases qBurr
#' @aliases rBurr
#' @aliases eBurr
#' @aliases lBurr
#' 
#' @details The Burr distribution is a special case of the Pareto(IV) distribution 
#' where the \code{location} parameter is equal \eqn{0} and \code{inequality} parameter is equal to \eqn{1/g}, Brazauskas (2003). \cr 
#' \cr
#' The \code{dBurr()}, \code{pBurr()}, \code{qBurr()},and \code{rBurr()} functions serve as wrappers of the
#' \code{\link[VGAM]{dparetoIV}}, \code{\link[VGAM]{pparetoIV}}, \code{\link[VGAM]{qparetoIV}}, and 
#' \code{\link[VGAM]{rparetoIV}} functions in the \pkg{{VGAM}} package. They allow for the parameters to be declared not only as 
#' individual numerical values, but also as a list so parameter estimation can be carried out. \cr  
#' \cr
#' The Burr distribution is most simply defined in terms of its cumulative distribution function (Johnson et.al p.576) 
#' \deqn{F(x) = [1 + (x/b)^g]^{-s}}
#' where \eqn{b}, \eqn{g} and \eqn{s} \eqn{> 0}. Parameter estimation can only be implemented numerically.\cr
#' \cr
#' The log-likelihood and score functions are as given by Watkins (1999) and the information matrix is as given by
#' Brazauskas (2003).
#' @param params A list that includes all named parameters.
#' @param x,q A vector of quantiles.
#' @param w An optional vector of sample weights.
#' @param p A vector of probabilities.
#' @param n Number of observations.
#' @param X Sample observations.
#' @param b Scale parameters.
#' @param g,s Shape parameters.
#' @param method Parameter estimation method.
#' @param logL logical; if TRUE, lBurr gives the log-likelihood, otherwise the likelihood is given. 
#' @param ... Additional parameters.
#'
#' @return dBurr gives the density, pBurr the distribution function,
#' qBurr the quantile function, rBurr generates random deviates, and 
#' eBurr estimate the distribution parameters. lBurr provides the log-likelihood function. 
#'
#' @references Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994) Continuous Univariate Distributions,
#'  volume 1, chapter 20, Wiley, New York.\cr
#'  \cr 
#' Brazauskas, V. (2003) Information matrix for Pareto(IV), Burr, and related distributions. Comm. Statist.
#' Theory and Methods 32, 315-325.\cr
#' \cr
#' Watkins A.J. (1999) An algorithm for maximum likelihood estimation in the three parameter Burr XII distribution,
#' Computational Statistics & Data Analysis, 32, 19-27.\cr
#' \cr
#' \href{http://au.mathworks.com/help/stats/burr-type-xii-distribution.html}{Mathworks: Matlab documentation for Burr Type XII distribution}
#' @seealso \pkg{\link{ExtDist}} for other standard distributions.
#' @author Haizhen Wu and A. Jonathan R. Godfrey. \cr
#' Updates and bug fixes by Sarah Pirikahu.
#' 
#' @examples
#' # Parameter estimation for a distribution of known shape parameters
#' X <- rBurr(n=500, b = 1, g = 2, s = 2)
#' est.par <- eBurr(X); est.par
#' plot(est.par)
#' 
#' # Fitted density curve and histogram
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dBurr(den.x, b=est.par$b, g=est.par$g, s=est.par$s)
#' hist(X, breaks=10, probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue")
#' lines(density(X), lty=2)
#' 
#' # Extracting shape or scale parameters
#' est.par[attributes(est.par)$par.type=="scale"]
#' est.par[attributes(est.par)$par.type=="shape"]
#' 
#' # Parameter Estimation for a distribution with unknown shape parameters
#' # Example from: Matlab Statistical Toolbox package 
#' # (See: http://au.mathworks.com/help/stats/burr-type-xii-distribution.html)
#' # Parameter estimates given are: b = 80.4515, g = 18.9251 and s = 0.4492.
#' QRS.duration <- c(91,81,138,100,88,100,77,78,84,89,102,77,78,91,77,75,82,70,91,82,83,90,71,75,82,
#'                  109,94,95,90,96,85,71,75,78,82,69,103,85,80,94,80,79,92,84,86,73,75,73,78,80,81,
#'                  83,103,92,88,77,79,90,91,83,80,78,76,82,81,80,82,71,73,87,76,101,93,90,87,88,94,
#'                  94,90,78,83,92,93,100,83,163,96,114,170,137,84,82,79,72,97,87,102,85,84,78,79,91,
#'                  98,86,72,97,82,78,97,94,82,78,79,87,93,75,106,96,88,90,74,85,90,71,75,77,87,95,
#'                  74,99,89,83,78,100,80,87,79,102,80,85,81,85,95,82,97,92,102,86,80,85,85,111,89,63,
#'                  70,92,75,93,83,84,91,81,113,92,81,74,78,80,82,95,106,95,100,90,88,71,78,77,87,79,
#'                  85,91,92,98,68,84,92,110,108,153,73,81,87,87,95,73,95,100,96,97,76,62,86,71,99,68,
#'                  90,146,86,80,90,93,91,111,89,79,77,73,92,98,78,87,98,84,82,90,85,71,84,85,77,93,
#'                  74,89,89,103,85,88,81,84,96,90,98,78,93,80,85,67,74,69,105,95,87,108,99,79,86,82,
#'                  91,93,80,84,90,81,90,78,98,81,90,85,79,61,90,79,83,84,78,86,72,87,91,102,80,82,104,
#'                  85,83,81,94,84,91,99,101,132,79,103,97,131,91,90,121,78,84,97,94,96,91,80,97,92,90,
#'                  90,123,105,85,77,83,92,85,96,69,88,84,107,91,74,89,109,80,83,92,100,113,105,99,84,
#'                  74,76,87,87,96,88,80,85,90,74,95,86,74,95,74,73,104,92,105,97,101,83,84,98,81,93,
#'                  84,102,94,91,100,92,94,98,146,84,77,82,84,76,106,70,87,118,86,82,96,89,93,82,97,
#'                  86,188,93,72,107,81,76,83,147,82,110,108,82,93,95,80,185,73,78,71,86,85,76,93,
#'                  87,96,86,78,87,80,98,75,78,82,94,83,94,140,87,55,133,83,77,123,79,88,80,88,79,
#'                  77,87,88,94,88,74,85,88,81,91,81,80,100,108,93,79)
#' est.par <- eBurr(QRS.duration); est.par
#' plot(est.par)
#' 
#' # log-likelihood function
#' lBurr(QRS.duration,param = est.par)
#' 
#' # Evaluation of the precision of the parameter estimates by the Hessian matrix
#' H <- attributes(est.par)$nll.hessian
#' var <- solve(H)
#' se <- sqrt(diag(var)); se

#' @rdname Burr
#' @export dBurr
dBurr <-function(x, b = 1, g = 2, s = 2, params = list(b = 1, g = 2, s = 2),...){
    if(!missing(params)){
      b = params$b; g = params$g; s = params$s
    }
    VGAM::dparetoIV(x, location = 0, scale = b, inequality = 1/g, shape = s)  }

#' @rdname Burr
#' @export pBurr
pBurr <- function(q, b = 1, g = 2, s = 2, params = list(b = 1, g = 2, s = 2),...){
    if(!missing(params)){
      b = params$b; g = params$g; s = params$s
    }
    q[q<=0] =0
    out = q
    c = (q/b)^g
    d = s*c
    out[q == Inf] =1  
    out[c>1e10 & q!=Inf] = 1- (q[c>1e10 & q!=Inf]/b)^(-g*s)
    out[c<1e-10 & q!=Inf] = 1- exp(-d[c<1e-10 & q!=Inf])
    out[c<=1e10 & c>=1e-10 & q!=Inf] = VGAM::pparetoIV(q[c<=1e10 & c>=1e-10 & q!=Inf], location = 0, scale = b, inequality = 1/g, shape = s)
    return(out)
  }

#' @rdname Burr
#' @export qBurr
qBurr <- function(p, b = 1, g = 2, s = 2, params = list(b = 1, g = 2, s = 2),...){
    if(!missing(params)){
      b = params$b; g = params$g; s = params$s
    }
    out = p
    c= (1/(1-p))^(1/s)
    d = -log(1-p)/s
    
    out[p== 0] =0
    
    c= (1/(1-p))^(1/s)
    out[c>1e10 & p!=0] = (1/(1-p[c>1e10 & p!=0]))^(1/(s*g))*b
    
    d = -log(1-p)/s
    out[d<1e-10 & p!=0] = b*(d[d<1e-10 & p!=0])^(1/g) 
    out[c<=1e10 & d>=(1e-10) & p!=0] = VGAM::qparetoIV(p[c<=1e10 & d>=(1e-10) & p!=0], location = 0, scale = b, inequality = 1/g, shape = s)
    return(out)
  }

#' @rdname Burr
#' @export rBurr
rBurr <- function(n, b = 1, g = 2, s = 2, params = list(b = 1, g = 2, s = 2),...){
    if(!missing(params)){
      b = params$b; g = params$g; s = params$s
    }
    VGAM::rparetoIV(n, location = 0, scale = b, inequality = 1/g, shape = s)
  }

#' @rdname Burr
#' @export eBurr
eBurr <- function(X,w, method ="numerical.MLE",...){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
{ if(method != "numerical.MLE") stop(paste("method ", method, " is not avaliable, use numerial.MLE instead."))  
  method = "numerical.MLE"  
  
  est.par <- wmle(X=X, w=w, distname = "Burr",
                  initial=list(b=2, g=2, s =2),
                  lower=list(b=0, g=0, s =0),
                  upper=list(b=Inf, g=Inf, s =Inf))
  
  est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
  if(class(est.par.se) == "try-error") {
    est.par.se <- rep(NA, length(est.par))
  } 
}

attributes(est.par)$ob <- X
attributes(est.par)$weights <- w
attributes(est.par)$distname <- "Burr"
attributes(est.par)$method <- method
attributes(est.par)$par.name <- c("b","g","s")
attributes(est.par)$par.type <- c("scale","shape","shape")
attributes(est.par)$par.vals <- c(est.par$b, est.par$g, est.par$s)
attributes(est.par)$par.s.e <-  est.par.se  

class(est.par) <- "eDist"

return(est.par)
  }

#' @rdname Burr
#' @export lBurr
## (weighted) (log) likelihood function
lBurr <- function(X, w, b = 1, g = 2, s = 2, params = list(b = 1, g = 2, s = 2), logL = TRUE,...){
    if(!missing(params)){
      b = params$b; g = params$g; s = params$s
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dBurr(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }