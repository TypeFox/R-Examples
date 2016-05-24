#' @title The Johnson SB distribution.
#' @description Density, distribution, quantile, random number 
#' generation, and parameter estimation functions for the Johnson SB (bounded support) distribution.
#' Parameter estimation can be based on a weighted or unweighted i.i.d. sample and can be performed numerically.
#' 
#' @rdname JohnsonSB
#' @name JohnsonSB
#' 
#' @aliases dJohnsonSB
#' @aliases pJohnsonSB
#' @aliases qJohnsonSB
#' @aliases rJohnsonSB
#' @aliases eJohnsonSB
#' @aliases lJohnsonSB
#' 
#' @details The Johnson system of distributions consists of families of distributions that, through specified transformations, can be 
#' reduced to the standard normal random variable. It provides a very flexible system for describing statistical distributions
#' and is defined by 
#' \deqn{z = \gamma + \delta f(Y)} 
#' with \eqn{Y = (X-xi)/lambda}. The Johnson SB distribution arises when \eqn{f(Y) = ln[Y/(1-Y)]}, where \eqn{0 < Y < 1}. 
#' This is the bounded Johnson family since the range of Y is \eqn{(0,1)}, Karian & Dudewicz (2011).\cr
#' \cr
#' The \code{dJohnsonSB()}, \code{pJohnsonSB()}, \code{qJohnsonSB()},and \code{rJohnsonSB()} functions serve as wrappers of the
#' \code{\link[SuppDists]{dJohnson}}, \code{\link[SuppDists]{pJohnson}}, \code{\link[SuppDists]{qJohnson}}, and 
#' \code{\link[SuppDists]{rJohnson}} functions in the \pkg{{SuppDists}} package. They allow for the parameters to be declared not only as 
#' individual numerical values, but also as a list so parameter estimation can be carried out. \cr  
#' \cr
#' The JohnsonSB distribution has probability density function
#' \deqn{p_X(x) = \frac{\delta lambda}{\sqrt{2\pi}(x-xi)(1- x + xi)}exp[-0.5(\gamma + \delta ln((x-xi)/(1-x+xi)))^2].} 
#' 
#' @param params A list that includes all named parameters.
#' @param x,q A vector of quantiles.
#' @param w An optional vector of sample weights.
#' @param p A vector of probabilities.
#' @param n Number of observations.
#' @param X Sample observations.
#' @param xi,lambda,a,b Location-scale parameters.
#' @param gamma,delta Shape parameters.
#' @param method Parameter estimation method.
#' @param logL logical, it is assumed that the log-likelihood is desired. Set to FALSE if the likelihood is wanted.
#' @param ... Additional parameters.
#'
#' @return dJohnsonSB gives the density, pJohnsonSB the distribution function,
#' qJohnsonSB gives quantile function, rJohnsonSB generates random deviates, and  
#' eJohnsonSB estimate the parameters. lJohnsonSB provides the log-likelihood function. The dJohnsonSB_ab provides 
#' an alternative parameterisation of the JohnsonSB distribution.
#'
#' @references Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994) Continuous Univariate Distributions,
#' volume 1, chapter 12, Wiley, New York.\cr
#' \cr
#' Kotz, S. and van Dorp, J. R. (2004). Beyond Beta: Other Continuous
#' Families of Distributions with Bounded Support and Applications. Appendix B.
#' World Scientific: Singapore.\cr
#' \cr
#' Z. A. Karian and E. J. Dudewicz (2011) Handbook of Fitting Statistical Distributions with R, Chapman & Hall.\cr 
#' 
#' @seealso \pkg{\link{ExtDist}} for other standard distributions.   
#'    
#' @author Haizhen Wu and A. Jonathan R. Godfrey.\cr
#' Updates and bug fixes by Sarah Pirikahu.
#'
#' @examples
#' # Parameter estimation for a distribution with known shape parameters
#' X <- rJohnsonSB(n=500, gamma=-0.5, delta=2, xi=-0.5, lambda=2)
#' est.par <- eJohnsonSB(X); est.par
#' plot(est.par)
#' 
#' #  Fitted density curve and histogram
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dJohnsonSB(den.x,params = est.par)
#' hist(X, breaks=10, probability=TRUE, ylim = c(0,1.2*max(den.y)))
#' lines(den.x, den.y, col="blue")
#' lines(density(X))
#'
#' # Extracting location, scale and shape parameters
#' est.par[attributes(est.par)$par.type=="location"]
#' est.par[attributes(est.par)$par.type=="scale"]
#' est.par[attributes(est.par)$par.type=="shape"]
#' 
#' # Parameter Estimation for a distribution with unknown shape parameters
#' # Example from Karian, Z.A and Dudewicz, E.J. (2011) p.647. 
#' # Original source of brain scan data Dudewich, E.J et.al (1989).
#' # Parameter estimates as given by Karian & Dudewicz using moments are:
#' # gamma =-0.2081, delta=0.9167, xi = 95.1280 and lambda = 21.4607 with log-likelihood = -67.03579
#' brain <- c(108.7, 107.0, 110.3, 110.0, 113.6, 99.2, 109.8, 104.5, 108.1, 107.2, 112.0, 115.5, 108.4, 
#'            107.4, 113.4, 101.2, 98.4, 100.9, 100.0, 107.1, 108.7, 102.5, 103.3)
#' est.par <- eJohnsonSB(brain); est.par
#' 
#' # Estimates calculated by eJohnsonSB differ from those given by Karian & Dudewicz (2011).
#' # However, eJohnsonSB's parameter estimates appear to be an improvement, due to a larger 
#' # log-likelihood of -66.35496 (as given by lJohnsonSB below).
#' 
#' # log-likelihood function
#' lJohnsonSB(brain, param = est.par)

#' @rdname JohnsonSB
#' @export dJohnsonSB
dJohnsonSB <- function(x, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2, params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2),...){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }
    out=SuppDists::dJohnson(x, parms = list(gamma, delta, xi, lambda, type= "SB"))
    out[is.nan(out)]=0
    return(out)
  }

#' @rdname JohnsonSB
#' @export dJohnsonSB_ab
dJohnsonSB_ab <-function(x, gamma = -0.5, delta = 2, a = -0.5, b = 1.5, params = list(gamma = -0.5, delta = 2, a = -0.5, b = 1.5),...){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      a <- params$a
      b <- params$b
    }
    xi <- a
    lambda <- (b-a)
    
    out = dJohnsonSB(x, params=list(gamma=gamma, delta=delta, xi=xi, lambda=lambda))
    return(out)
  }

#' @rdname JohnsonSB
#' @export pJohnsonSB
pJohnsonSB <- function(q, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2, params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2),...){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }
    out = NULL
    for( i in 1:length(q)){
      out[i] = try(SuppDists::pJohnson(q[i], parms = list(gamma, delta, xi, lambda, type= "SB")),silent=T)
      if((!is.numeric(out[i]))&(q[i]<=xi)) out[i]= 0
      if((!is.numeric(out[i]))&(q[i]>=xi+lambda)) out[i]= 1
    }
    return(as.numeric(out))
  }

#' @rdname JohnsonSB
#' @export qJohnsonSB
qJohnsonSB <- function(p, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2, params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2),...){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }
    out = lambda/(exp(-(qnorm(p)-gamma)/delta)+1)+xi
    return(out)
  }

#' @rdname JohnsonSB
#' @export rJohnsonSB
rJohnsonSB <- function(n, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2, params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2),...){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }
    out <- SuppDists::rJohnson(n, parms = list(gamma=gamma, delta=delta, xi=xi, lambda=lambda, type= "SB"))
    return(out)
  }

#' @rdname JohnsonSB
#' @export eJohnsonSB
eJohnsonSB <- function(X,w, method ="numerical.MLE",...){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    # xi and lambda cannot be fitted directly.
    
{if(method != "numerical.MLE") warning(paste("method ", method, " is not avaliable, use numerial.MLE instead."))  
 method = "numerical.MLE"  
 
#     l = min(X)
#     u = max(X)
#     Y=(X-l)/(u-l)
# 	
#  out <- wmle(X=Y, w=w, distname = "JohnsonSB_ab",
#              initial=list(gamma = -0.5, delta = 2, a = -1, b = 2),
#              lower=list(gamma = -Inf, delta = 0, a = -Inf, b = 1),
#              upper=list(gamma = Inf, delta = Inf, a = 0, b = Inf)
#  )
#     out$a = l+out$a*(u-l)
#     out$b = l+out$b*(u-l)
# 	
# 	deriv = c(1,1,1/(u-l),1/(u-l))
#     if(!is.null(attr(out,"nll.hessian"))) {    
#       attr(out,"nll.hessian") = attr(out,"nll.hessian") * deriv%*%t(deriv)
#     }
	
 out <- wmle(X=X, w=w, distname = "JohnsonSB_ab",
             initial=list(gamma = -0.5, delta = 2, a = min(X)-1, b = max(X)+1),
             lower=list(gamma = -Inf, delta = 0, a = -Inf, b = max(X)),
             upper=list(gamma = Inf, delta = Inf, a = min(X), b = Inf)
 )

 est.par = list(gamma = out$gamma , delta = out$delta, xi = out$a, lambda = out$b-out$a)
 trans.matrix = diag(rep(1,4))
 trans.matrix[4,3]=1 
 attributes(est.par)$nll.hessian = trans.matrix %*% attributes(out)$nll.hessian %*% t(trans.matrix)
 
 est.par.se <- try(suppressWarnings(sqrt(diag(solve(attributes(est.par)$nll.hessian )))),silent=TRUE)
 if(class(est.par.se) == "try-error" | any(is.nan(est.par.se))) {
   est.par.se <- rep(NA, length(est.par))
 } 
}

attributes(est.par)$ob <- X
attributes(est.par)$weights <- w
attributes(est.par)$distname <- "JohnsonSB"
attributes(est.par)$method <- method
attributes(est.par)$par.name <- c("gamma","delta","xi","lambda")
attributes(est.par)$par.type <- c("shape","shape","location","scale")
attributes(est.par)$par.vals <- c(est.par$gamma, est.par$delta, est.par$xi, est.par$lambda)
attributes(est.par)$par.s.e <-  est.par.se  

class(est.par) <- "eDist"

return(est.par)
  }

#' @rdname JohnsonSB
#' @export lJohnsonSB
## (weighted) (log) likelihood function
lJohnsonSB <- function(X, w, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2,  params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2), logL = TRUE,...){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dJohnsonSB(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }
