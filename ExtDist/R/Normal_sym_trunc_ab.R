#' @title The symmetric truncated normal distribution.
#' @description Density, distribution, quantile, random number 
#' generation and parameter estimation functions for the symmetric truncated normal distribution with parameters, \code{sigma}, 
#' \code{a} and \code{b} which represent the lower and upper truncation points respectively.  
#' Parameter estimation can be based on a weighted or unweighted i.i.d sample and can be carried out numerically.
#' 
#' @rdname Normal_sym_trunc_ab
#' @name Normal_sym_trunc_ab
#' @aliases dNormal_sym_trunc_ab 
#' @aliases pNormal_sym_trunc_ab 
#' @aliases qNormal_sym_trunc_ab 
#' @aliases rNormal_sym_trunc_ab 
#' @aliases eNormal_sym_trunc_ab
#' @aliases lNormal_sym_trunc_ab
#'
#' @details The normal symmetric truncated distribution is a special case of the trucated normal distribution. 
#' See \code{\link{Normal_trunc_ab}}.
#' 
#' @param params A list that includes all named parameters.
#' @param x,q A vector of quantiles.
#' @param w An optional vector of sample weights.
#' @param p A vector of probabilities.
#' @param n Number of observations.
#' @param X Sample observations.
#' @param a,b Boundary parameters.
#' @param mu,sigma Shape parameters.
#' @param method Parameter estimation method.
#' @param logL logical;if TRUE, lNormal_sym_trunc_ab gives the log-likelihood, otherwise the likelihood is given.
#' @param ... Additional parameters
#'
#' @return dNormal_sym_trunc_ab gives the density, pNormal_sym_trunc_ab the distribution function,
#' qNormal_sym_trunc_ab the quantile function, rNormal_sym_trunc_ab generates random deviates,and 
#' eNormal_sym_trunc_ab estimates the parameters. lNormal_sym_trunc_ab provides the log-likelihood function. 
#' @seealso \pkg{\link{ExtDist}} for other standard distributions.
#' @author Haizhen Wu and A. Jonathan R. Godfrey.
#'

#' @rdname Normal_sym_trunc_ab
#' @export dNormal_sym_trunc_ab
dNormal_sym_trunc_ab <- function(x, sigma=0.3, a= 0, b=1, params = list(sigma, a, b),...){
    if(!missing(params)){
      mu <- params$mu; sigma <- params$sigma; a <- params$a; b <- params$b
    }
    out <- truncdist::dtrunc( x, spec="norm", a=a, b=b, mean = (a+b)/2, sd = sigma)
    return(out)
  }

#' @rdname Normal_sym_trunc_ab
#' @export pNormal_sym_trunc_ab
pNormal_sym_trunc_ab <- function(q, sigma=0.3, a= 0, b=1, params = list(mu=2, sigma = 5, a = 0, b = 1),...){
    if(!missing(params)){
      mu <- params$mu; sigma <- params$sigma; a <- params$a; b <- params$b
    }
    out <- truncdist::ptrunc( q, spec="norm", a=a, b=b, mean = (a+b)/2, sd = sigma)
    return(out)
  }

#' @rdname Normal_sym_trunc_ab
#' @export qNormal_sym_trunc_ab
qNormal_sym_trunc_ab <- function(p, sigma=0.3, a= 0, b=1, params = list(mu=2, sigma = 5, a = 0, b = 1),...){
    if(!missing(params)){
      mu <- params$mu; sigma <- params$sigma; a <- params$a; b <- params$b
    }
    out <- truncdist::qtrunc( p, spec="norm", a=a, b=b, mean = (a+b)/2, sd = sigma)
    return(out)
  }

#' @rdname Normal_sym_trunc_ab
#' @export rNormal_sym_trunc_ab
rNormal_sym_trunc_ab <- function(n, mu=2, sigma=3, a = 0, b = 1, params = list(sigma, a, b),...){
    if(!missing(params)){
      mu <- params$mu; sigma <- params$sigma; a <- params$a; b <- params$b
    }
    out <- truncdist::rtrunc(n, spec="norm", a=a, b=b, mean = (a+b)/2, sd = sigma)
    return(out)
  }

#' @rdname Normal_sym_trunc_ab
#' @export eNormal_sym_trunc_ab
eNormal_sym_trunc_ab <- function(X,w, method ="numerical.MLE",...){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
{
  if(method != "numerical.MLE") warning(paste("method ", method, " is not avaliable, use numerial.MLE instead."))  
  method = "numerical.MLE"  
  
  d <- (max(X)-min(X))
  est.par <- wmle(X=X, w=w, distname = "Normal_sym_trunc_ab",
                  initial=list(sigma=sd(X),a=min(X)-0.1*d,b=max(X)+0.1*d),
                  lower=list(sigma=0,a=-Inf,b=max(X)),
                  upper=list(sigma=Inf,a=min(X),b=Inf)
  )
  
  est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
  if(class(est.par.se) == "try-error") {
    est.par.se <- rep(NA, length(est.par))
  } 
} 

attributes(est.par)$ob <- X
attributes(est.par)$weights <- w
attributes(est.par)$distname <- "Normal_sym_trunc_ab"
attributes(est.par)$method <- method
attributes(est.par)$par.name <- c("mu","sigma","a","b")
attributes(est.par)$par.type <- c("shape","shape","boundary","boundary")
attributes(est.par)$par.vals <- c(est.par$mu, est.par$sigma, est.par$a, est.par$b)
attributes(est.par)$par.s.e <-  est.par.se  

class(est.par) <- "eDist"

return(est.par)
  }

#' @rdname Normal_sym_trunc_ab
#' @export lNormal_sym_trunc_ab
## (weighted) (log) likelihood function
lNormal_sym_trunc_ab <- function(X, w, mu=2, sigma =3, a = 0, b = 1,  params = list(sigma, a, b), logL = TRUE,...){
    if(!missing(params)){
      mu <- params$mu; sigma <- params$sigma; a <- params$a; b <- params$b
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dNormal_sym_trunc_ab(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }