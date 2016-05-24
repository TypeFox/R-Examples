#' @title The symmetric-reflected truncated beta (SRTB) distribution.
#' @description Density, distribution, quantile, random number 
#' generation and parameter estimation functions for the SRTB distribution.  
#' Parameter estimation can be based on a weighted or unweighted i.i.d. sample and can be carried out numerically.
#' 
#' @rdname SRTB_ab
#' @name SRTB_ab
#' 
#' @aliases dSRTB_ab
#' @aliases pSRTB_ab
#' @aliases qSRTB_ab
#' @aliases rSRTB_ab
#' @aliases eSRTB_ab
#' @aliases lSRTB_ab 
#' @aliases sSRTB_ab
#' 
#' @details No details as of yet.
#' 
#' @param params A list that includes all named parameters.
#' @param x,q A vector of quantiles.
#' @param w An optional vector of sample weights.
#' @param p A vector of probabilities.
#' @param n Number of observations.
#' @param X Sample observations.
#' @param shape1,shape2 Shape parameters.
#' @param a,b Boundary parameters.
#' @param method Parameter estimation method.
#' @param logL logical; if TRUE, lSRTB_ab gives the log-likelihood, otherwise the likelihood is given.
#' @param ... Additional parameters.
#'
#' @return dSRTB_ab gives the density, pSRTB_ab the distribution function,
#' qSRTB_ab gives the quantile function, rSRTB_ab generates random variables, and
#' eSRTB_ab estimates the parameters. lSRTB_ab provides the log-likelihood function and sSRTB_ab the score function. 
#' @seealso \pkg{\link{ExtDist}} for other standard distributions.
#' @author Haizhen Wu.
#'
#' @examples 
#' # Parameter estimation for a distribution with known shape parameters
#' X <- rSRTB_ab(n=500, shape1=2, shape2=10, a=1, b=2)
#' est.par <- eSRTB_ab(X)
#' plot(est.par) 
#' 
#' # Extracting boundary and shape parameters
#' est.par[attributes(est.par)$par.type=="boundary"]
#' est.par[attributes(est.par)$par.type=="shape"]
#' 
#' # log-likelihood function
#' lSRTB_ab(X,param = est.par)


#' @rdname SRTB_ab
#' @export dSRTB_ab
dSRTB_ab <-function(x, shape1=2, shape2=3, a = 0, b=1, params = list(shape1, shape2, a, b),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
      a <- params$a
      b <- params$b
    }
    out = dSSRTB((x-a)/(b-a),shape1,shape2)/(b-a)
    return(out)
  }

#' @rdname SRTB_ab
#' @export pSRTB_ab
pSRTB_ab <- 
  function(q, shape1=2, shape2=3, a = 0, b=1, params = list(shape1=2, shape2 = 5, a = 0, b = 1),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
      a <- params$a
      b <- params$b
    }
    out = pSSRTB((q-a)/(b-a),shape1,shape2)
    return(out)
  }

#' @rdname SRTB_ab
#' @export qSRTB_ab
qSRTB_ab <- 
  function(p, shape1=2, shape2=3, a = 0, b=1, params = list(shape1=2, shape2 = 5, a = 0, b = 1),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
      a <- params$a
      b <- params$b
    }
    out = (b-a)*qSSRTB(p,shape1,shape2) + a
    return(out)
  }

#' @rdname SRTB_ab
#' @export rSRTB_ab
rSRTB_ab <- 
  function(n, shape1=2, shape2=3, a = 0, b = 1, params = list(shape1, shape2, a, b),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
      a <- params$a
      b <- params$b
    }
    x = rSSRTB(n,shape1,shape2)
    out = y = (b-a)*x + a
    return(out)
  }

#' @rdname SRTB_ab
#' @export eSRTB_ab
eSRTB_ab <-     
  function(X,w, method ="numerical.MLE",...){
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
  est.par <- wmle(X=X, w=w, distname = "SRTB_ab",
                  initial=list(shape1=3,shape2=3,a=min(X)-0.1*d,b=max(X)+0.1*d),
                  lower=list(shape1=1,shape2=1,a=-Inf,b=max(X)),
                  upper=list(shape1=Inf,shape2=Inf,a=min(X),b=Inf))
  
  est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
  if(class(est.par.se) == "try-error") {
    est.par.se <- rep(NA, length(est.par))
  } 
} 

attributes(est.par)$ob <- X
attributes(est.par)$weights <- w
attributes(est.par)$distname <- "SRTB_ab"
attributes(est.par)$method <- method
attributes(est.par)$par.name <- c("shape1","shape2","a","b")
attributes(est.par)$par.type <- c("shape","shape","boundary","boundary")
attributes(est.par)$par.vals <- c(est.par$shape1, est.par$shape2, est.par$a, est.par$b)
attributes(est.par)$par.s.e <-  est.par.se  

class(est.par) <- "eDist"

return(est.par)
  }

#' @rdname SRTB_ab
#' @export lSRTB_ab
## (weighted) (log) likelihood function
lSRTB_ab <- 
  function(X, w, shape1=2, shape2 =3, a = 0, b = 1,  params = list(shape1, shape2, a, b), logL = TRUE,...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
      a <- params$a
      b <- params$b
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dSRTB_ab(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }