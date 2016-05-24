#' @title The standard symmetric-reflected truncated beta (SSRTB) distribution.
#' @description Density, distribution, quantile, random number 
#' generation and parameter estimation functions for the SSRTB distribution. Parameter estimation can be
#' based on a weighted or unweighted i.i.d sample and can be carried out numerically. 
#' 
#' @rdname SSRTB
#' @name SSRTB
#' 
#' @aliases dSSRTB
#' @aliases pSSRTB
#' @aliases qSSRTB
#' @aliases rSSRTB 
#' @aliases eSSRTB
#' @aliases lSSRTB 
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
#' @param method Parameter estimation method.
#' @param logL logical; if TRUE, lSSRTB gives the log-likelihood, otherwise the likelihood is given.
#' @param ... Additional parameters.
#' 
#' @return dSSRTB gives the density, pSSRTB the distribution function,
#' qSSRTB the quantile function, rSSRTB generates random variables, 
#' eSSRTB estimates the parameters and lSSRTB provides the log-likelihood.
#'
#' @seealso \pkg{\link{ExtDist}} for other standard distributions.
#' @author Haizhen Wu.
#'
#' @examples
#' # Parameter estimation for a distribution with known shape parameters
#' X <- rSSRTB(n=500, shape1=2, shape2=10)
#' est.par <- eSSRTB(X); est.par
#' plot(est.par)
#' 
#' #  Fitted density curve and histogram
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dSSRTB(den.x,shape1=est.par$shape1,shape2=est.par$shape2)
#' hist(X, breaks=10, probability=TRUE, ylim = c(0,1.2*max(den.y)))
#' lines(den.x, den.y, col="blue")
#' lines(density(X), lty=2)
#' 
#' # Extracting shape parameters
#' est.par[attributes(est.par)$par.type=="shape"]
#' 
#' # log-likelihood function
#' lSSRTB(X,param = est.par)


#' @rdname SSRTB
#' @export dSSRTB
dSSRTB <-
  function(x, shape1=2, shape2 =3, params = list(shape1, shape2),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    m = (1-1/shape1)/(1+(shape2-2)/shape1)
    out = m*dbeta(2*m*(1/2-abs(x-1/2)),shape1,shape2)/pbeta(m,shape1,shape2)

    return(out)
  }


#' @rdname SSRTB
#' @export pSSRTB
pSSRTB <- 
  function(q, shape1=2, shape2 =3, params = list(shape1, shape2),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    m = (1-1/shape1)/(1+(shape2-2)/shape1)
    out= ifelse(q<1/2, 
                pbeta(2*m*q,shape1,shape2)/2/pbeta(m,shape1,shape2), 
                1-pbeta(2*m*(1-q),shape1,shape2)/2/pbeta(m,shape1,shape2))

    return(out)
  }

#' @rdname SSRTB
#' @export qSSRTB
qSSRTB <- 
  function(p, shape1=2, shape2 =3, params = list(shape1, shape2),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    m = (1-1/shape1)/(1+(shape2-2)/shape1)
    out = ifelse(p<1/2, qbeta((p<1/2)*p*2*pbeta(m,shape1,shape2),shape1,shape2)/2/m,
                 1-qbeta((p>=1/2)*(1-p)*2*pbeta(m,shape1,shape2),shape1,shape2)/2/m)
    return(out)
  }

#' @rdname SSRTB
#' @export rSSRTB
rSSRTB <- 
  function(n, shape1=2, shape2 =3, params = list(shape1, shape2),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    qSSRTB(runif(n),shape1,shape2)
  }

#' @rdname SSRTB
#' @export eSSRTB
eSSRTB <-     
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
      
      est.par <- wmle(X=X, w=w, distname = "SSRTB",
                      initial=list(shape1 = 5, shape2 = 5),
                      lower=list(shape1 = 1.1, shape2 = 1),
                      upper=list(shape1 = Inf, shape2 = Inf))
      
      est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
      if(class(est.par.se) == "try-error") {
        est.par.se <- rep(NA, length(est.par))
      } 
    }
    
    attributes(est.par)$ob <- X
    attributes(est.par)$weights <- w
    attributes(est.par)$distname <- "SSRTB"
    attributes(est.par)$method <- method
    attributes(est.par)$par.name <- c("shape1","shape2")
    attributes(est.par)$par.type <- c("shape","shape")
    attributes(est.par)$par.vals <- c(est.par$shape1, est.par$shape2)
    attributes(est.par)$par.s.e <-  est.par.se  
    
    class(est.par) <- "eDist"
    
    return(est.par)
  }

#' @rdname SSRTB
#' @export lSSRTB
## (weighted) (log) likelihood function
lSSRTB <- 
  function(X, w, shape1=2, shape2 =3, params = list(shape1, shape2), logL = TRUE,...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dSSRTB(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }