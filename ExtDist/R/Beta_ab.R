#' @title The four-parameter beta distribution.
#' @description Density, distribution, quantile, random number 
#' generation, and parameter estimation functions for the 4-parameter beta distribution. 
#' Parameter estimation can be based on a weighted or unweighted i.i.d sample and can be performed
#' numerically. 
#' @rdname Beta_ab
#' @name Beta_ab
#' @aliases dBeta_ab 
#' @aliases pBeta_ab
#' @aliases qBeta_ab
#' @aliases rBeta_ab
#' @aliases eBeta_ab
#' @aliases lBeta_ab
#' @aliases sBeta_ab
#' 
#' @details The \code{dBeta_ab()}, \code{pBeta_ab()}, \code{qBeta_ab()},and \code{rBeta_ab()} functions serve as wrappers of the standard
#' \code{\link[stats]{dbeta}}, \code{\link[stats]{pbeta}}, \code{\link[stats]{qbeta}} and \code{\link[stats]{rbeta}} functions 
#' in the \pkg{\link{stats}} package.They allow for the parameters to be declared not only as 
#' individual numerical values, but also as a list so parameter estimation can be carried out. \cr
#' \cr
#' The four-parameter beta distribution with parameters \code{shape1}=p, \code{shape2}=q, \code{a} = \eqn{a} and \code{b}=\eqn{b} 
#' has probability density function
#' \deqn{f(x) = \frac{1}{B(p,q)} \frac{(x-a)^{(p-1)})(b-x)^{(q-1)}}{((b-a)^{(p+q-1)}))}}  
#' with \eqn{p >0}, \eqn{q > 0}, \eqn{a \leq x \leq b} and where B is the \link[base]{beta} function, Johnson et.al (p.210). \cr
#' \cr
#' The log-likelihood function of the four-parameter beta distribution is 
#' \deqn{l(p,q,a,b| x) = -ln B(p,q) + ((p-1) ln (x-a) + (q-1) ln (b-x)) - (p + q -1) ln (b-a).} 
#' Johnson et.al (p.226) provides the Fisher's information matrix of the four-parameter beta distribution in
#' the regular case where \eqn{p,q > 2}.
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
#' @param logL logical; if TRUE, lBeta_ab gives the log-likelihood, otherwise the likelihood is given.
#' @param ... Additional parameters.
#'
#' @return dBeta_ab gives the density, pBeta_ab the distribution function,
#' qBeta_ab the quantile function, rBeta_ab generates random deviates, and
#' eBeta_ab estimates the parameters. lBeta_ab provides the log-likelihood function, sBeta_ab the observed score function
#' and iBeta_ab the observed information matrix. 
#'
#' @references Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions,
#'  volume 2, chapter 25, Wiley, New York.\cr
#'  \cr
#'  Bury, K. (1999) Statistical Distributions in Engineering, Chapter 14, pp.261-262, 
#'  Cambridge University Press.
#' @seealso \pkg{\link{ExtDist}} for other standard distributions.
#' @author Haizhen Wu and A. Jonathan R. Godfrey \cr
#' Updates and bug fixes by Sarah Pirikahu.
#'
#' @examples
#' # Parameter estimation for a distribution with known shape parameters
#' X <- rBeta_ab(n=500, shape1=2, shape2=5, a=1, b=2)
#' est.par <- eBeta_ab(X); est.par
#' plot(est.par)
#' 
#' # Fitted density curve and histogram
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dBeta_ab(den.x,params = est.par)
#' hist(X, breaks=10, probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue")   # Original data
#' lines(density(X), lty=2)          # Fitted density curve
#' 
#' # Extracting boundary and shape parameters
#' est.par[attributes(est.par)$par.type=="boundary"]
#' est.par[attributes(est.par)$par.type=="shape"]
#' 
#' # Parameter Estimation for a distribution with unknown shape parameters
#' # Example from: Bury(1999) pp.261-262, parameter estimates as given by Bury are 
#' # shape1 = 4.088, shape2 = 10.417, a = 1.279 and b = 2.407. 
#' # The log-likelihood for this data and Bury's parameter estimates is 8.598672.
#' data <- c(1.73, 1.5, 1.56, 1.89, 1.54, 1.68, 1.39, 1.64, 1.49, 1.43, 1.68, 1.61, 1.62)
#' est.par <- eBeta_ab(X=data, method="numerical.MLE");est.par
#' plot(est.par)
#' 
#' # Estimates calculated by eBeta_ab differ from those given by Bury(1999).
#' # However, eBeta_ab's parameter estimates appear to be an improvement, due to a larger 
#' # log-likelihood of 9.295922 (as given by lBeta_ab below).
#' 
#' # log-likelihood and score functions 
#' lBeta_ab(data,param = est.par)
#' sBeta_ab(data,param = est.par)

#' @rdname Beta_ab
#' @export dBeta_ab
dBeta_ab <-function(x, shape1=2, shape2=3, a = 0, b=1, params = list(shape1, shape2, a, b),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
      a <- params$a
      b <- params$b
    }
    out <- (x>=a & x<=b) * dbeta((x-a)/(b-a),shape1,shape2)/(b-a)
    return(out)
  }

#' @rdname Beta_ab
#' @export pBeta_ab
pBeta_ab <- function(q, shape1=2, shape2=3, a = 0, b=1, params = list(shape1=2, shape2 = 5, a = 0, b = 1),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
      a <- params$a
      b <- params$b
    }
    out <- pbeta((q-a)/(b-a),shape1,shape2)
    return(out)
  }

#' @rdname Beta_ab
#' @export qBeta_ab
qBeta_ab <- function(p, shape1=2, shape2=3, a = 0, b=1, params = list(shape1=2, shape2 = 5, a = 0, b = 1),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
      a <- params$a
      b <- params$b
    }
    out <- (b-a)*qbeta(p,shape1,shape2) + a
    return(out)
  }

#' @rdname Beta_ab
#' @export rBeta_ab
rBeta_ab <- function(n, shape1=2, shape2=3, a = 0, b = 1, params = list(shape1, shape2, a, b),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
      a <- params$a
      b <- params$b
    }
    X <- rbeta(n,shape1,shape2)
    out <- (b-a)*X + a
    return(out)
  }

#' @rdname Beta_ab
#' @export eBeta_ab
eBeta_ab <- function(X,w, method ="numerical.MLE",...){
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
  est.par <- wmle(X=X, w=w, distname = "Beta_ab",
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
attributes(est.par)$distname <- "Beta_ab"
attributes(est.par)$method <- method
attributes(est.par)$par.name <- c("shape1","shape2","a","b")
attributes(est.par)$par.type <- c("shape","shape","boundary","boundary")
attributes(est.par)$par.vals <- c(est.par$shape1, est.par$shape2, est.par$a, est.par$b)
attributes(est.par)$par.s.e <-  est.par.se  

class(est.par) <- "eDist"

return(est.par)
  }

#' @rdname Beta_ab
#' @export lBeta_ab
## (weighted) (log) likelihood function
lBeta_ab <- function(X, w, shape1=2, shape2 =3, a = 0, b = 1,  params = list(shape1, shape2, a, b), logL = TRUE,...){
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
    
    #     ll <- sum(w*((shape1-1)*log(X-a)+(shape2-1)*log(b-X)-log(beta(shape1,shape2))-(shape1+shape2-1)*log(b-a)))
    ll <- sum(w*log(dBeta_ab(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }

#' @rdname Beta_ab
#' @export sBeta_ab
## (weighted) score vectors
sBeta_ab <- function(X, w, shape1=2, shape2 =3, a = 0, b = 1,  params = list(shape1, shape2, a, b),...){
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
    
    score1 <- sum(w*(digamma(shape1+shape2)-digamma(shape1)+log(X-a)-log(b-a)))
    score2 <- sum(w*(digamma(shape1+shape2)-digamma(shape2)+log(b-X)-log(b-a)))
    score3 <- sum(w*((shape1+shape2-1)/(b-a)-(shape1-1)/(X-a)))
    score4 <- sum(w*((shape2-1)/(b-X)-(shape1+shape2-1)/(b-a)))
    
    score <- c(score1,score2,score3,score4)
    names(score) <- c("shape1","shape2","a","b")
    return(score)
  }



