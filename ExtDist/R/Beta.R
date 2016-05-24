#' @title The Standard Beta Distribution.
#' @description Density, distribution, quantile, random number 
#' generation, and parameter estimation functions for the beta distribution with parameters \code{shape1} and \code{shape2}.
#' Parameter estimation can be based on a weighted or unweighted i.i.d. sample and can be carried out analytically or numerically. 
#' @rdname Beta
#' @name Beta
#' @aliases dBeta 
#' @aliases pBeta 
#' @aliases qBeta 
#' @aliases rBeta 
#' @aliases eBeta 
#' @aliases lBeta 
#' @aliases sBeta 
#' @aliases iBeta
#' 
#' @param params A list that includes all named parameters.
#' @param x,q Vector of quantiles.
#' @param w Optional vector of sample weights.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param X Sample observations.
#' @param shape1,shape2 Shape parameters.
#' @param method Parameter estimation method.
#' @param logL logical, if TRUE lBeta gives the log-likelihood, otherwise the likelihood is given.  
#' @param ... Additional parameters.
#' @details The \code{dBeta()}, \code{pBeta()}, \code{qBeta()},and \code{rBeta()} functions serve as wrappers of the standard
#' \code{\link[stats]{dbeta}}, \code{\link[stats]{pbeta}}, \code{\link[stats]{qbeta}}, and \code{\link[stats]{rbeta}} functions 
#' in the \pkg{\link{stats}} package. They allow for the shape parameters to be declared not only as 
#' individual numerical values, but also as a list so parameter estimation can be carried out. \cr
#' \cr
#' The beta distribution with parameters \code{shape1}=\eqn{\alpha} and \code{shape2}=\eqn{\beta} is given by   
#' \deqn{f(x) = \frac{x^{\alpha-1} (1-x)^{\beta-1}}{B(\alpha,\beta)}}
#' where \eqn{0 \le x \le 1}, \eqn{\alpha>0}, \eqn{\beta>0}, and \eqn{B} is the \link[base]{beta} function. \cr
#' \cr
#' Analytical parameter estimation is conducted using the method of moments. The parameter
#' estimates for \eqn{\alpha} and \eqn{\beta} are as given in the \href{http://www.itl.nist.gov/div898/handbook/eda/section3/eda366h.htm}{Engineering Statistics Handbook}.\cr
#' \cr
#' The log-likelihood function of the beta distribution is given by 
#' \deqn{l(\alpha, \beta | x) = (\alpha-1)\sum_{i} ln(x_i) + (\beta-1)\sum_{i}
#' ln(1-x_i) - ln B(\alpha,\beta).}
#' Aryal & Nadarajah (2004) derived the score function and Fisher's information matrix for the 4-parameter beta function, from
#' which the 2-parameter cases can be obtained.
#'
#' @return dBeta gives the density, pBeta the distribution function,
#' qBeta the quantile function, rBeta generates random deviates, and 
#' eBeta estimates the parameters. lBeta provides the log-likelihood function, sBeta the observed score function,
#' and iBeta the observed information matrix. 
#'
#' @author Haizhen Wu and A. Jonathan R. Godfrey.\cr
#' Updates and bug fixes by Sarah Pirikahu.
#'
#' @references Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions,
#'  volume 2, chapter 25, Wiley, New York.\cr
#'  \cr
#'  \href{http://www.itl.nist.gov/div898/handbook/eda/section3/eda366h.htm}{Engineering Statistics Handbook}\cr
#'  \cr
#'  Bury, K. (1999) Statistical Distributions in Engineering, Chapter 14, pp.253-255, 
#'  Cambridge University Press.\cr
#'  \cr
#'  Aryal, G. and Nadarajah, S. (2004) Information Matrix for Beta Distributions, Serdica Math. J. 30, 513-526.
#'  
#' @seealso \pkg{\link{ExtDist}} for other standard distributions. 
#' 
#' @examples
#' # Parameter estimation for a distribution with known shape parameters
#' x <- rBeta(n=500, params=list(shape1=2, shape2=2))
#' est.par <- eBeta(x); est.par
#' plot(est.par)
#' 
#' # Fitted density curve and histogram
#' dens <- dBeta(x=seq(0,1,length=100), params=list(shape1=2, shape2=2)) 
#' hist(x, breaks=10, probability=TRUE, ylim = c(0,1.2*max(dens)))
#' lines(seq(0,1,length=100), dens, col="blue")
#' lines(density(x), lty=2)
#' 
#' # Extracting shape parameters
#' est.par[attributes(est.par)$par.type=="shape"]
#' 
#' # Parameter estimation for a distribution with unknown shape parameters
#' # Example from; Bury(1999) pp.253-255, parameter estimates as given by Bury are
#' # shape1 = 4.222 and shape2 = 6.317
#' data <- c(0.461, 0.432, 0.237, 0.113, 0.526, 0.278, 0.275, 0.309, 0.67, 0.428, 0.556,
#' 0.402, 0.472, 0.226, 0.632, 0.533, 0.309, 0.417, 0.495, 0.241)
#' est.par <- eBeta(X=data, method="numerical.MLE"); est.par
#' plot(est.par)
#' 
#' # Log-likelihood, score function, and observed information matrix  
#' lBeta(data, param=est.par)
#' sBeta(data, param=est.par)
#' iBeta(data, param=est.par)
#' 
#' # Evaluating the precision of parameter estimation by the Hessian matrix.
#' H <- attributes(est.par)$nll.hessian;H
#' var <- solve(H)  
#' se <- sqrt(diag(var)); se

 
#' @rdname Beta
#' @export dBeta
dBeta <-function(x, shape1=2, shape2 =3, params = list(shape1, shape2),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    out <- dbeta(x, shape1, shape2)
    return(out)
  }

#' @rdname Beta
#' @export pBeta
pBeta <- function(q, shape1=2, shape2 =3, params = list(shape1, shape2),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    out <- pbeta(q, shape1, shape2)
    return(out)
  }

#' @rdname Beta
#' @export qBeta
qBeta <- function(p, shape1=2, shape2 =3, params = list(shape1, shape2),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    qbeta(p, shape1 = shape1, shape2 = shape2)
  }

#' @rdname Beta
#' @export rBeta
rBeta <- function(n, shape1=2, shape2 =3, params = list(shape1, shape2),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    rbeta(n, shape1 = shape1, shape2 = shape2)
  }

#' @rdname Beta
#' @export eBeta
eBeta <-function(X,w, method =c("MOM","numerical.MLE"),...){
  method <- match.arg(method)
  # Accounting for weights
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    # Method of moments parameter estimation
    if(method == "MOM") {
      sample.mean <- mean(w*X)
      sample.var <- (mean(w*X^2)-sample.mean^2)*n/(n-1)
      
      v <- sample.mean*(1-sample.mean)
      if(sample.var<v){
        shape1 <- sample.mean*(v/sample.var-1)
        shape2 <- (1-sample.mean)*(v/sample.var-1)
      } else {
        shape2 <- sample.mean*(v/sample.var-1)
        shape1 <- (1-sample.mean)*(v/sample.var-1)      
      } 
      # Need to find analytic solns for standard errors of parameter estimates. 
      est.par <- list(shape1 = shape1, shape2 = shape2)
      est.par.se <- rep(NA, length(est.par))
    } else {  
      # Numerical parameter estimation
      method = "numerical.MLE"
      est.par <- wmle(X=X, w=w, distname = "Beta",
                      initial=list(shape1 = 3, shape2 = 3),
                      lower=list(shape1 = 1, shape2 = 1),
                      upper=list(shape1 = Inf, shape2 = Inf))
      
      est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
      if(class(est.par.se) == "try-error") {
        est.par.se <- rep(NA, length(est.par))
      } 
    }
    
    attributes(est.par)$ob <- X
    attributes(est.par)$weights <- w
    attributes(est.par)$distname <- "Beta"
    attributes(est.par)$method <- method
    attributes(est.par)$par.name <- c("shape1","shape2")
    attributes(est.par)$par.type <- c("shape","shape")
    attributes(est.par)$par.vals <- c(est.par$shape1, est.par$shape2)
    attributes(est.par)$par.s.e <-  est.par.se  
    
    class(est.par) <- "eDist"
    
    return(est.par)
  }

#' @rdname Beta
#' @export lBeta
## (weighted) (log) likelihood function
lBeta <- function(X, w, shape1=2, shape2 =3, params = list(shape1, shape2), logL = TRUE,...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    # Accounting for weights
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*((shape1-1)*log(X)+(shape2-1)*log(1-X)-log(beta(shape1,shape2))))
    #     ll <- sum(w*log(dBeta(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }

#' @rdname Beta
#' @export sBeta
## (weighted) score vectors
sBeta <- function(X, w, shape1=2, shape2 =3, params = list(shape1, shape2),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    # Accounting for weights
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    score1 <- sum(w*(digamma(shape1+shape2)-digamma(shape1)+log(X)))
    score2 <- sum(w*(digamma(shape1+shape2)-digamma(shape2)+log(1-X)))
    
    score <- c(score1,score2)
    names(score) <- c("shape1","shape2")
    return(score)
  }

#' @rdname Beta
#' @export iBeta
## (weighted) (observed) information matrix
iBeta <- function(X, w, shape1=2, shape2 =3, params = list(shape1, shape2),...){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    # Accounting for weights
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    info11 <- -sum(w*(trigamma(shape1+shape2)-trigamma(shape1)))
    info12 <- -sum(w*trigamma(shape1+shape2))
    info21 <- info12
    info22 <- -sum(w*(trigamma(shape1+shape2)-trigamma(shape2)))
    info <- matrix(c(info11,info12,info21,info22), nrow=2,ncol=2)
    rownames(info) <- colnames(info) <- c("shape1","shape2")
    return(info)
  }
