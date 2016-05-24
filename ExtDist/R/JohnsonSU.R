#' @title The Johnson SU distribution.
#' @description Density, distribution, quantile, random number 
#' generation and parameter estimation functions for the Johnson SU (unbounded support) distribution.
#' Parameter estimation can be based on a weighted or unweighted i.i.d sample and can be carried out numerically.
#' 
#' @rdname JohnsonSU
#' @name JohnsonSU
#' @aliases dJohnsonSU
#' @aliases pJohnsonSU
#' @aliases qJohnsonSU
#' @aliases rJohnsonSU
#' @aliases eJohnsonSU
#' @aliases lJohnsonSU
#' 
#' @details The Johnson system of distributions consists of families of distributions that, through specified transformations, can be 
#' reduced to the standard normal random variable. It provides a very flexible system for describing statistical distributions
#' and is defined by 
#' \deqn{z = \gamma + \delta f(Y)} 
#' with \eqn{Y = (X-xi)/lambda}. The Johnson SB distribution arises when \eqn{f(Y) = archsinh(Y)}, where \eqn{-\infty < Y < \infty}. 
#' This is the unbounded Johnson family since the range of Y is \eqn{(-\infty,\infty)}, Karian & Dudewicz (2011).\cr
#' \cr
#' The JohnsonSU distribution has probability density function 
#' \deqn{p_X(x) = \frac{\delta}{\sqrt{2\pi((x-xi)^2 + lambda^2)}}exp[-0.5(\gamma + \delta ln(\frac{x-xi + \sqrt{(x-xi)^2 + lambda^2}}{lambda}))^2].}
#' Parameter estimation can only be carried out numerically.
#' @param params A list that includes all named parameters.
#' @param x,q A vector of quantiles.
#' @param w An optional vector of sample weights.
#' @param p A vector of probabilities.
#' @param n Number of observations.
#' @param X Sample observations.
#' @param xi,lambda Location-scale parameters.
#' @param gamma,delta Shape parameters.
#' @param method Parameter estimation method.
#' @param logL logical; if TRUE, lJohnsonSU gives the log-likelihood, otherwise the likelihood is given.
#' @param ... Additional parameters.
#'
#' @return dJohnsonSU gives the density, pJohnsonSU the distribution function,
#' qJohnsonSU gives the quantile function, rJohnsonSU generates random variables, and  
#' eJohnsonSU estimates the parameters. lJohnsonSU provides the log-likelihood function.
#' @seealso \pkg{\link{ExtDist}} for other standard distributions.
#' @author Haizhen Wu and A. Jonathan R. Godfrey.\cr
#' Updates and bug fixes by Sarah Pirikahu.
#'
#' @references Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994) Continuous Univariate Distributions,
#'  volume 1, chapter 12, Wiley, New York.\cr
#'  \cr
#'  Bowman, K.O., Shenton, L.R. (1983). Johnson's system of distributions. In: Encyclopedia
#'  of Statistical Sciences, Volume 4, S. Kotz and N.L. Johnson (eds.), pp. 303-314. John
#'  Wiley and Sons, New York.\cr
#'  \cr
#'   Z. A. Karian and E. J. Dudewicz (2011) Handbook of Fitting Statistical Distributions with R, Chapman & Hall.\cr 
#'  
#' @examples
#' # Parameter estimation for a known distribution
#' X <- rJohnsonSU(n=500, gamma=-0.5, delta=2, xi=-0.5, lambda=2)
#' est.par <- eJohnsonSU(X); est.par
#' plot(est.par)
#' 
#' # Fitted density curve and histogram
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dJohnsonSU(den.x,params = est.par)
#' hist(X, breaks=10, probability=TRUE, ylim = c(0,1.2*max(den.y)))
#' lines(den.x, den.y, col="blue")
#' lines(density(X), lty=2) 
#' 
#' # Extracting shape and boundary parameters
#' est.par[attributes(est.par)$par.type=="shape"]
#' est.par[attributes(est.par)$par.type=="boundary"]
#' 
#' # Parameter Estimation for a distribution with unknown shape parameters
#' # Example from Karian, Z.A and Dudewicz, E.J. (2011) p.657.
#' # Parameter estimates as given by Karian & Dudewicz are:
#' # gamma =-0.2823, delta=1.0592, xi = -1.4475 and lambda = 4.2592  with log-likelihood = -277.1543
#' data <- c(1.99, -0.424, 5.61, -3.13, -2.24, -0.14, -3.32, -0.837, -1.98, -0.120,
#'          7.81, -3.13, 1.20, 1.54, -0.594, 1.05, 0.192, -3.83, -0.522, 0.605,
#'          0.427, 0.276, 0.784, -1.30, 0.542, -0.159, -1.66, -2.46, -1.81, -0.412,
#'          -9.67, 6.61, -0.589, -3.42, 0.036, 0.851, -1.34, -1.22, -1.47, -0.592,
#'          -0.311, 3.85, -4.92, -0.112, 4.22, 1.89, -0.382, 1.20, 3.21, -0.648,
#'          -0.523, -0.882, 0.306, -0.882, -0.635, 13.2, 0.463, -2.60, 0.281, 1.00,
#'          -0.336, -1.69, -0.484, -1.68, -0.131, -0.166, -0.266, 0.511, -0.198, 1.55,
#'          -1.03, 2.15, 0.495, 6.37, -0.714, -1.35, -1.55, -4.79, 4.36, -1.53, 
#'          -1.51, -0.140, -1.10, -1.87, 0.095, 48.4, -0.998, -4.05, -37.9, -0.368,
#'          5.25, 1.09, 0.274, 0.684, -0.105, 20.3, 0.311, 0.621, 3.28, 1.56)
#' est.par <- eJohnsonSU(data); est.par
#' plot(est.par)
#'
#' # Estimates calculated by eJohnsonSU differ from those given by Karian & Dudewicz (2011).
#' # However, eJohnsonSU's parameter estimates appear to be an improvement, due to a larger 
#' # log-likelihood of -250.3208 (as given by lJohnsonSU below).
#'
#' # log-likelihood function
#' lJohnsonSU(data, param = est.par)
#'
#' # Evaluation of the precision using the Hessian matrix
#' H <- attributes(est.par)$nll.hessian
#' var <- solve(H)
#' se <- sqrt(diag(var)); se

#' @rdname JohnsonSU
#' @export dJohnsonSU
dJohnsonSU <-function(x, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2, params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2),...){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }

    u = (x - xi)/lambda
    out = x
    out[abs(u)<1e-3] = dnorm(gamma+delta*u[abs(u)<1e-3])*delta/lambda  # good approximation when abs(u) is small. log(u+sqrt(1+u^2)) ~ u.
    out[u < (-1e10)] = dnorm(gamma+delta*(-Inf)) *(delta/lambda/sqrt(u[u<(-1e10)]^2+1))  # for u<-1e10, log(u+sqrt(1+u^2)) ~ -Inf . but if u = -1e200, R gives Inf
    out[abs(u)>=1e-3 & u>(-1e10)] = dnorm(gamma+delta*log(u[abs(u)>=1e-3 & u>(-1e10)]+sqrt(1+u[abs(u)>=1e-3 & u>(-1e10)]^2))) *(delta/lambda/sqrt(u[abs(u)>=1e-3 & u>(-1e10)]^2+1))    
    return(out)
  }

#' @rdname JohnsonSU
#' @export pJohnsonSU
pJohnsonSU <- function(q, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2, params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2),...){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }

    u = (q - xi)/lambda
    out = pnorm(delta*asinh(u)+gamma)
    return(out)
  }

#' @rdname JohnsonSU
#' @export qJohnsonSU
qJohnsonSU <- function(p, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2, params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2),...){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }
    u = (qnorm(p)-gamma)/delta
    out = lambda*sinh((qnorm(p)-gamma)/delta)+xi
    return(out)
  }

#' @rdname JohnsonSU
#' @export rJohnsonSU
rJohnsonSU <- function(n, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2, params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2),...){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }
    out <- SuppDists::rJohnson(n, parms = list(gamma=gamma, delta=delta, xi=xi, lambda=lambda, type= "SU"))
    return(out)
  }

#' @rdname JohnsonSU
#' @export eJohnsonSU
eJohnsonSU <- function(X,w, method ="numerical.MLE",...){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
{if(method != "numerical.MLE") warning(paste("method ", method, " is not avaliable, use numerial.MLE instead."))  
 method = "numerical.MLE"  
	
 est.par <- wmle(X=X, w=w, distname = "JohnsonSU",
          initial=list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2),
          lower=list(gamma = -Inf, delta = 1, xi = -Inf, lambda = 0),
          upper=list(gamma = Inf, delta = Inf, xi = Inf, lambda = Inf)
 )

 est.par.se <- try(suppressWarnings(sqrt(diag(solve(attributes(est.par)$nll.hessian )))),silent=TRUE)
 if(class(est.par.se) == "try-error" | any(is.nan(est.par.se))) {
   est.par.se <- rep(NA, length(est.par))
 } 
}

attributes(est.par)$ob <- X
attributes(est.par)$weights <- w
attributes(est.par)$distname <- "JohnsonSU"
attributes(est.par)$method <- method
attributes(est.par)$par.name <- c("gamma","delta","xi","lambda")
attributes(est.par)$par.type <- c("shape","shape","boundary","boundary")
attributes(est.par)$par.vals <- c(est.par$gamma, est.par$delta, est.par$xi, est.par$lambda)
attributes(est.par)$par.s.e <-  est.par.se  

class(est.par) <- "eDist"

return(est.par)
  }

#' @rdname JohnsonSU
#' @export lJohnsonSU
## (weighted) (log) likelihood function
lJohnsonSU <- function(X, w, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2,  params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2), logL = TRUE,...){
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
    
    ll <- sum(w*log(dJohnsonSU(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }
