#' @title The truncated normal distribution.
#' @description Density, distribution, quantile, random number 
#' generation and parameter estimation functions for the truncated normal distribution with parameters \code{mean}, \code{sd} and 
#' \code{a} and \code{b} which represent the lower and upper truncation points respectively. 
#' Parameter estimation can be based on a weighted or unweighted i.i.d. sample and is performed numerically. 
#' 
#' @rdname Normal_trunc_ab
#' @name Normal_trunc_ab
#' @aliases dNormal_trunc_ab
#' @aliases pNormal_trunc_ab 
#' @aliases qNormal_trunc_ab 
#' @aliases rNormal_trunc_ab
#' @aliases eNormal_trunc_ab 
#' @aliases lNormal_trunc_ab
#' 
#' @details If the \code{mean}, \code{sd}, \code{a} or \code{b} are not specified they assume the default values of 0, 1, 0, 1 respectively.\cr
#' \cr
#' The \code{dNormal_trunc_ab()}, \code{pNormal_trunc_ab()}, \code{qNormal_trunc_ab()},and \code{rNormal_trunc_ab()} functions serve
#'  as wrappers of the \code{\link[truncdist]{dtrunc}}, \code{\link[truncdist]{ptrunc}}, \code{\link[truncdist]{qtrunc}}, and
#' \code{\link[truncdist]{rtrunc}} functions in the \pkg{{truncdist}} package. They allow for the parameters to be declared not only as 
#' individual numerical values, but also as a list so parameter estimation can be carried out. \cr
#' \cr
#' The probability density function of the doubly truncated normal distribution is given by 
#' \deqn{f(x) = \sigma^{-1} Z(x-\mu/\sigma)[\Phi(b-\mu/\sigma) - \Phi(a-\mu/\sigma)]^{-1}}
#' where \eqn{\infty <a \le x \le b < \infty}. The degrees of truncation are \eqn{\Phi((a-\mu)/\sigma)} from below and \eqn{1-\Phi((a-\mu)/\sigma)}
#' from above. If a is replaced by \eqn{-\infty}, or b by \eqn{\infty}, the distribution is singly truncated, (Johnson et.al, p.156).
#' The upper and lower limits of truncation \eqn{a} and \eqn{b} are normally known
#' parameters whereas \eqn{\mu} and \eqn{\sigma} may be unknown. Crain (1979) discusses parameter estimation for the truncated normal 
#' distribution and the method of numerical maximum likelihood estimation is used for parameter estimation in \code{eNormal_trunc_ab}.\cr 
#' \cr  
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
#' @param logL logical;if TRUE, lNormal_trunc_ab gives the log-likelihood, otherwise the likelihood is given.
#' @param ... Additional parameters.
#'
#' @return dNormal_trunc_ab gives the density, pNormal_trunc_ab the distribution function,
#' qNormal_trunc_ab the quantile function, rNormal_trunc_ab generates random variables, and 
#' eNormal_trunc_ab estimates the parameters. lNormal_trunc_ab provides the log-likelihood function. 
#'
#' @references Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994) Continuous Univariate Distributions,
#'  volume 1, chapter 13, Wiley, New York.\cr
#' \cr
#' Crain, B.R (1979). Estimating the parameters of a truncated normal distribution, Applied Mathematics and Computation,
#' vol 4, pp. 149-156
#' @seealso \pkg{\link{ExtDist}} for other standard distributions.
#' @author Haizhen Wu and A. Jonathan R. Godfrey.\cr
#' Updates and bug fixes by Sarah Pirikahu.

#' @examples
#' # Parameter estimation for a distribution with known shape parameters
#' X <- rNormal_trunc_ab(n= 500, mu= 2, sigma = 5, a = 1, b = 2)
#' est.par <- eNormal_trunc_ab(X); est.par
#' plot(est.par)
#' 
#' #  Fitted density curve and histogram
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dNormal_trunc_ab(den.x,params = est.par)
#' hist(X, breaks=10, probability=TRUE, ylim = c(0,1.2*max(den.y)))
#' lines(den.x, den.y, col="blue")
#' lines(density(X), lty = 2)
#' 
#' # Extracting boundary and shape parameters
#' est.par[attributes(est.par)$par.type=="boundary"]
#' est.par[attributes(est.par)$par.type=="shape"]
#' 
#' # log-likelihood function
#' lNormal_trunc_ab(X,param = est.par)


#' @rdname Normal_trunc_ab
#' @export dNormal_trunc_ab
dNormal_trunc_ab <-function(x, mu=0, sigma=1, a = 0, b=1, params = list(mu, sigma, a, b),...){
    if(!missing(params)){
      mu <- params$mu; sigma <- params$sigma; a <- params$a; b <- params$b
    }
    out <- truncdist::dtrunc( x, spec="norm", a=a, b=b, mean = mu, sd = sigma)
    return(out)
  }

#' @rdname Normal_trunc_ab
#' @export pNormal_trunc_ab
pNormal_trunc_ab <- function(q, mu=0, sigma=1, a = 0, b=1, params = list(mu=2, sigma = 5, a = 0, b = 1),...){
    if(!missing(params)){
      mu <- params$mu; sigma <- params$sigma; a <- params$a; b <- params$b
    }
    out <- truncdist::ptrunc( q, spec="norm", a=a, b=b, mean = mu, sd = sigma)
    return(out)
  }

#' @rdname Normal_trunc_ab
#' @export qNormal_trunc_ab
qNormal_trunc_ab <- function(p, mu=0, sigma=1, a = 0, b=1, params = list(mu=2, sigma = 5, a = 0, b = 1),...){
    if(!missing(params)){
      mu <- params$mu; sigma <- params$sigma; a <- params$a; b <- params$b
    }
    out <- truncdist::qtrunc( p, spec="norm", a=a, b=b, mean = mu, sd = sigma)
    return(out)
  }

#' @rdname Normal_trunc_ab
#' @export rNormal_trunc_ab
rNormal_trunc_ab <- function(n, mu=0, sigma=1, a = 0, b = 1, params = list(mu, sigma, a, b),...){
    if(!missing(params)){
      mu <- params$mu; sigma <- params$sigma; a <- params$a; b <- params$b
    }
    out <- truncdist::rtrunc(n, spec="norm", a=a, b=b, mean = mu, sd = sigma)
    return(out)
  }

#' @rdname Normal_trunc_ab
#' @export eNormal_trunc_ab
eNormal_trunc_ab <- function(X,w, method ="numerical.MLE",...){
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
  est.par <- wmle(X=X, w=w, distname = "Normal_trunc_ab",
                  initial=list(mu=mean(min(X),max(X)),sigma=1,a=min(X)-0.1*d,b=max(X)+0.1*d),
                  lower=list(mu=-Inf,sigma=0,a=-Inf,b=max(X)),
                  upper=list(mu=Inf,sigma=Inf,a=min(X),b=Inf)
				  )
  
  est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
  if(class(est.par.se) == "try-error") {
    est.par.se <- rep(NA, length(est.par))
  } 
} 

attributes(est.par)$ob <- X
attributes(est.par)$weights <- w
attributes(est.par)$distname <- "Normal_trunc_ab"
attributes(est.par)$method <- method
attributes(est.par)$par.name <- c("mu","sigma","a","b")
attributes(est.par)$par.type <- c("shape","shape","boundary","boundary")
attributes(est.par)$par.vals <- c(est.par$mu, est.par$sigma, est.par$a, est.par$b)
attributes(est.par)$par.s.e <-  est.par.se  

class(est.par) <- "eDist"

return(est.par)
  }

#' @rdname Normal_trunc_ab
#' @export lNormal_trunc_ab
## (weighted) (log) likelihood function
lNormal_trunc_ab <- function(X, w, mu=0, sigma =1, a = 0, b = 1,  params = list(mu, sigma, a, b), logL = TRUE,...){
    if(!missing(params)){
      mu <- params$mu; sigma <- params$sigma; a <- params$a; b <- params$b
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dNormal_trunc_ab(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }