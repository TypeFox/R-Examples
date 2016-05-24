#' @title The Laplace Distribution.
#' @description Density, distribution, quantile, random number generation 
#' and parameter estimation functions for the Laplace distribution with \code{location} parameter \eqn{\mu} and 
#' \code{scale} parameter \eqn{b}. Parameter estimation can for the Laplace distribution can be carried out numerically
#' or analytically but may only be based on an unweighted i.i.d. sample.
#'
#' @rdname Laplace
#' @name Laplace
#' @aliases dLaplace 
#' @aliases pLaplace
#' @aliases qLaplace 
#' @aliases rLaplace 
#' @aliases eLaplace 
#' @aliases lLaplace
#' @aliases sLaplace 
#' @aliases iLaplace
#'
#' @details The \code{dLaplace()}, \code{pLaplace()}, \code{qLaplace()},and \code{rLaplace()} functions allow for the parameters to be declared not only as 
#' individual numerical values, but also as a list so parameter estimation can be carried out. \cr
#' \cr
#' The Laplace distribution with parameters \code{location} = \eqn{\mu} and \code{scale}=\eqn{b} has probability density 
#' function 
#' \deqn{f(x) = (1/2b) exp(-|x-\mu|/b)}
#' where \eqn{-\infty < x < \infty} and \eqn{b > 0}. The cumulative distribution function for \code{pLaplace} is defined 
#' by Johnson et.al (p.166). \cr
#' \cr
#' Parameter estimation can be carried out analytically via maximum likelihood estimation, see Johnson et.al (p.172). Where the population 
#' mean, \eqn{\mu}, is estimated using the sample median and \eqn{b} by the mean of \eqn{|x-b|}.\cr
#' \cr 
#' Johnson et.al (p.172) also provides the log-likelihood function for the Laplace distribution 
#' \deqn{l(\mu, b | x) = -n ln(2b) - b^{-1} \sum |xi-\mu|.} 
#' 
#'  
#' @param params A list that includes all named parameters
#' @param x,q A vector of quantiles.
#' @param w Optional vector of sample weights.
#' @param p A vector of probabilities.
#' @param n Number of observations.
#' @param X Sample observations.
#' @param mu Location parameter.
#' @param b Scale parameter.
#' @param method Parameter estimation method.
#' @param logL logical; if TRUE, lLaplace gives the log-likelihood, otherwise the likelihood is given.
#' @param ... Additional parameters. 
#'
#' @return dLaplace gives the density, pLaplace the distribution function,
#' qLaplace the quantile function, rLaplace generates random deviates, and 
#' eLaplace estimates the distribution parameters. lLaplace provides the log-likelihood function, sLaplace the score function, 
#' and iLaplace the observed information matrix. 
#'
#' @note The estimation of the population mean is done using the median of the sample. Unweighted 
#' samples are not yet catered for in the eLaplace() function.
#'
#' @author A. Jonathan R. Godfrey and Haizhen Wu. \cr
#' Updates and bug fixes by Sarah Pirikahu
#' 
#' @references  Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions,
#'  volume 2, chapter 24, Wiley, New York.\cr
#'  \cr
#'  Best, D.J., Rayner, J.C.W. and Thas O. (2008) Comparison of some tests of fit for the Laplace distribution, 
#'  Computational Statistics and Data Analysis, Vol. 52, pp.5338-5343.\cr 
#'  \cr
#'  Gumbel, E.J., Mustafi, C.K., 1967. Some analytical properties of bivariate extremal distributions.
#'  J. Am. Stat. Assoc. 62, 569-588
#'  
#' @seealso \pkg{\link{ExtDist}} for other standard distributions.
#' 
#' @examples
#' # Parameter estimation for a distribution with known shape parameters
#' X <- rLaplace(n=500, mu=1, b=2)
#' est.par <- eLaplace(X, method="analytic.MLE"); est.par
#' plot(est.par)
#' 
#' #  Fitted density curve and histogram
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dLaplace(den.x, location = est.par$location, scale= est.par$scale)
#' hist(X, breaks=10, probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue")
#' lines(density(X), lty=2)
#' 
#' # Extracting location or scale parameters
#' est.par[attributes(est.par)$par.type=="location"]
#' est.par[attributes(est.par)$par.type=="scale"]
#' 
#' # Parameter estimation for a distribution with unknown shape parameters
#' # Example from Best et.al (2008). Original source of flood data from Gumbel & Mustafi.
#' # Parameter estimates as given by Best et.al mu=10.13 and  b=3.36
#' flood <- c(1.96, 1.96, 3.60, 3.80, 4.79, 5.66, 5.76, 5.78, 6.27, 6.30, 6.76, 7.65, 7.84, 7.99,
#'            8.51, 9.18, 10.13, 10.24, 10.25, 10.43, 11.45, 11.48, 11.75, 11.81, 12.34, 12.78, 13.06, 
#'            13.29, 13.98, 14.18, 14.40, 16.22, 17.06)
#' est.par <- eLaplace(flood, method="numerical.MLE"); est.par
#' plot(est.par) 
#' 
#' #log-likelihood function
#' lLaplace(flood,param=est.par)
#' 
#' # Evaluating the precision by the Hessian matrix
#' H <- attributes(est.par)$nll.hessian
#' var <- solve(H)
#' se <- sqrt(diag(var));se


#' @rdname Laplace
#' @export dLaplace
dLaplace <- function(x, mu=0, b=1, params=list(mu, b),...){
  if(!missing(params)){
    mu <- params$mu
    b <- params$b
  }
  # pdf Laplace dist Johnson et.al, Vol 2, chapter 24, p.164.
  d <- exp(-abs(x-mu)/b) / (2*b)
}

#' @rdname Laplace
#' @export pLaplace
pLaplace <- function(q, mu=0, b=1,params=list(mu, b),...){
  if(!missing(params)){
    mu <- params$mu
    b <- params$b
  }
  # simplification of the Laplace dist cdf Johnson et.al, Vol 2, chapter 24, p.166
  x <- q-mu
  0.5 + 0.5*sign(x)*(1-exp(-abs(x)/b))
}

#' @rdname Laplace
#' @export qLaplace
qLaplace <- function(p, mu=0, b=1,params=list(mu, b),...){
  if(!missing(params)){
    mu <- params$mu
    b <- params$b
  }
  
  x <- p-0.5
  mu-b*sign(x)*log(1-2*abs(x))
}

#' @rdname Laplace
#' @export rLaplace
rLaplace <- function(n, mu=0, b=1,params=list(mu, b),...){
  if(!missing(params)){
    mu <- params$mu
    b <- params$b}
  
  u<-runif(n)-0.5
  x<-mu-b*sign(u)*log(1-2*abs(u))
}

#' @rdname Laplace
#' @export eLaplace
eLaplace <-function(X,w, method =c("analytic.MLE","numerical.MLE"),...){
    n<-length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    if(method == "analytic.MLE") {
      
      mu <- median(X, na.rm=T)
      b <- mean(abs(X-mu))
      est.par <- list(mu=mu, b=b) 
      est.par.se <- rep(NA, length(est.par))
    } else{
      method == "numerical.MLE"  
      
      est.par <- wmle(X=X, w=w, distname="Laplace",
                      initial=list(mu=0, b=1),
                      lower=list(mu=-Inf, b=0),
                      upper=list(mu=Inf, b=Inf))
      
      est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
      if(class(est.par.se) == "try-error") {
        est.par.se <- rep(NA, length(est.par))
      } 
    } 
    
    attributes(est.par)$ob <- X
    attributes(est.par)$weights <- w
    attributes(est.par)$distname <- "Laplace"
    
    attributes(est.par)$method <- method
    attributes(est.par)$par.name <- c("mu","b")
    attributes(est.par)$par.type <- c("location", "scale")
    attributes(est.par)$par.vals <- c(est.par$mu, est.par$b)
    attributes(est.par)$par.s.e <-  est.par.se  
    
    class(est.par) <- "eDist"
    
    return(est.par)
  }

#' @rdname Laplace
#' @export lLaplace
## (weighted) (log) likelihood function
lLaplace <- function(x, w=1, mu=0, b=1, params=list(mu, b), logL=TRUE,...){
  if(!missing(params)){
    mu <- params$mu
    b <- params$b}
  # log-likelihood as given by Johnson et.al (p.172) 
  ll <- sum((-abs(x-mu)/b) - log(2*b))
  if(logL){return(ll)} else{return(exp(ll))}
}
