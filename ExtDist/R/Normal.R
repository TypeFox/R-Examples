#' @title The Normal Distribution.
#' @description Density, distribution, quantile, random number
#' generation and parameter estimation functions for the normal distribution. Parameter estimation can be based on a weighted or unweighted i.i.d. sample and can be carried out 
#' analytically or numerically.  
#' @rdname Normal
#' @name Normal
#' @aliases dNormal 
#' @aliases pNormal 
#' @aliases qNormal 
#' @aliases rNormal 
#' @aliases eNormal 
#' @aliases lNormal 
#' @aliases sNormal 
#' @aliases iNormal
#' @param params A list that includes all named parameters.
#' @param x,q Vector of quantiles.
#' @param w Optional vector of sample weights.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param X Sample observations.
#' @param mean Location parameter.
#' @param sd Scale parameter.
#' @param method Parameter estimation method.
#' @param logL logical; if TRUE, lNormal gives the log-likelihood, otherwise the likelihood is given.
#' @param ... Additional parameters.
#' 
#' @details If the \code{mean} or \code{sd} are not specified they assume the default values of 0 and 1, respectively.\cr 
#' \cr
#' The \code{dNormal()}, \code{pNormal()}, \code{qNormal()},and \code{rNormal()} functions serve as wrappers of the standard
#' \code{\link[stats]{dnorm}}, \code{\link[stats]{pnorm}}, \code{\link[stats]{qnorm}}, and \code{\link[stats]{rnorm}} functions 
#' in the \pkg{\link{stats}} package. They allow for the parameters to be declared not only as 
#' individual numerical values, but also as a list so parameter estimation can be carried out. \cr
#' \cr
#' The normal distribution has probability density function
#' \deqn{f(x) = \frac{1}{\sqrt{2 \pi} \sigma} e^{-\frac{(x-\mu)^2}{2\sigma^2}}} 
#' where \eqn{\mu} is the mean of the distribution and \eqn{\sigma} is the standard deviation.
#' The analytical unbiased parameter estimations are as given by Johnson et.al (Vol 1, pp.123-128).\cr
#' \cr
#' The log-likelihood function of the normal distribution is given by
#' \deqn{l(\mu, \sigma| x) = \sum_{i}[-0.5 ln(2\pi) - ln(\sigma) - 0.5\sigma^{-2}(x_i-\mu)^2].}
#' The score function and observed information matrix are as given by Casella & Berger (2nd Ed, pp.321-322). 
#' 
#' @return dNormal gives the density, pNormal gives the distribution function,
#' qNormal gives the quantiles, rNormal generates random deviates, and 
#' eNormal estimates the parameters. lNormal provides the log-likelihood function, sNormal the score function, 
#' and iNormal the observed information matrix.  
#'
#' @author Haizhen Wu and A. Jonathan R. Godfrey.\cr
#' Updates and bug fixes by Sarah Pirikahu.
#' 
#' @references Johnson, N. L., Kotz, S. and Balakrishnan, N. (1994) Continuous Univariate Distributions,
#'  volume 1, chapter 13, Wiley, New York.\cr
#'  \cr
#'  Casella, G. and Berger R. L. (2002) Statistical Inference, 2nd Ed, pp.321-322, Duxbury.\cr
#'  \cr
#'  Bury, K. (1999) Statistical Distributions in Engineering, Chapter 10, p.143, 
#'  Cambridge University Press.
#'  
#' @seealso \pkg{\link{ExtDist}} for other standard distributions. 
#' 
#' @examples
#' # Parameter estimation for a distribution with known shape parameters
#' x <- rNormal(n=500, params=list(mean=1, sd=2))
#' est.par <- eNormal(X=x, method="unbiased.MLE"); est.par
#' plot(est.par)
#' 
#' #  Fitted density curve and histogram
#' den.x <- seq(min(x),max(x),length=100)
#' den.y <- dNormal(den.x, mean = est.par$mean, sd = est.par$sd)
#' hist(x, breaks=10, probability=TRUE, ylim = c(0,1.2*max(den.y)))
#' lines(lines(den.x, den.y, col="blue")) # Original data
#' lines(density(x), col="red")           # Fitted curve
#' 
#' # Extracting location and scale parameters
#' est.par[attributes(est.par)$par.type=="location"] 
#' est.par[attributes(est.par)$par.type=="scale"]
#'
#' # Parameter Estimation for a distribution with unknown shape parameters 
#' # Example from: Bury(1999) p.143, parameter estimates as given by Bury are
#' # mu = 11.984 and sigma = 0.067 
#' data <- c(12.065, 11.992, 11.992, 11.921, 11.954, 11.945, 12.029, 11.948, 11.885, 11.997, 
#'          11.982, 12.109, 11.966, 12.081, 11.846, 12.007, 12.011)
#' est.par <- eNormal(X=data, method="numerical.MLE"); est.par
#' plot(est.par)
#' 
#' # log-likelihood, score function and observed information matrix
#' lNormal(data, param = est.par)
#' sNormal(data, param = est.par)
#' iNormal(data, param = est.par)
#'
#' # Evaluating the precision of the parameter estimates by the Hessian matrix
#' H <- attributes(est.par)$nll.hessian; H
#' var <- solve(H)
#' se <- sqrt(diag(var)); se

#' @rdname Normal
#' @export dNormal
dNormal <-function(x, mean=0, sd =1, params = list(mean, sd), ...){
    if(!missing(params)){
      mean <- params$mean
      sd <- params$sd
    }
    dnorm(x, mean = mean, sd = sd)
  }

#' @rdname Normal
#' @export pNormal
pNormal <- function(q, mean=0, sd =1, params = list(mean, sd), ...){
    if(!missing(params)){
      mean <- params$mean
      sd <- params$sd
    }
    pnorm(q, mean = mean, sd = sd)
  }

#' @rdname Normal
#' @export qNormal
qNormal <- function(p, mean=0, sd =1, params = list(mean, sd), ...){
    if(!missing(params)){
      mean <- params$mean
      sd <- params$sd
    }
    qnorm(p, mean = mean, sd = sd)
  }

#' @rdname Normal
#' @export rNormal
rNormal <- function(n, mean=0, sd =1, params = list(mean, sd), ...){
    if(!missing(params)){
      mean <- params$mean
      sd <- params$sd
    }
    rnorm(n, mean = mean, sd = sd)
  }

#' @rdname Normal
#' @export eNormal
eNormal <- function(X,w, method =c("unbiased.MLE", "analytical.MLE", "numerical.MLE"), ...){
  method <- match.arg(method)  
  # Accounting weights
  n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    if(method == "analytical.MLE") {
      # Weighted parameter estimation un-adjusted for bias  
      mu <- mean(w*X)                 
      sigma = sqrt(mean(w*X^2)-mu^2)
      
      # Standard errors of estimates 
      SE.mu <- sigma/sqrt(n)               
      SE.sigma <- sigma/sqrt(2*n) 
      
      est.par <- list(mean = mu, sd = sigma)
      est.par.se <- c(SE.mu, SE.sigma)
    } else if(method == "unbiased.MLE"){
      
      # Weighted parameter estimation adjusted for bias
      mu <- mean(w*X)                     
      sigma <- sd(w*X)   # sigma = sqrt((mean(w*X^2)-mu^2)*n/(n-1))
      
      # Standard errors of parameter estimates 
      SE.mu <- sigma/sqrt(n)          # (Johnson et.al, p.123)         
      SE.sigma <- sigma/sqrt(2*(n-1)) # (Johnson et.al, p.128)
      
      est.par <- list(mean = mu, sd = sigma)
      est.par.se <- c(SE.mu, SE.sigma)
    } else {
      # Numerical parameter estimation.
      method <- "numerical.MLE"
      est.par <- wmle(X=X, w=w, distname = "Normal",
                      initial=list(mean = 0, sd = 1),
                      lower=list(mean = -Inf, sd = 0),
                      upper=list(mean = Inf, sd = Inf))
      
      # SE calculation by taking the square root of the inverted Hessian matrix. 
      est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
      if(class(est.par.se) == "try-error") {
        est.par.se <- rep(NA, length(est.par))
      } 
    }
    
    attributes(est.par)$ob <- X
    attributes(est.par)$weights <- w
    attributes(est.par)$distname <- "Normal"
    attributes(est.par)$method <- method
    attributes(est.par)$par.name <- c("mean","sd")
    attributes(est.par)$par.type <- c("location","scale")
    attributes(est.par)$par.vals <- c(est.par$mean, est.par$sd)
    attributes(est.par)$par.s.e <-  est.par.se
    
    class(est.par) <- "eDist"
    
    return(est.par)
  }

#' @rdname Normal
#' @export lNormal
## (weighted) (log) likelihood function
lNormal <- function(X, w, mean=0, sd =1, params = list(mean, sd), logL = TRUE, ...){ 
   if(!missing(params)){
      mean <- params$mean
      sd <- params$sd
    }
    # Accounting for weights
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    # Likelihood function given by Bury(1999) p.129 taking logs gives the 
    # following log-likelihood function.
    ll <- sum(w*(-log(2*pi)/2-log(sd)-(X-mean)^2/(2*sd^2)))
    #     ll <- sum(w*log(dNormal(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }

#' @rdname Normal
#' @export sNormal
## (weighted) score vectors
sNormal <- function(X, w, mean=0, sd =1, params = list(mean, sd), ...){
    if(!missing(params)){
      mean <- params$mean
      sd <- params$sd
    }
    # Accounting for weights
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    # Score function given by Casella & Berger (p.321)
    score1 <- sum(w*(X-mean)/sd^2)
    score2 <- sum(w*((X-mean)^2-sd^2)/sd^3)
    
    score <- c(score1,score2)
    names(score) <- c("mean","sd")
    return(score)
  }

#' @rdname Normal
#' @export iNormal
## (weighted) (observed) information matrix
iNormal <- function(X, w, mean=0, sd =1, params = list(mean, sd), ...){
    if(!missing(params)){
      mean <- params$mean
      sd <- params$sd
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    # Information matrix given by Casella & Berger (p.322)
    info11 <- sum(w*rep(1/sd^2,n))
    info12 <- sum(w*2*(X-mean)/sd^3)
    info21 <- info12
    info22 <- sum(w*(3*(X-mean)^2-sd^2)/sd^4)
    
    info <- matrix(c(info11,info12,info21,info22), nrow=2,ncol=2)
    rownames(info) <- colnames(info) <- c("mean","sd")
    return(info)
  }
