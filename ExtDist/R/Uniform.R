#' @title The Uniform Distribution.
#' @description Density, distribution, quantile, random number 
#' generation and parameter estimation functions for the uniform distribution on the interval 
#' \eqn{[a,b]}.
#' Parameter estimation can be based on an unweighted i.i.d. sample only and can be performed analytically or  
#' numerically. 
#' @rdname Uniform
#' @name Uniform
#' @aliases dUniform 
#' @aliases pUniform 
#' @aliases qUniform
#' @aliases rUniform 
#' @aliases eUniform
#' @aliases lUniform
#' @details If \code{a} or \code{b} are not specified they assume the default values of 0 and 1, respectively.\cr
#' \cr
#' The \code{dUniform()}, \code{pUniform()}, \code{qUniform()},and \code{rUniform()} functions serve as wrappers of the standard
#' \code{\link[stats]{dunif}}, \code{\link[stats]{punif}}, \code{\link[stats]{qunif}}, and \code{\link[stats]{runif}} functions 
#' in the \pkg{\link{stats}} package. They allow for the parameters to be declared not only as 
#' individual numerical values, but also as a list so parameter estimation can be carried out. \cr
#' \cr
#' The uniform distribution has probability density function
#' \deqn{p_x(x) = 1/(b-a)}
#' for \eqn{a \le x \le b}. The analytic maximum likelihood parameter estimates are as given by 
#' \href{http://www.itl.nist.gov/div898/handbook/eda/section3/eda3662.htm}{Engineering Statistics Handbook}.
#' The method of moments parameter estimation option is also avaliable and the estimates are as given by Forbes et.al (2011), p.179.\cr
#' \cr 
#' The log-likelihood function for the uniform distribution is given by
#' \deqn{l(a,b|x) = -n log(b-a)}  
#' 
#' @param params A list that includes all named parameters.
#' @param x,q A vector of quantiles.
#' @param w An optional vector of sample weights.
#' @param p A vector of probabilities.
#' @param n Number of observations.
#' @param X Sample observations.
#' @param a,b Boundary parameters.
#' @param method Parameter estimation method.
#' @param logL logical;if TRUE, lUniform gives the log-likelihood, otherwise the likelihood is given.
#' @param ... Additional parameters.
#'
#' @return dUniform gives the density, pUniform the distribution function,
#' qUniform the quantile function, rUniform generates random deviates, and 
#' eUniform estimates the parameters. lUniform provides the log-likelihood function. 
#' @seealso \pkg{\link{ExtDist}} for other standard distributions.
#' @author Haizhen Wu and A. Jonathan R. Godfrey.\cr
#' Updates and bugfixes by Sarah Pirikahu.
#' 
#' @references Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions,
#'  volume 2, chapter 26, Wiley, New York.\cr
#'  \cr
#'  \href{http://www.itl.nist.gov/div898/handbook/eda/section3/eda3662.htm}{Engineering Statistics Handbook}\cr
#'  \cr
#'  Forbes, C. Evans, M. Hastings, N. & Peacock, B. (2011) Statistical Distributions, 4th Ed, chapter 40, Wiley, New Jersey.\cr
#'  
#' @note The analytical maximum likelihood estimation of the parameters \eqn{a} and \eqn{b} is calculated using the range and
#' mid-range of the sample. Therefore, only unweighted samples are catered for in the eUniform distribution when the method 
#' \code{analytic.MLE} is selected.
#'  
#' @examples
#' # Parameter estimation for a distribution with known shape parameters
#' X <- rUniform(n=500, a=0, b=1)
#' est.par <- eUniform(X, method="analytic.MLE"); est.par
#' plot(est.par)
#' 
#' # Histogram and fitted density
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dUniform(den.x,a=est.par$a,b=est.par$b)
#' hist(X, breaks=10, probability=TRUE, ylim = c(0,1.2*max(den.y)))
#' lines(den.x, den.y, col="blue")  # Original data
#' lines(density(X), lty=2)         # Fitted curve
#' 
#' # Extracting boundary parameters
#' est.par[attributes(est.par)$par.type=="boundary"]
#' 
#' # log-likelihood
#' lUniform(X,param = est.par)
#' 
#' # Example of parameter estimation for a distribution with 
#' # unknown parameters currently been sought after.

#' @rdname Uniform
#' @export dUniform
dUniform <-
  function(x, a=0, b =1, params = list(a, b),...){
    if(!missing(params)){
      a <- params$a
      b <- params$b
    }
    out <- dunif(x, min = a, max = b)
    return(out)
  }


#' @rdname Uniform
#' @export pUniform
pUniform <- 
  function(q, a=0, b =1, params = list(a, b),...){
    if(!missing(params)){
      a <- params$a
      b <- params$b
    }
    out <- punif(q, min = a, max = b)
    return(out)
  }

#' @rdname Uniform
#' @export qUniform
qUniform <- 
  function(p, a=0, b =1, params = list(a, b),...){
    if(!missing(params)){
      a <- params$a
      b <- params$b
    }
    out <- qunif(p, min = a, max = b)
    return(out)
  }

#' @rdname Uniform
#' @export rUniform
rUniform <- 
  function(n, a=0, b =1, params = list(a, b),...){
    if(!missing(params)){
      a <- params$a
      b <- params$b
    }
    out <- runif(n, min = a, max = b)
    return(out)
  }

#' @rdname Uniform
#' @export eUniform
eUniform <-     
  function(X,w, method =c("analytic.MLE", "moments", "numerical.MLE"),...){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    if(method=="analytic.MLE"){    
        # Analytic maximum likelihood estimates as given by http://www.itl.nist.gov/div898/handbook/eda/section3/eda3662.htm
        # Weighted parameter estimation can not be carried out for the analytic MLE method as the parameter estimates are 
        # estimated by the range and mid-range.
    a <- ((max(X) + min(X))/2) - 0.5*(abs(max(X)-min(X)))   
    b <- ((max(X) + min(X))/2) + 0.5*(abs(max(X)-min(X)))
    
    est.par <- list(a = a, b = b)
    est.par.se <- rep(NA, length(est.par))
    }
    else if(method == "moments"){
      # Parameter estimates using the method of moments estimates as given by Forbes et.al (2011), chapeter 40, p.179
      a <- mean(X)-sqrt(3)*sd(X)
      b <- mean(X)+sqrt(3)*sd(X) 
        
      est.par <- list(a = a, b = b)
      est.par.se <- rep(NA, length(est.par))
    }
    else
    {method <- "numerical.MLE"  
      # Numerical parameter estimation
	  d <- (max(X)-min(X))
    est.par <- wmle(X=X, w=w, distname = "Uniform",
                      initial=list(a = min(X)-0.1*d, b = max(X)+0.1*d),
                      lower=list(a = -Inf, b = max(X)),
                      upper=list(a = min(X), b = Inf))
      
      est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
      if(class(est.par.se) == "try-error") {
        est.par.se <- rep(NA, length(est.par))
      } 
    }
    
    attributes(est.par)$ob <- X
    attributes(est.par)$weights <- w
    attributes(est.par)$distname <- "Uniform"
    attributes(est.par)$method <- method
    attributes(est.par)$par.name <- c("a","b")
    attributes(est.par)$par.type <- c("boundary","boundary")
    attributes(est.par)$par.vals <- c(est.par$a, est.par$b)
    attributes(est.par)$par.s.e <-  est.par.se  
    
    class(est.par) <- "eDist"
    
    return(est.par)
  }

#' @rdname Uniform
#' @export lUniform
## (weighted) (log) likelihood function
lUniform <- 
  function(X, w, a=0, b =1, params = list(a, b), logL = TRUE,...){
    if(!missing(params)){
      a <- params$a
      b <- params$b
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dUniform(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }