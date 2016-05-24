#' Obtain log-likelihood and parameter estimates for a given break point.
#' 
#' Takes a time series with values \code{x} obtained at time \code{t} and a time break \code{tbreak}, and returns the estimates of \eqn{\mu}, \eqn{\sigma} and \eqn{\tau} (or \eqn{\rho}) as well as the negative log-likelihood of those estimates before and after the break. Mostly for use internally within \code{\link{GetBestBreak}}.
#' 
#' @param x  vector of time series values.
#' @param t  vector of times of measurements associated with x.
#' @param tbreak breakpoint to test (in terms of the INDEX within "t" and "x", not actual time value).
#' @param ... additional parameters to pass to \code{\link{GetRho}}.
#' 
#' @return a vector containing the parameters and the negative log-likelihoods in order: \code{mu1, sigma1, tau1, LL1, mu2, sigma2, tau2, LL2}
#' @seealso \code{\link{GetBestBreak}} uses this function, while this function uses \code{\link{GetRho}}
#' @author Eliezer Gurarie


GetDoubleL <-
function(x,t,tbreak, ...)
  {  
    x1 <- x[1:tbreak]  
    x2 <- x[tbreak:length(x)]
    
    t1 <- t[1:tbreak]
    t2 <- t[tbreak:length(t)]

    o1<-GetRho(x1,t1, ...)
    o2<-GetRho(x2,t2, ...)
    
    # remember: o returns the estimate AND the log-likelihood
    
		return(c(mu1 = mean(x1), s1 = sd(x1), o1, 
						 mu2 = mean(x2), s2 = sd(x2), o2))
  }