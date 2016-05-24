#' Characteristic time / auto-correlation for irregular time series
#' 
#' Estimates characteristic time \eqn{\tau} or auto-correlation \eqn{\rho} from a gappy time series dataset.  It works first by estimating the mean and standard deviation directly from the time series X, using these to standardize the time series, and then optimizes for the likelihood of the value of \eqn{\tau} or \eqn{\rho}. 
#' 
#' @aliases GetRho2
#' @param x  vector of time series values.
#' @param t  vector of times of measurements associated with x.
#' @param tau whether or not to estimate time scale \eqn{\tau} (preferred) or autocorrelation \eqn{\rho}.
#' @return Returns a vector of length two: the estimate and the negative log-likelihood
#' @author Eliezer Gurarie

GetRho <-
function(x,t, tau = TRUE)
  {
    getL <- function(rho)
      GetL(x, t, rho, tau)
    
    if(!tau)
      o <- optimize(getL,lower=0,upper=1,tol=.001, maximum=TRUE)
    else
      o <- optimize(getL,lower=0,upper=t[length(t)] - t[1],tol=.001, maximum=TRUE)
    
    # MUST REMEMBER: 1) rhohat  2) LL
    return(c(rho.hat = o$maximum, LL = o$objective))
  }

GetRho2 <-
  function(x,t, tau = TRUE)
  {
    if(length(x)>2)
    {
      getL <- function(rho)
        GetL(x, t, rho, tau)
      
      if(!tau)
        o <- optimize(getL,lower=0,upper=1,tol=.001, maximum=TRUE)
      else
        o <- optimize(getL,lower=0,upper=t[length(t)] - t[1],tol=.001, maximum=TRUE)
      
      # MUST REMEMBER: 1) rhohat  2) LL
      return(rho.hat = o$maximum)
    } else return(NA)
  }
