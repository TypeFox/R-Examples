#' Sample cross extremogram
#' 
#' @description The function estimates the sample cross extremogram and creates an extremogram plot.
#' @param a Bivariate time series (n by 2 matrix).
#' @param quant1 Quantile of the first time series to indicate an extreme event (a number between 0 and 1).
#' @param quant2 Quantile of the second time series to indicate an extreme event (a number between 0 and 1).
#' @param maxlag Number of lags to include in the extremogram (an integer).
#' @param type If type=1, the upper tail extremogram is estimated - P(Y>y,X>x). 
#'             If type=2, the lower tail extremogram is estimated - P(Y<y,X<x). 
#'             If type=3, the extremogram is estimated for a lower tail extreme value in the 
#'             first time series and an upper tail extreme value in the second time series -  P(Y>y,X<x). 
#'             If type=4, the extremogram is estimated for a lower tail extreme value in the 
#'             second time series and an upper tail extreme value in the first time series -  P(Y<y,X>x). 
#' @param ploting An extremogram plot. If ploting = 1, a plot is created (default). If ploting = 0,
#'                no plot is created.
#' @param cutoff The cutoff of the y-axis on the plot (a number between 0 and 1, default is 1).
#' @param start The lag that the extremogram plots starts at (an integer not greater than \code{maxlag}, default is 0).
#' @param ... further arguments: plot and axis names.
#' @return Cross extremogram values and a plot (if requested).
#' @references \enumerate{
#'             \item Davis, R. A., Mikosch, T., & Cribben, I. (2012). Towards estimating extremal 
#'             serial dependence via the bootstrapped extremogram. Journal of Econometrics,170(1), 
#'             142-152.
#'             \item Davis, R. A., Mikosch, T., & Cribben, I. (2011). Estimating extremal 
#'             dependence in univariate and multivariate time series via the extremogram.arXiv 
#'             preprint arXiv:1107.5592.}
#' @examples
#' # generate a GARCH(1,1) process
#' omega  = 1
#' alpha1 = 0.1
#' beta1  = 0.6
#' alpha2 = 0.11
#' beta2  = 0.78
#' n      = 1000
#' quant  = 0.95
#' type   = 1
#' maxlag = 70
#' df     = 3
#' G1     = extremogram:::garchsim(omega,alpha1,beta1,n,df)
#' G2     = extremogram:::garchsim(omega,alpha2,beta2,n,df)
#' data   = cbind(G1, G2)
#' 
#' extremogram2(data, quant, quant, maxlag, type, 1, 1, 0)
#' @export 

extremogram2 = function(a, quant1, quant2, maxlag, type, ploting=1, cutoff=1, start=0, ...) {
  
  x=a[,1]
  y=a[,2]
  
  level1 = quantile(a[,1],prob = quant1);level2 = quantile(a[,2],prob = quant2)
  n      = length(a[,1]); rhohat = rep(0,maxlag)
  if (type == 1) { 
    for ( i in 1:maxlag) {
      rhohat[i]=length((1:(n-i))[x[1:(n-i+1)] > level1 & y[i:n]> level2])
      rhohat[i]=rhohat[i]/length((1:n)[x[1:(n-i+1)]>level1])
    }
  }
  else
    if (type == 2) {
      for ( i in 1:maxlag) {
        rhohat[i]=length((1:(n-i))[x[1:(n-i+1)] < level1 & y[i:n]< level2])
        rhohat[i]=rhohat[i]/length((1:n)[x[1:(n-i+1)] < level1])
      }
    }
  else
    if (type == 3) {
      for ( i in 1:maxlag) {
        rhohat[i]=length((1:(n-i))[x[1:(n-i+1)] < level1 & y[i:n]> level2])
        rhohat[i]=rhohat[i]/length((1:n)[x[1:(n-i+1)] < level1])
      }
    }
  else
    if (type == 4) {
      for ( i in 1:maxlag) {
        rhohat[i]=length((1:(n-i))[x[1:(n-i+1)] > level1 & y[i:n]< level2])
        rhohat[i]=rhohat[i]/length((1:n)[x[1:(n-i+1)] > level1])
      }
    }
  
  
  
  if (ploting == 1) {
    plot.extr(rhohat, start = start, maxlag = maxlag, cutoff = cutoff, ...) 
  }
  
  return(rhohat)	
}
