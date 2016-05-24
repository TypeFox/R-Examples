#' Sample univariate extremogram
#' 
#' @description The function estimates the sample univariate extremogram and creates an 
#'              extremogram plot.
#' @param x Univariate time series (a vector).
#' @param quant Quantile of the time series to indicate an extreme event (a number between 0 and 1).
#' @param maxlag Number of lags to include in the extremogram (an integer).
#' @param type Extremogram type. If type = 1, the upper tail extremogram is estimated. 
#'             If type = 2, the lower tail extremogram is estimated.
#' @param ploting An extremogram plot. If ploting = 1, a plot is created (default). If ploting = 0,
#'                no plot is created.
#' @param cutoff The cutoff of the y-axis on the plot (a number between 0 and 1, default is 1).
#' @param start The lag that the extremogram plots starts at (an integer not greater than \code{maxlag}, default is 0).
#' @param ... further arguments: plot and axis names.
#' @return Extremogram values and a plot (if requested).
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
#' alpha  = 0.1
#' beta   = 0.6
#' n      = 1000
#' quant  = 0.95
#' type   = 1
#' maxlag = 70
#' df     = 3
#' G = extremogram:::garchsim(omega,alpha,beta,n,df)
#' 
#' extremogram1(G, quant, maxlag, type, 1, 1, 0)
#' @export 

extremogram1 = function(x, quant, maxlag, type, ploting=1, cutoff=1, start=0, ...) {
  
  level = quantile(x,prob = quant)
  
  n     = length(x);      rhohat = rep(0,maxlag, ...)
  if (type == 1)
  { rhohat[1] = 1
    for ( i in 1:(maxlag-1)){
      rhohat[i+1]=length((1:(n-i))[x[1:(n-i)] > level & x[(i+1):n]> level])
      rhohat[i+1]=rhohat[i+1]/length((1:n)[x[1:(n-i)]>level])
    }
  }
  else
    if (type == 2)
    { rhohat[1] = 1
      for ( i in 1:(maxlag-1)){
        rhohat[i+1]=length((1:(n-i))[x[1:(n-i)] < level & x[(i+1):n]< level])
        rhohat[i+1]=rhohat[i+1]/length((1:n)[x[1:(n-i)]<level])
      }
    }
  if (ploting == 1)
  {
    plot.extr(rhohat, start = start, maxlag = maxlag, cutoff = cutoff, ...)
  }
  return(rhohat)  
}

plot.extr = function(x, start, maxlag, cutoff = cutoff, xlab = "lag", ylab = "extremogram",
                      main = "Extremogram plot"){
  
  plot((start:(maxlag-1)), x[(start+1):maxlag], type="n", 
         xlab = xlab, ylab = ylab , ylim=c(0,cutoff), main = main)
  lines((start:(maxlag-1)), x[(start+1):maxlag], col=1, lwd=1, type="h")
  abline((0:(maxlag-1)),0,col=1,lwd=1)
  
}