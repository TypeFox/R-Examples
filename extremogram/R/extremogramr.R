#' Sample return time extremogram 
#' 
#' @description The function estimates the sample return time extremogram and creates an extremogram plot.
#' @param x Univariate time series (a vector).
#' @param type Extremogram type. If type = 1, the upper tail extremogram is estimated. If type = 2,
#'             the lower tail extremogram is estimated. If type = 3, both upper and lower tail 
#'             extremogram is estimated.
#' @param maxlag Number of lags to include in the extremogram (an integer).
#' @param uplevel Quantile of the time series to indicate a upper tail extreme event 
#'                (a number between 0 and 1, default is 1).
#' @param lowlevel Quantile of the time series to indicate a lower tail extreme event 
#'                 (a number between 0 and 1, default is 0).
#' @param histogram An extremogram plot. If histogram = 1, a plot is created (default). If histogram = 0,
#'                 no plot is created.
#' @param cutoff The cutoff of the y-axis on the plot (a number between 0 and 1, default is 1).
#' @param ... further arguments: plot and axis names.
#' @return Extremogram values, return time for extreme events, mean return time and a plot (if requested).
#' @references \enumerate{
#'             \item Davis, R. A., Mikosch, T., & Cribben, I. (2012). Towards estimating extremal 
#'             serial dependence via the bootstrapped extremogram. Journal of Econometrics,170(1), 
#'             142-152.
#'             \item Davis, R. A., Mikosch, T., & Cribben, I. (2011). Estimating extremal 
#'             dependence in univariate and multivariate time series via the extremogram.arXiv 
#'             preprint arXiv:1107.5592.}
#' @examples
#' # generate a GARCH(1,1) process
#' omega    = 1
#' alpha    = 0.1
#' beta     = 0.6
#' n        = 1000
#' uplevel  = 0.95
#' lowlevel = 0.05
#' type     = 3
#' maxlag   = 70
#' df       = 3
#' G = extremogram:::garchsim(omega,alpha,beta,n,df)
#' 
#' extremogramr(G, type, maxlag, uplevel, lowlevel, 1, 1)
#' @export 
extremogramr = function(x, type, maxlag, uplevel=1, lowlevel=0, histogram=1, cutoff=1, ...) {  
  
  n   = length(x); mat = mat = rep(0,n)
  if (type == 1){ 
    uplevel  = quantile(x,prob = uplevel)
    for (i in 1:n){
      mat[i] = ifelse(x[i] > uplevel,1,0)
    }
  }
  else
    if (type == 2){ 
      lowlevel  = quantile(x,prob = lowlevel)
      for (i in 1:n){
        mat[i] = ifelse(x[i] < lowlevel,1,0)
      }
    }
  else
    if (type == 3) { 
      uplevel  = quantile(x,prob = uplevel)
      lowlevel = quantile(x,prob = lowlevel)
      for (i in 1:n){
        mat[i] = ifelse(x[i] > uplevel || x[i] < lowlevel,1,0)
      }
    }
  sequence = seq(1,n); gar = cbind(sequence,mat)
  junk = matrix(0,ncol=2,nrow=n)
  for  (i in 1:n){
    if (mat[i] == 1){
      junk[i,] = gar[i,]
    }
  }
  ind <- rowSums(junk == 0) != ncol(junk)
  junk = junk[ind, ]
  n = dim(junk)[1];return_time = rep(0,n-1)
  for (i in 1:n-1){
    return_time[i] = (junk[i+1,1] - junk[i,1])
  }
  if (histogram == 1) {	
    plot.return(return_time, maxlag = maxlag, cutoff = cutoff,...)
  }
  
  aa    = as.matrix(table(return_time))
  for (j in 1:dim(aa)[1]){
    aa[j,1]=aa[j,1]/length(return_time)
  }
  aa = as.double(aa)[1:maxlag]
  aa[is.na(aa)] <- 0
  
  return(list(aa, return_time, mean(return_time)))
  
}


plot.return = function(x, maxlag = maxlag, cutoff = cutoff, ylab = "extremogram", 
                       xlab = "lag", main = "Return time extremogram plot"){
  MASS::truehist(x, 
                 nbins = max(x),
                 xlim=c(0,(maxlag-1)),
                 col=0,
                 prob=TRUE,
                 ylim=c(0,cutoff),
                 main = main,
                 ylab = ylab,
                 xlab = xlab)
}