#' Confidence bands for the sample cross extremogram
#' 
#' @description The function estimates confidence bands for the sample cross extremogram using the 
#'              stationary bootstrap.
#' @param x Bivariate time series (n by 2 matrix). 
#' @param R Number of bootstrap replications (an integer).
#' @param l Mean block size for stationary bootstrap or mean of the geometric distribution 
#'               used to generate resampling blocks (an integer that is not longer than the 
#'               length of the time series).
#' @param maxlag Number of lags to include in the extremogram (an integer).
#' @param quant1 Quantile of the first time series to indicate an extreme event (a number between 0 and 1).
#' @param quant2 Quantile of the second series to indicate an extreme event (a number between 0 and 1).
#' @param type Extremogram type (see function \code{\link{extremogram2}}).
#' @param par If par = 1, the bootstrap replication procedure will be parallelized. If par = 0, 
#'            no parallelization will be used.
#' @param start The lag that the extremogram plots starts at (an integer not greater than \code{maxlag}, default is 1).
#' @param cutoff The cutoff of the y-axis on the plot (a number between 0 and 1, default is 1).
#' @param alpha Significance level for the confidence bands (a number between 0 and 1, default is 0.05).
#' @param ... further arguments: plot and axis names.
#' @return Returns a plot of the confidence bands for the sample cross extremogram.
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
#' R      = 10
#' l      = 30
#' par    = 0
#' G1     = extremogram:::garchsim(omega,alpha1,beta1,n,df)
#' G2     = extremogram:::garchsim(omega,alpha2,beta2,n,df)
#' data   = cbind(G1, G2)
#' 
#' extremogram2(data, quant, quant, maxlag, type, 1, 1, 0)  
#' bootconf2(data, R, l, maxlag, quant, quant, type, par, 1, 1, 0.05)  
#' @export 
bootconf2 = function(x, R, l, maxlag, quant1, quant2, type, par, start=1, cutoff=1, alpha=0.05, ...) {

if (par == 1){  
  n = parallel::detectCores()
  boot  = boot::tsboot(x, 
                       permboot2, 
                       R, 
                       l=l, 
                       sim="geom", 
                       endcorr=TRUE, 
                       maxlag=maxlag, 
                       quant1=quant1,
                       quant2=quant2, 
                       type=type,
                       parallel="snow",
                       ncpus=n)
  
  tmp = boot[[2]]; 
  mat = tmp[,(2:maxlag)]
  
  k = dim(mat)[2]
  pocket = matrix(0,ncol=3,nrow=k)
  
  for (i in 1:k){
    pocket[i,1] = quantile(mat[,i],prob = (alpha/2))
    pocket[i,2] = mean(mat[,i])
    pocket[i,3] = quantile(mat[,i],prob = (1-alpha/2))
  }
  
  plot.boot(pocket, start = start, cutoff = cutoff, k = k, ...)
}

else{
  boot  = boot::tsboot(x, 
                       permboot2, 
                       R, 
                       l=l, 
                       sim="geom", 
                       endcorr=TRUE, 
                       maxlag=maxlag, 
                       quant1=quant1,
                       quant2=quant2, 
                       type=type)
  
  tmp = boot[[2]]; 
  mat = tmp[,(2:maxlag)]
  
  k = dim(mat)[2]
  pocket = matrix(0,ncol=3,nrow=k)
  
  for (i in 1:k){
    pocket[i,1] = quantile(mat[,i],prob = (alpha/2))
    pocket[i,2] = mean(mat[,i])
    pocket[i,3] = quantile(mat[,i],prob = (1-alpha/2))
  }
  
  plot.boot(pocket, start = start, cutoff = cutoff, k = k, ...)
}
  
}


