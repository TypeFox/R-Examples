#' Confidence bands for the sample return time extremogram
#' 
#' @description The function estimates empirical confidence bands for the sample returt time extremogram 
#'              via a permutation procedure under the assumption that the data are independent.
#' @param x Univariate time series (a vector).
#' @param m Number of permutations (an integer).
#' @param type Type of confidence bands. If type=1, it adds all permutations to the sample 
#'             extremogram plot. If type=2, it adds the \code{alpha}/2 and (1-\code{alpha})/2 empirical 
#'             confidence bands for each lag. If type=3, it calculates the lag 1 \code{alpha}/2 and 
#'             (1-\code{alpha})/2 empirical confidence bands lag and uses them for all of the lags.
#' @param exttype  Extremogram type (see  \code{\link{extremogramr}}).                         
#' @param maxlag Number of lags to include in the extremogram (an integer).
#' @param uplevel Quantile of the time series to indicate a upper tail extreme event 
#'                (a number between 0 and 1, default is 1).
#' @param lowlevel Quantile of the time series to indicate a lower tail extreme event 
#'                 (a number between 0 and 1, default is 0).
#' @param start The lag that the extremogram plots starts at (an integer not greater than \code{maxlag}, default is 1).
#' @param alpha Significance level for the confidence bands (a number between 0 and 1, default is 0.05).
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
#' exttype  = 3
#' maxlag   = 70
#' type     = 3
#' m        = 10
#' df       = 3
#' G = extremogram:::garchsim(omega,alpha,beta,n,df)
#' 
#' extremogramr(G, type, maxlag, uplevel, lowlevel, 1, 1)
#' permfnr(G, m, type, exttype, maxlag, uplevel, lowlevel, 1, 0.05)
#' @export 
permfnr = function(x, m, type, exttype, maxlag, uplevel=1, lowlevel=0, start=1, alpha = 0.05 ){
    
  if (type == 1){
    for (i in 1:m){
      pBACC = sample(x)
      cc    = permbootr(pBACC, exttype, uplevel, lowlevel, maxlag)
      lines((start:(maxlag-1)),cc[(start+1):maxlag],col=1,lwd=1)
      
    }
  }
  if (type == 2){
    cc = matrix(0,ncol = m, nrow = maxlag)
    for (i in 1:m){
      pBACC  = sample(x)
      bb = permbootr(pBACC, exttype, uplevel, lowlevel, maxlag)
      cc[,i] = bb[1:maxlag]
    }
    k = dim(cc)[1]
    pocket = matrix(0,ncol=3,nrow=k)
    
    for (i in 1:k){
      pocket[i,1] = quantile(cc[i,],prob = (alpha/2))
      pocket[i,3] = quantile(cc[i,],prob = (1-alpha/2))
    }
    lines(start:(k-1+start),pocket[start:(k-1+start),2],col=1,lwd=2)
    lines(start:(k-1+start),pocket[start:(k-1+start),3],col=1,lwd=2)
  }
  if (type == 3){
    cc = matrix(0,ncol = m, nrow = (start+1))
    for (i in 1:m){
      pBACC  = sample(x)
      bb = permbootr(pBACC, exttype, uplevel, lowlevel,(start+1))
      cc[,i] = bb[1:(start+1)]
    }
    dde = as.numeric()
    gge = as.numeric()
    for (i in 1:maxlag){
      dde[i]  = quantile(cc[(start+1),], prob = (alpha/2))
      gge[i]  = quantile(cc[(start+1),], prob = (1-alpha/2))
    }
    lines((start:(maxlag-1)),dde[(start+1):maxlag],col=1,lwd=2)
    lines((start:(maxlag-1)),gge[(start+1):maxlag],col=1,lwd=2)
  }
}



##================================================================================

# This function calculates the simplified return time extremogram for a univariate 
# time series

permbootr = function(x, type, uplevel=1, lowlevel=0, maxlag) {
  
  # x         = time series
  # type      = 1   extreme event is only upper tail event
  #           = 2   extreme event is only lower tail event
  #           = 3   extreme event is lower and upper tail event
  # uplevel   = upper quantile for extreme event
  # lowlevel  = lower quantile for extreme event
  # maxlag    = the number of lags in the extremogram
  
  
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
  
  aa    = as.matrix(table(return_time))
  for (j in 1:dim(aa)[1]){
    aa[j,1]=aa[j,1]/length(return_time)
  }
  aa = as.double(aa)[1:maxlag]
  aa[is.na(aa)] <- 0
  
  return(aa)
  
}