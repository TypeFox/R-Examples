#' Confidence bands for the sample univariate extremogram
#' 
#' @description The function estimates empirical confidence bands for the sample univariate extremogram 
#'              via a permutation procedure under the assumption that the data are independent.
#' @param x Univariate time series (a vector).
#' @param p Quantile of the time series to indicate an extreme event (a number between 0 and 1).
#' @param m Number of permutations (an integer).
#' @param type Type of confidence bands. If type=1, it adds all permutations to the sample 
#'             extremogram plot. If type=2, it adds the \code{alpha}/2 and (1-\code{alpha})/2 empirical 
#'             confidence bands for each lag. If type=3, it calculates the lag 1 \code{alpha}/2 and 
#'             (1-\code{alpha})/2 empirical confidence bands  lag and uses them for all of the lags.
#' @param exttype Extremogram type (see  \code{\link{extremogram1}}).
#' @param maxlag Number of lags to include in the extremogram (an integer).
#' @param alpha Significance level for the confidence bands (a number between 0 and 1, default is 0.05).
#' @param start The lag that the extremogram plots starts at (an integer not greater than \code{maxlag}, default is 1).
#' @return The empirical confidence bands are added to the sample univariate extremogram plot.
#' @references \enumerate{
#'             \item Davis, R. A., Mikosch, T., & Cribben, I. (2012). Towards estimating extremal 
#'             serial dependence via the bootstrapped extremogram. Journal of Econometrics,170(1), 
#'             142-152.
#'             \item Davis, R. A., Mikosch, T., & Cribben, I. (2011). Estimating extremal 
#'             dependence in univariate and multivariate time series via the extremogram.arXiv 
#'             preprint arXiv:1107.5592.}
#' @examples
#' # generate a GARCH(1,1) process
#' omega   = 1
#' alpha   = 0.1
#' beta    = 0.6
#' n       = 1000
#' quant   = 0.95
#' exttype = 1
#' maxlag  = 70
#' df      = 3
#' type    = 3
#' m       = 10
#' G = extremogram:::garchsim(omega,alpha,beta,n,df)
#' 
#' extremogram1(G, quant, maxlag, exttype, 1, 1, 0)
#' permfn1(G, quant, m, type, exttype, maxlag, 1, 0.05)
#' @export 
permfn1 = function(x, p, m, type, exttype, maxlag, start=1, alpha = 0.05){
    
  if (type == 1){
    for (i in 1:m){
      pBACC = sample(x)
      cc    = permboot1(pBACC,p,maxlag=maxlag,type=exttype)
      lines((start:(maxlag-1)),cc[(start+1):maxlag],col=1,lwd=1)
    }
  }
  if (type == 2){
    cc = matrix(0,ncol = m, nrow = maxlag)
    for (i in 1:m){
      pBACC  = sample(x)
      cc[,i] = permboot1(pBACC,p,maxlag=maxlag,type=exttype)
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
      cc[,i] = permboot1(pBACC,p,maxlag=(start+1),type=exttype)
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


##==================================================================================

#the simplified function for the univariate time series ot calculate extremogram, some of the imput is omited

permboot1 = function(x, quant, maxlag, type){
  
  #x       = univariate time series
  #quant   = quantile of the series the extremogram is calculated for e.g. 0.95
  #maxlag  = number of lags
  #type    = 1 if interested in the upper quantile
  #        = 2 if interested in the lower qunatile
  
  
  
  level = quantile(x,prob = quant)
  
  n     = length(x);      rhohat = rep(0,maxlag)
  if (type == 1) { 
    rhohat[1] = 1
    for ( i in 1:(maxlag-1)){
      rhohat[i+1]=length((1:(n-i))[x[1:(n-i)] > level & x[(i+1):n]> level])
      rhohat[i+1]=rhohat[i+1]/length((1:n)[x[1:(n-i)]>level])
    }
  }
  else
    if (type == 2){ 
      rhohat[1] = 1
      for ( i in 1:(maxlag-1)){
        rhohat[i+1]=length((1:(n-i))[x[1:(n-i)] < level & x[(i+1):n]< level])
        rhohat[i+1]=rhohat[i+1]/length((1:n)[x[1:(n-i)]<level])
      }
    }
  
  return(rhohat)
}
