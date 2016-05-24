#' Confidence bands for the sample cross extremogram
#' 
#' @description The function estimates empirical confidence bands for the sample cross extremogram 
#'              via a permutation procedure under the assumption that the data are independent.
#' @param x Bivariate time series (n by 2 matrix).
#' @param p1 Quantile of the first time series to indicate an extreme event (a number between 0 and 1).
#' @param p2 Quantile of the second time series to indicate an extreme event (a number between 0 and 1).
#' @param m Number of permutations (an integer).
#' @param type Type of confidence bands. If type=1, it adds all permutations to the sample 
#'             extremogram plot. If type=2, it adds the \code{alpha}/2 and (1-\code{alpha})/2 empirical 
#'             confidence bands for each lag. If type=3, it calculates the lag 1 \code{alpha}/2 and 
#'             (1-\code{alpha})/2 empirical confidence bands  lag and uses them for all of the lags.
#' @param exttype Extremogram type (see  \code{\link{extremogram2}}).
#' @param maxlag Number of lags to include in the extremogram (an integer).
#' @param alpha Significance level for the confidence bands (a number between 0 and 1, default is 0.05).
#' @param start The lag that the extremogram plots starts at (an integer not greater than \code{maxlag}, default is 1).
#' @return The empirical confidence bands are added to the sample cross extremogram plot.
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
#' alpha1  = 0.1
#' beta1   = 0.6
#' alpha2  = 0.11
#' beta2   = 0.78
#' n       = 1000
#' quant   = 0.95
#' exttype = 1
#' maxlag  = 70
#' df      = 3
#' type    = 3
#' m       = 10
#' G1      = extremogram:::garchsim(omega,alpha1,beta1,n,df)
#' G2      = extremogram:::garchsim(omega,alpha2,beta2,n,df)
#' data    = cbind(G1, G2)
#' 
#' extremogram2(data, quant, quant, maxlag, type, 1, 1, 0)
#' permfn2(data, quant, quant, m, type, exttype, maxlag, 1, 0.05)
#' @export 
permfn2 = function(x, p1, p2, m, type, exttype, maxlag, start=1, alpha = 0.05){
  
  if (type == 1){
    for (i in 1:m){
      pBACC = permatrix(x)
      cc = permboot2(pBACC, p1, p2, maxlag, exttype)
      lines((start:(maxlag-1)),cc[(start+1):maxlag],col=1,lwd=1)
    }
  }
  
  if (type == 2){
    cc = matrix(0,ncol = m, nrow = maxlag)
    for (i in 1:m){
      pBACC  = permatrix(x)
      cc[,i] = permboot2(pBACC, p1, p2, maxlag, exttype)
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
      pBACC  = permatrix(x)
      cc[,i] = permboot2(pBACC, p1, p2, (start+1), exttype)
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

#the simplified function for the bivariate time series ot calculate cross-extremogram, some of the imput is omited

permboot2 = function(a,quant1, quant2, maxlag,type) {
  
  # x1       = bivariate time series (n by 2 matrix)
  # quant1   = quantile of the 1st series
  # quant2   = quantile of the 2nd series
  # maxlag   = number of lags
  # type     = 1 if interested in the upper quantile - P(Y > y | X > x)
  #          = 2 if interested in the lower quantile - P(Y < y | X < x)
  #          = 3 if interested in opposite quantiles - P(Y > y | X < x)
  #          = 4 if interested in opposite quantiles - P(Y < y | X > x)
  
  
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
  
  return(rhohat)
}

##==================================================================================

# Thie function permute the rows of a 2-column or 3-column matrix
# It is used to create premutations for the cross-extremograms of bivariate
# and trivariate time series

permatrix = function(x){
  # x       = matrix of data to be permuted (matrix must be nx2 or nx3)
  
  # Returns: the permuted data
  
  x         = as.matrix(x);nrow1     = dim(x)[1]
  ncol1     = dim(x)[2];mat       = rep(0,ncol1*nrow1)
  sequence  = seq(1,nrow1);junk      = sample(sequence)
  for (i in 1:nrow1){
    mat[(i*ncol1-(ncol1-1)):(i*ncol1)] = x[junk[i],]
  }
  mat = matrix(mat,ncol=ncol1,byrow=T)
  return(mat)
}