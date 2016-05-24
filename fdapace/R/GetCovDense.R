# This function obtains the sample covariance matrix at observed grid
# for dense regular functional data

######
# Input:
######  
#  ymat: n by p matrix of dense regular functional data
#  mu: p-dim vector, estimated cross-sectional mean
#  optns: options for FPCA function
#  y: list of amplitude information
#  t: list of time information
######
# Output: 
######
#  a SmoothCov object containing: 
#    - p by p matrix of sample cov surface estimation on observed grid
#    - NULL for all other entires
##########################################################################

GetCovDense <- function(ymat, mu, optns) {
  if(!(optns$dataType %in% c('Dense', 'DenseWithMV'))){
    stop('Sample Covariance is only applicable for option: dataType = "Dense" or "DenseWithMV"!')
  }
  # if( optns$muCovEstMethod == 'cross-sectional' ){
    n = nrow(ymat)
    m = ncol(ymat)
    if( !is.null(optns$userMu) ){
      ymat = ymat - matrix(rep(times= nrow(ymat), mu), ncol= ncol(ymat), byrow=TRUE)
      K = matrix( rep(0,m^2), m)
      for( i in (1:m)){
        for( j in (1:m)){
          XcNaNindx = which(is.na(ymat[,i]));
          YcNaNindx = which(is.na(ymat[,j]));
          NaNrows = union(XcNaNindx, YcNaNindx);
          # Find inner product of the columns with no NaN values
          indx = setdiff( 1:n, NaNrows)
          K[i,j] =  sum(ymat[indx,i] * ymat[indx,j]) * (1/(n-1-length(NaNrows)));  
        }
      }    
    } else {
      K = cov(ymat, use = 'pairwise.complete.obs') # sample variance using non-missing data
    }
    K = 0.5 * (K + t(K)) # ensure that K is symmetric
  
    if (optns[['error']] == TRUE) {
      # 2nd order difference method for finding sigma2
      if (!is.null(optns[['userSigma2']])) {
        sigma2 <- optns[['userSigma2']]
      } else {
        ord <- 2 
        sigma2 <- mean(diff(t(ymat), differences=ord)^2, na.rm=TRUE) / 
                  choose(2 * ord, ord)
        diag(K) <- diag(K) - sigma2
      }
    } else {
      sigma2 <- NULL
    }
  
    ret = list('rawCov' = NULL, 'smoothCov' = K, 'bwCov' = NULL,
     'sigma2' = sigma2, outGrid = NULL)
    class(ret) = "SmoothCov"
    # Garbage Collection
    gc()
    return(ret)
  # } else {
    # stop('optns$muCovEstMethod is unknown!\n')
  # }
}
