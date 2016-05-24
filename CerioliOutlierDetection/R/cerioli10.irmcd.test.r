cerioli2010.irmcd.test <- 
# 
# implements the iterated reweighted MCD
# outlier detection method presented in Cerioli (2010)
#
# Author: Christopher G. Green
# Date: 2011-06-23
#
function( datamat, mcd.alpha=max.bdp.mcd.alpha(n,v), 
  signif.alpha=0.05, nsamp = 500, nmini = 300, trace=FALSE ) 
{

  datamat <- as.matrix(datamat)
  if ( any(is.na(datamat)) ) 
    stop("datamat cannot have missing values.")
  n <- nrow(datamat) # number of observations
  v <- ncol(datamat) # dimension

  # steps 1-4: compute the FSRMCD of Cerioli (2010) 
  fsout    <- cerioli2010.fsrmcd.test( datamat, mcd.alpha=mcd.alpha, 
    signif.alpha=signif.alpha, nsamp=nsamp, nmini=nmini, trace=trace )
  # test whether any point is an outlier
  n.sigalf <- length(signif.alpha)
  gamma    <- 1. - ((1. - signif.alpha)^(1./n))
  outliers <- fsout$outliers
  for ( i in 1:n.sigalf ) {
    if ( any(outliers[,i]) ) {
      # test each mahalanobis distance at the gamma[i] level
      outliers[,i] <- fsout$mahdist.rw[,i] > fsout$critvalfcn(gamma[i])
    } else {
      # accept null hypothesis: no outliers
      # don't need to do anything to outliers[,i]
    }
  }

  list(outliers=outliers, mahdist.rw=fsout$mahdist.rw)
}
