partial_correlation <- function(mat, method, verbose=FALSE)
{ 
  # calculate inverse covariance matrix
  if (method=="glasso")
  {
    invcov <- as.matrix(huge::huge.select(huge(t(mat),method="glasso",verbose=verbose),criterion="ebic",verbose=verbose)$opt.icov)
  }
  else if (method=="shrinkage")
  {
    # use auto-selected lambda shrinkage parameter
    invcov <- corpcor::invcov.shrink(t(mat), verbose=verbose)
  }
  else if (method=="ridge")
  {
    # use auto-selected lambda penalty parameter based on approximate leave one out cross validation
    invcov <- rags2ridges::optPenalty.LOOCVauto(t(mat), lambdaMin=1e-3,lambdaMax=1e4)$optPrec
  }
  else if (method=="exact")
  {
    # estimate inverse covariance from 
    if (verbose==TRUE) {cat('Calculating exact inverse covariance...\n')}
    invcov <- solve(cov(t(mat)))
  }
  
  # convert inverse covariances to partial correlation coefficients
  if (method!="correlation")
  {
    pcor <- decompose.invcov(invcov)$pr
  }
  else
  {
    if (verbose==TRUE) {cat('Calculating correlations...\n')}
    pcor <- cor(t(mat))
  }
  
  return(pcor)
}