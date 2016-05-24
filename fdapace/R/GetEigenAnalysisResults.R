# phi: a nRegGrid * no_FVE
# The input smoothCov is possibly truncated.

GetEigenAnalysisResults <- function(smoothCov, regGrid, optns, muWork = NULL) {

  maxK <- optns$maxK
  FVEthreshold <- optns$FVEthreshold
  verbose <- optns$verbose
  
  gridSize <- regGrid[2] - regGrid[1]
  numGrids <- nrow(smoothCov)
  
  eig <- eigen(smoothCov)

  positiveInd <- eig[['values']] >= 0
  d <- eig[['values']][positiveInd]
  eigenV <- eig[['vectors']][, positiveInd, drop=FALSE]

  if (maxK < length(d)) {
    if (optns[['verbose']]) {
      cat(sprintf("At most %d number of PC can be selected, thresholded by `maxK` = %d. \n", length(d), maxK)) 
    }
    
    d <- d[1:maxK]
    eigenV <- eigenV[, 1:maxK, drop=FALSE]
  }

  # thresholding for corresponding FVE option 
  #(not before to avoid not being able to reach the FVEthreshold when pos eigenvalues > maxk)
  # i.e. default FVE 0.9999 outputs all components remained here.
  FVE <- cumsum(d) / sum(d) * 100  # cumulative FVE for all available eigenvalues from fitted cov
  no_opt <- min(which(FVE >= FVEthreshold * 100)) # final number of component chosen based on FVE
  
  # normalization
  if (is.null(muWork)) {
    muWork = 1:dim(eigenV)[1]
  }
  
  phi <- apply(eigenV, 2, function(x) {
                    x <- x / sqrt(trapzRcpp(regGrid, x^2)) 
                    if ( 0 <= cov(x, muWork) )
                      return(x)
                    else
                      return(-x)
  })
  lambda <- gridSize * d;

  fittedCov <- phi %*% diag(x=lambda, nrow = length(lambda)) %*% t(phi)

  return(list(lambda = lambda[1:no_opt], phi = phi[,1:no_opt, drop=FALSE], 
              cumFVE = FVE, kChoosen=no_opt, fittedCov=fittedCov))
}
