# This function makes the output object for function FPCA

######
# Input:
######  
#  optns: options for FPCA function
#  smcObj: smooth mean curve object
#  scsObj: smooth cov surface object
#  eigObj: eigen analysis object
#  scoresObj: FPC scores object
#  obsGrid: observed time grid
#  workGrid: time grid for smoothed covariance surface
#  rho: regularization parameter for sigma2
#  fitLambda: eigenvalues by least squares fit method
######
# Output: 
######
#  ret: FPCA output object
##########################################################################

MakeResultFPCA <- function(optns, smcObj, mu, scsObj, eigObj, 
                           scoresObj, obsGrid, workGrid, rho=NULL, fitLambda=NULL, inputData){

  if (optns$methodXi == 'CE') {
    xiEst <- t(do.call(cbind, scoresObj[1, ])) 
    xiVar <- scoresObj[2, ]
  } else if (optns$methodXi == 'IN') {
    xiEst <- scoresObj$xiEst
    xiVar <- scoresObj$xiVar
  }
  ret <- list(sigma2 = scsObj$sigma2, 
              lambda = eigObj$lambda, 
              phi = eigObj$phi, 
              xiEst = xiEst, 
              xiVar = xiVar, 
              obsGrid = obsGrid, 
              workGrid = workGrid, 
              mu = spline(x= obsGrid, y= mu, xout=  workGrid)$y, 
              smoothedCov = scsObj$smoothCov, 
              FVE = eigObj$cumFVE[eigObj$kChoosen], 
              cumFVE =  eigObj$cumFVE, 
              fittedCov = eigObj$fittedCov, 
              optns = optns, 
              bwMu = smcObj$bw_mu, 
              bwCov = scsObj$bwCov)

  if (optns$methodXi == 'CE') {
    ret$rho <- rho
  }

  if (optns$fitEigenValues) {
    ret$fitLambda <- fitLambda
  }
  
  if(!optns$lean){
    ret$inputData <- inputData;
  } else {
    ret$inputData <- NULL
  }
  class(ret) <- 'FPCA'
  return(ret)
}
