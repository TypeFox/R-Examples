GetUserCov <- function(optns, obsGrid, cutRegGrid, buff, ymat) {
# If covariance function is provided


  rangeUser <- range(optns$userCov$t)
  rangeCut <- range(cutRegGrid)
  if( rangeUser[1] > rangeCut[1] + buff || 
      rangeUser[2] < rangeCut[2] - buff   ) {
    stop('The range defined by the user provided covariance does not cover the support of the data.')
  }
  
  bwCov  = NULL
  smoothCov = ConvertSupport(fromGrid = optns$userCov$t, cutRegGrid, Cov =  optns$userCov$cov)
  
  if (optns$error) { # error == TRUE
    if (!is.null(optns[['userSigma2']])) {
      sigma2 <- optns[['userSigma2']]
    } else if (optns$dataType %in% c('Dense', 'DenseWithMV')) {
      ord <- 2
      sigma2 <- mean(diff(t(ymat), differences=ord)^2, na.rm=TRUE) / choose(2 * ord, ord)
    } else {
      stop('Use GetSmoothedCovarSurface instead!')
    }
  } else { # error == FALSE
    sigma2 <- NULL
  }

  res <- list(rawCov = NULL,
              smoothCov = (smoothCov + t(smoothCov)) / 2, 
              bwCov = NULL, 
              sigma2 = sigma2, 
              outGrid = cutRegGrid)
  class(res) <- "SmoothCov"  
  return(res)
}
