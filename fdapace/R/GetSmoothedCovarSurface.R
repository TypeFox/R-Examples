# The output outGrid of this function is the (potentially) truncated greid.
GetSmoothedCovarSurface <- function(y, t, mu, obsGrid, regGrid, optns, useBinnedCov=FALSE) {

  dataType <- optns$dataType
  error <- optns$error
  kern <- optns$kernel
  userBwCov <- optns$userBwCov
  methodBwCov <- optns$methodBwCov
  verbose <- optns$verbose
  rotationCut <- optns$rotationCut

# get the truncation of the output grids.
  outPercent <- optns$outPercent
  buff <- .Machine$double.eps * max(abs(obsGrid)) * 10
  rangeGrid <- range(regGrid)
  minGrid <- rangeGrid[1]
  maxGrid <- rangeGrid[2]
  cutRegGrid <- regGrid[regGrid > minGrid + diff(rangeGrid) * outPercent[1] -
                        buff & 
                        regGrid < minGrid + diff(rangeGrid) * outPercent[2] +
                        buff]

  # Get raw covariance, unless user covariance/sigma2 are specified.
  if (is.null(optns[['userCov']]) || 
       (is.null(optns[['userSigma2']]) && error)) {
       
    rcov <- GetRawCov(y, t, obsGrid, mu, dataType, error)
    if (useBinnedCov && methodBwCov == 'CV') {
      stop('If methodBwCov == \'CV\' then we must use the unbinned rcov.')
    }
    
    if (useBinnedCov) {
      rcov <- BinRawCov(rcov)
    }
  } else {
    rcov <- NULL
  }

  # Obtain smoothed covariance.
  if( !is.null(optns$userCov)) { # If covariance function is provided
    rangeUser <- range(optns$userCov$t)
    rangeCut <- range(cutRegGrid)
    if( rangeUser[1] > rangeCut[1] + buff || 
        rangeUser[2] < rangeCut[2] - buff   ) {
      stop('The range defined by the user provided covariance does not cover the support of the data.')
    }
    
    bwCov  = NULL
    smoothCov = ConvertSupport(fromGrid = optns$userCov$t, cutRegGrid, Cov =  optns$userCov$cov)
  
  } else { # estimate the smoothed covariance
    
    if (userBwCov == 0) { # bandwidth selection
      if (methodBwCov %in% c('GCV', 'GMeanAndGCV')) { # GCV
        gcvObj <- GCVLwls2DV2(obsGrid, regGrid, kern=kern, rcov=rcov, verbose=verbose, t=t)
        bwCov <- gcvObj$h
        if (methodBwCov == 'GMeanAndGCV') {
          bwCov <- sqrt(bwCov * gcvObj$minBW)
        }  
      } else if (methodBwCov == 'CV') { # CV 10 fold
        gcvObj <- GCVLwls2DV2(obsGrid, regGrid, kern=kern, rcov=rcov, t=t,
                            verbose=optns$verbose, 
                            CV=optns[['kFoldMuCov']])
        bwCov <- gcvObj$h
      }
    } else if (userBwCov != 0) {
      bwCov <- userBwCov
    }

    if (!useBinnedCov) {
      smoothCov <- Lwls2D(bwCov, kern, xin=rcov$tPairs, yin=rcov$cxxn,
                          xout1=cutRegGrid, xout2=cutRegGrid)
    } else { 
      smoothCov <- Lwls2D(bwCov, kern, xin=rcov$tPairs, yin=rcov$meanVals,
                          win=rcov$count, xout1=cutRegGrid, xout2=cutRegGrid)
    }
  }
  
  # Obtain the error sigma2.
  if (error) {
    if (!is.null(optns[['userSigma2']])) {
      sigma2 <- optns[['userSigma2']]
    } else if (!is.null(optns[['userCov']])) {
      a0 = min(regGrid)
      b0 = max(regGrid)
      lint = b0 - a0
      middleCutRegGrid <- cutRegGrid > a0 + lint * rotationCut[1] - buff & 
                          cutRegGrid < a0 + lint * rotationCut[2] + buff
      if (useBinnedCov) {
        diagT <- rcov[['tDiag']]
        diagVal <- rcov[['diagMeans']]
      } else {
        diagTV <- aggregate(rcov[['diag']][, 2], list(rcov[['diag']][, 1]), mean)
        diagT <- diagTV[, 1]
        diagVal <- diagTV[, 2]
      }
      diagEst <- approx(diagT, diagVal, cutRegGrid[middleCutRegGrid])[['y']]
      sigma2 <- mean(diagEst - diag(smoothCov)[middleCutRegGrid])
      
    } else { # has to estimate sigma2 from scratch
      sigma2 <- PC_CovE(obsGrid, regGrid, bwCov, rotationCut=rotationCut, kernel=kern, rcov=rcov)$sigma2
    }
  
    if(sigma2 < 0) {
      warning("Estimated sigma2 is negative and thus is reset to 1e-6.")
      sigma2 <- 1e-6
    }
    
  } else { # error=FALSE
    sigma2 <- NULL
  }
  
  res <- list(rawCov = rcov, 
              smoothCov = (smoothCov + t(smoothCov)) / 2, 
              bwCov = bwCov, 
              sigma2 = sigma2, 
              outGrid = cutRegGrid)
  class(res) <- "SmoothCov"  
  # Garbage Collection
  gc()
  return(res)
}
