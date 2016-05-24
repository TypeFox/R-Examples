GCVLwls2DV2 <- function(obsGrid, regGrid, ngrid=NULL, dataType=rcov$dataType, error=rcov$error, kern, rcov, h0=NULL, verbose=FALSE, CV=FALSE, t) {

# TODO:? get the residual values only within truncated regGrid

# Returns: a list of length 2, containing the optimal bandwidth and the gcv score.
# obsGrid: observation points. 
# ngrid: I think this should not be used in the gcv function.
# CV: whether to use CV rather than GCV. Default to FALSE as not using CV. If CV is used use an integer value to specify the number of cross-validation folds. 

# This function computes the optimal bandwidth choice for the covariance surface. 
# function use GCV method by pooling the longitudinal data together. 
# verbose is unused for now
# this follows exactly the matlab 2D gcv selector.
  
  r <- diff(range(obsGrid)) * sqrt(2) # sqrt(2) because the window is circular.

  minBW <- GetMinb(t, dataType=rcov$dataType, obsGrid=obsGrid)

  if (missing(h0)) {
    h0 <- minBW
  }
  
  if (kern == 'gauss') {
    if (is.null(h0))
      stop('Not implemented')
    
    h0 = h0 * 0.2;
  }

  if (is.null(h0))
    stop('the data is too sparse, no suitable bandwidth can be found! Try Gaussian Kernel instead!\n')

  h0 <- min(h0, r/4)
  if (h0 < r/4) {    
    q <- (r / (4 * h0)) ^ (1/9)
  } else if (h0 < r/2) {
    q <- (r / (2 * h0)) ^ (1/9)
  } else if (h0 < r) {
    q <- (r / h0) ^ (1/9)
  } else {
    stop('Data is too sparse. The minimal bandwidth is the range of data')
  }

  bw <- (q ^ (0:9)) * h0 # from h0 to r / 4

  opth <- h0

  leave <- FALSE
  iter <- 0
  maxIter <- 1
  
  if (CV != FALSE) {
# We partition the raw covariance rather than partition the individuals.
    fold <- CV
    partition <- CreateFolds(1:nrow(rcov$tPairs), k=fold)
  }
  
  minBWInvalid <- FALSE
  while (!leave && iter < maxIter) {
    if (minBWInvalid)
      minBW <- bw[1]

    # if (CV == FALSE) {
    
    Scores <- rep(Inf, length(bw))
# try the bandwidths large to small in order to save time due to sparseness in the windows.
    for (i in rev(seq_along(bw))) {
      h <- bw[i]

      if (class(rcov) == 'BinnedRawCov') {
        if (CV == FALSE) # GCV
          Scores[i] <- getGCVscoresV2(h, kern, rcov$tPairs, rcov$meanVals, win=rcov$count, regGrid=regGrid, RSS=rcov$RSS)
        else # CV
          Scores[i] <- getCVscoresV2(partition, h, kern, rcov$tPairs, rcov$meanVals, win=rcov$count, regGrid=regGrid, RSS=rcov$RSS)
      } else {
        if (CV == FALSE) # GCV
          Scores[i] <- getGCVscoresV2(h, kern, rcov$tPairs, rcov$cxxn, regGrid=regGrid)
        else # CV
          Scores[i] <- getCVscoresV2(partition, h, kern, rcov$tPairs, rcov$cxxn, regGrid=regGrid)
      }
      
      if (is.infinite(Scores[i])) {
        minBWInvalid <- TRUE
        if (i < length(bw)) {
          if (minBWInvalid) {
            minBW <- bw[i + 1]
            minBWInvalid <- FALSE
          }
        }
        break
      }
    }
      
    # } else if (CV != FALSE) { 
      # Scores <- sapply(bw, getCVscoresV2, partition=partition, kern=kern, xin=rcov$tPairs, yin=rcov$cxxn) 
    # }
    
    optInd <- which.min(Scores)
    opth <- bw[optInd]
    optgcv <- Scores[optInd]
    
    if (opth >= r - 1e-12) {
      minBW <- r
      leave <- TRUE
      stop('Data is too sparse. The optimal bandwidth equals to the range of input time points. Try Gaussian kernel.')
    }
    # else if (opth < r) {
      
      if (optInd != length(bw) && !is.infinite(optgcv))
        leave <- TRUE            
      else if (is.infinite(optgcv)) {
        if (verbose)
        warning('Data is too sparse, retry with larger bandwidths!')
        h0 <- bw[10] * 1.01
      } else if (opth == bw[length(bw)]) {
        warning('Optimal bandwidth not found in the candidate bandwidths. Retry with larger bandwidths')
        h0 <- bw[9] 
      }
      
      if (!leave) {
        newr <- seq(0.5, 1, by=0.05) * r # ??? this can be quite slow
        ind <- which(newr > h0)[1]
        q <- (newr[ind] / h0) ^ (1/9)
        bw <- q ^ (0:9) * h0
        if (verbose) {
          cat('New bwuserCov candidates:\n')
          print(bw)
        }

        iter <- iter + 1
      }
  # }
    
  }

  ret <- list(h=opth, gcv=optgcv, minBW=minBW)
  if (CV != FALSE)
    names(ret)[2] <- 'cv'
  
  return(ret)

}


getGCVscoresV2 <- function(bw, kern, xin, yin, win=NULL, regGrid, RSS=NULL) {
# ...: passed on to Lwls2D
# RSS: for implementing GCV of binned rcov.
  # browser() 
  if (is.null(win))
    win <- rep(1, length(yin))
    
  fit <- tryCatch(Lwls2D(bw, kern, xin=xin, yin=yin, win=win, xout1=regGrid, xout2=regGrid), error=function(err) {
               warning('Invalid bandwidth. Try enlarging the window size.')
               return(Inf)
  })

  # Catch
  if (is.infinite(fit[1]))
    return(Inf)
    
  # workaround for degenerate case.
  if (any(is.nan(fit)))
    return(Inf)

  obsFit <- interp2lin(regGrid, regGrid, fit, xin[, 1], xin[, 2])
  
  # residual sum of squares
  res <- sum((yin - obsFit) ^ 2 * win)
  if (!is.null(RSS))
    res <- res + sum(RSS)
  
  # kernel at x=0
  k0 <- KernelAt0(kern)
  N <- sum(win)
  r <- diff(range(xin[, 1]))
  bottom <- 1 - (1 / N) * (r * k0 / bw)^2
  GCV <- res / bottom^2

  return(GCV)
}


# k-fold CV
# partition: a list of testset observation indices, returned by caret::createFolds
# ...: passed on to Lwls2D
getCVscoresV2 <- function(partition, bw, kern, xin, yin, win=NULL, regGrid, RSS=NULL) {

  if (is.null(win))
    win <- rep(1, length(yin))
    
  n <- length(yin)
  
  # browser()
  cvSubSum <- sapply(partition, function(testSet) {
    # browser()
    fit <- tryCatch(Lwls2D(bw, kern, xin=xin, yin=yin, win=win, xout1=regGrid, xout2=regGrid, subset=-testSet), error=function(err) {
                 warning('Invalid bandwidth. Try enlarging the window size.')
                 return(Inf)
    })

    # Catch
    if (is.infinite(fit[1]))
      return(Inf)
      
    # workaround for degenerate case.
    if (any(is.nan(fit)))
      return(Inf)

    obsPred <- interp2lin(regGrid, regGrid, fit, xin[testSet, 1], xin[testSet, 2])
    tmpSum <- sum((yin[testSet] - obsPred) ^ 2 * win[testSet])

    # residual sum of squares
    if (!is.null(RSS))
      tmpSum <- tmpSum + sum(RSS[testSet])
      
    return(tmpSum)
  })
  
  return(sum(cvSubSum))
}

# ??? Why is this not the 2D kernel
KernelAt0 <- function(kern) {
  if (kern == 'quar')
    k0 <- 0.9375
  else if (kern == 'epan')
    k0 <- 0.75
  else if (kern == 'rect')
    k0 <- 0.5
  else if (kern == 'gausvar')
    k0 <- 0.498677850501791
  else if (kern == 'gauss')
    k0 <- 0.398942280401433
  else
    stop('Unknown kernel')
    
  return(k0)
}
