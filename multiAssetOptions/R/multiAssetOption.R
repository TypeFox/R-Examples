multiAssetOption <- function(X) {
  
  # check length of input vectors
  if (X$opt$nAsset < 1 | length(X$opt$pcFlag) != X$opt$nAsset
      | length(X$opt$strike) != X$opt$nAsset | length(X$opt$vol) 
      != X$opt$nAsset | length(X$opt$q) != X$opt$nAsset | 
      length(X$fd$kMult) != X$opt$nAsset | length(X$fd$leftBound) 
      != X$opt$nAsset | length(X$fd$density) != X$opt$nAsset | 
      length(X$fd$kShift) != X$opt$nAsset | length(X$fd$m) != X$opt$nAsset) {
    stop("invalid input vector length. All input vectors must be same length")
  }
  
  # calculate penalty size
  if (X$opt$exerType == 0) {
    large <- 0
  } else if (X$opt$exerType == 1) {
    large <- 1 / X$fd$tol
  } else {
    stop("invalid exercise type. Must be 0 or 1")
  }
  
  # determine timestep size
  mT <- ttm <- X$opt$ttm
  if (X$time$tsType == 0) {
    dt <- ttm / X$time$N
  } else if (X$time$tsType == 1) {
    dt <- X$time$dtInit
  } else {
    stop("invalid timestep type. Must be 0 or 1")
  }
  
  # create underlying spatial grid and calculate length of S
  S <- dimS <- rightBound <- c()
  for (i in 1:X$opt$nAsset) {
    if (X$fd$kMult[i] == 0) {
      rightBound <- c(rightBound, X$opt$strike[i] * exp(sqrt(2 * 
                      X$opt$vol[i]^2 * ttm * abs(log(1/1000)))))
    } else if (X$fd$kMult[i] > 0) {
      rightBound <- c(rightBound, X$fd$kMult[i]*X$opt$strike[i])
    } else {
      stop("invalid kMult parameter. Must be greater than or equal to zero")
    }
    S <- c(S, list(nodeSpacer(X$opt$strike[i], X$fd$leftBound[i],
            rightBound[i], X$fd$m[i]+1, X$fd$density[i],
            X$fd$kShift[i])))
    dimS <- c(dimS, length(S[[i]]))
  }
  
  # calculate option payoff array (initial condition in backwards time)
  aPayoff <- payoff(X$opt$payType, X$opt$pcFlag, X$opt$strike, S)
  
  # create state and payoff vectors from payoff array (classically ordered)
  mV <- vV <- vPayoff <- cbind(c(aPayoff))
  
  # create sparse difference matrix
  mO <- matrixFDM(S, X$opt$rf, X$opt$q, X$opt$vol, X$opt$rho)

  # create Rannacher smoothing index
  smoothIndex <- 0
  
  # loop until ttm = 0
  while (ttm > 0) {
    
    # create matrix M
    mM <- mO * dt
    
    # set penalty iteration vector
    vViter <- vV
    
    # Rannacher smoothing
    if (smoothIndex < X$fd$maxSmooth) {
      tempTheta <- 1
      smoothIndex <- smoothIndex + 1
    } else {
      tempTheta <- X$fd$theta
    }
    
    # precalculate (I + (1 - X$fd$theta)M)Vn
    vB <- (bandSparse(prod(dimS), prod(dimS), c(0), matrix(1, prod(dimS), 1)) +
           (1 - tempTheta) * mM) %*% vV
    
    # reset iteration checks
    checkP <- checkV <- checkLoop <- checkIter <- iterIndex <- 0
    
    # begin penalty iteration
    while (checkIter < 1) {
      
      # calculate old penalty matrix
      P0 <- large * (vPayoff - vViter > 0)
      
      # set up system Av=b
      mA <- bandSparse(prod(dimS), prod(dimS), c(0), matrix(1, prod(dimS), 1)) -
              tempTheta * mM + bandSparse(prod(dimS), prod(dimS), c(0), 
              matrix(P0, prod(dimS), 1))
      vB1 <- vB + P0 * vPayoff
      
      # solve system
      vVnew <- Matrix::solve(mA, vB1, sparse=TRUE)
      
      # calculate new penalty matrix
      P1 <- large * (vPayoff - vVnew > 0)
      
      # check if old and new penalty matrices the same
      if (sum(P1==P0) == prod(dimS)) {
        checkP <- 1
      } else {
        checkP <- 0
      }
      
      # check if L-infinity norm is less than tolerance
      if (max(abs(vVnew - vViter) / pmax(1, abs(vVnew))) < X$fd$tol) {
        checkV <- 1
      } else {
        checkV <- 0
      }
      
      # check if loop count exceeds loop threshold
      if (iterIndex >= X$fd$maxIter - 1) {
        checkLoop <- 1
      } else {
        checkLoop <- 0
      }
      
      # reset iteration state vector as new state vector
      vViter <- vVnew
      
      # advance iteration loop count
      iterIndex <- iterIndex + 1
      
      # check for iteration loop exit
      checkIter <- checkP + checkV + checkLoop
      
      # end iteration loop
    }
    
    # decay time to maturity by dt
    ttm <- ttm - dt
    
    # adaptive timestep calculation
    if (X$time$tsType == 1) {
      dt <- min(X$time$dNorm / (abs(vVnew - vV) / pmax(X$time$D, 
              abs(vVnew), abs(vV)))) * dt
    }

    # check for ttm = 0
    if (ttm - dt < 0) {
      dt <- ttm
    }

    # reset state vector vV from vector vVnew
    vV <- vVnew
    
    # append new state vector and ttm to history
    mV <- cbind(mV, vV)
    mT <- c(mT, ttm)
    
    # update user on time
    message(paste("time = ", round(ttm, digits = 4)))

    # end time loop
  }

  # return option value matrix
  list("value" = mV, "S" = S, "dimS" = dimS, "time" = mT)
  
}