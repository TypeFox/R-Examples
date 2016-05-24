# Random fourier transformation
randomFourierTrans <- function (X, Dim, sigma, seedW=NULL) {
  
  # Draw weights of multivariate normal distribution
  n <- dim(X) [1]
  d <- dim(X) [2]
  X <- t(X)
  set.seed(seedW)
  rW <- rmvnorm (n=d, sigma=diag(Dim) / sigma)
  
  # Calculate Z matrix
  inner <- crossprod(rW, X)
  Z <- rbind(cos(inner), sin(inner)) / sqrt(Dim)
  Output <- list(Z=Z, rW=rW)
  return(Output)
}

# Used for prediction of new inputs given weight matrix
fourierTransPredict <- function (newx, rW) {
  newx <- t(newx)
  inner <- crossprod(rW, newx)
  Z <- rbind(cos(inner), sin(inner)) / sqrt(dim(rW)[2])
  return(Z)
}

# Robust Standardization
# X: Design Matrix
# Standardizes a matrix with median and median absolute deviation
robustStandard <- function (X) {
  colIndex <- 1:dim(X) [2]
  MedianVal <- sapply(colIndex, function (j) median(x=X [, j]))
  MadValue <- sapply(colIndex, function (j) mad (x=X [, j], constant=1))
  TypeScaling <- ifelse(MadValue!=0, "Robust", "Standard")
  for(j in 1:length(TypeScaling)) {
    if(TypeScaling[j]=="Standard") {
      MedianVal[j] <- mean(X[, j])
      MadValue[j] <- sd(X[, j])
    }
    X [, j] <- (X [, j] - MedianVal [j]) / MadValue [j]
  }
  attr(X, which="StandValues") <- list(Location=MedianVal, 
                                       Scale=MadValue, 
                                       Type=TypeScaling)
  return(X)
}

# Fit KDSN
fitKDSN <- function (y, X, levels=1, 
                           Dim=rep(round(sqrt(dim(X) [1]) / 2), levels), 
                           sigma=rep(1, levels), lambda=rep(0, levels), 
                           alpha=rep(0, levels), 
                           info=FALSE, seedW=NULL, standX=TRUE) {
  
  # Conversion to match format
  Y <- matrix(y, ncol=1)
  
  # Checks of Inputs
  stopifnot (is.matrix(Y) & is.matrix(X))
  stopifnot (dim(Y)[2]==1)
  stopifnot (length(sigma) == levels)
  stopifnot (length(lambda) == levels)
  stopifnot (length(alpha) == levels)
  stopifnot (length(seedW) == levels || is.null(seedW))
  
  # Preparation
  rftX <- vector("list", levels)
  linWeights <- vector("list", levels)
  StandValues <- list(Design=list(Location=rep(NA, dim(X)[2]), 
                                  Scale=rep(NA, dim(X)[2]), 
                                  Type=rep(NA, dim(X)[2])), 
                      LinPreds=list(Location=rep(NA, levels-1), 
                                    Scale=rep(NA, levels-1), 
                                    Type=rep(NA, levels-1)))
  
  # Standardization with median and median absolute deviation
  if(standX) {
    X <- robustStandard (X=X)
    StandValues$Design <- attr(X, which="StandValues")
  }
  
  ################
  # Only one level
  if(levels==1) {
    
    # Apply random Fourier transformation
    rftX [[1]] <- randomFourierTrans (X=X, Dim=Dim, sigma=sigma, seedW=seedW)
    
    # Calculate y_new and X_new
    y_new <- rftX [[1]] $Z %*% Y
    X_new <- tcrossprod(rftX [[1]] $Z)
    
    # Estimate linear weights
    linWeights [[1]] <- coef(glmnet (x=X_new, y=y_new, 
                              alpha=alpha, lambda=lambda, 
                              standardize = FALSE), s=lambda [1])
    linWeights [[1]] <- c(as.matrix(linWeights [[1]]))
    
    # Save input
    Input <- list(Y=Y, levels=levels, Dim=Dim, sigma=sigma, 
                  lambda=lambda, alpha=alpha, 
                  info=info, seedW=seedW, standX=standX)
    
    # Arrange output
    Output <- list(rftX=rftX, linWeights=linWeights, StandValues=StandValues)
    All <- list(Output=Output, Input=Input)
    class(All) <- "KDSN"
    return (All)
  }
  
  ###################
  # Two or more levels
  else {
    
    for(l in 1:(levels-1)) {
      
      # Apply random fourier transformation
      rftX [[l]] <- randomFourierTrans (X=X, Dim=Dim [l], 
                                        sigma=sigma [l], seedW=seedW [l])
      
      # Calculate y_new and X_new
      y_new <- rftX [[l]] $Z %*% Y
      X_new <- tcrossprod(rftX [[l]] $Z)
      
      # Estimate linear weights
      linWeights [[l]] <- coef(glmnet (x=X_new, y=y_new, 
                                         alpha=alpha [l], 
                                         lambda=lambda [l], standardize = FALSE), 
                                 s=lambda [l])
      linWeights [[l]] <- c(as.matrix(linWeights [[l]])) 
        
      # Predict linear outputs
      linPreds <- cbind(1, t(rftX [[l]] $Z)) %*% linWeights [[l]]

      # Expand Input Matrix
      if(standX) {
        StandValues$LinPreds$Location[l] <- median(x=linPreds)
        StandValues$LinPreds$Scale[l] <- mad(x=linPreds, constant=1)
        StandValues$LinPreds$Type[l] <- ifelse(StandValues$LinPreds$Scale[l]!=0, "Robust", "Standard")
        if(StandValues$LinPreds$Type[l]=="Standard") {
          StandValues$LinPreds$Location[l] <- mean(linPreds)
          StandValues$LinPreds$Scale[l] <- sd(linPreds)
        }
        linPreds <- (linPreds - StandValues$LinPreds$Location[l]) / StandValues$LinPreds$Scale[l]
      }
      X <- cbind(X, linPreds)
      if(info) {cat("level", l, "fit", "\n")}
      
    }
    
    ################
    # Final level
    
    # Apply random Fourier transformation
    rftX [[levels]] <- randomFourierTrans (X=X, Dim=Dim [levels], 
                                           sigma=sigma [levels], seedW=seedW [levels])
    
    # Calculate y_new and X_new
    y_new <- rftX [[levels]] $Z %*% Y
    X_new <- tcrossprod(rftX [[levels]] $Z)

    # Use glmnet as final level
    linWeights [[levels]] <- coef(glmnet (x=X_new, y=y_new, 
                                     alpha=alpha [levels], 
                                     lambda=lambda [levels], standardize = FALSE),
                                  s=lambda [levels])
    linWeights [[levels]] <- c(as.matrix(linWeights [[levels]]))
    if(info) {cat("level", levels, "fit", "\n")}
    
    # Output
    # Save input
    Input <- list(Y=Y, levels=levels, Dim=Dim, sigma=sigma, lambda=lambda, 
                  alpha=alpha, info=info, seedW=seedW, 
                  standX=standX)
    
    # Arrange output
    Output <- list(rftX=rftX, linWeights=linWeights, StandValues=StandValues)
    All <- list(Output=Output, Input=Input)
    class(All) <- "KDSN"
    return (All)
  }
}

# Predict method for KDSN
predict.KDSN <- function (object, newx, ...) {
  
  noLevels <- object$Input$levels
  if(object$Input$standX) {
    for(j in 1:dim(newx)[2]) {
      newx[, j] <- (newx[, j] - object$Output$StandValues$Design$Location[j]) / 
        object$Output$StandValues$Design$Scale[j]
    }
  }
  
  ############
  # One level
  
  if(noLevels == 1) {
    
    # Transform new data to Fourier space
    tZ <- fourierTransPredict (newx=newx, rW=object$Output$rftX [[1]] $rW)
    
    # Predict linear outputs
    preds <- cbind(1, t(tZ)) %*% object$Output$linWeights [[1]]
    return(preds)
  }
  
  #############
  # More levels
  
  else {
    
    # Predict hidden levels
    
    for(l in 1:(noLevels-1)) {
      
      # Transform new data to Fourier space
      tZ <- fourierTransPredict (newx=newx, rW=object$Output$rftX [[l]] $rW)
      
      # Predict output values
      preds <- cbind(1, t(tZ)) %*% object$Output$linWeights [[l]]

      # Expand inputs
      if(object$Input$standX) {
        preds <- (preds - object$Output$StandValues$LinPreds$Location[l]) / 
          object$Output$StandValues$LinPreds$Scale[l]
      }
      newx <- cbind(newx, preds)
    }
    # Predict output level
    
    # Transform new data to Fourier space
    tZ <- fourierTransPredict (newx=newx, rW=object$Output$rftX [[noLevels]] $rW)

    # Predict output values
    preds <- cbind(1, t(tZ)) %*% object$Output$linWeights [[noLevels]]
    return(c(preds))
  }
}
