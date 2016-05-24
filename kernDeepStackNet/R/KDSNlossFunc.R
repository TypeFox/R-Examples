# Calculate loss function of KDSN
lossKDSN <- function (parOpt, y, X, gammaPar=1, seedW) {
  
  # Checks
  lenParOpt <- length(parOpt)
  stopifnot((lenParOpt %% 3) ==0)
  levels <- lenParOpt / 3
  stopifnot(length(seedW)==levels)

  # Extract parameters
  Dim <- round (parOpt [seq(1, lenParOpt, 3)])
  sigma <- parOpt [seq(2, lenParOpt, 3)]
  lambda <- parOpt [seq(3, lenParOpt, 3)]
  
  # Fit model
  KDSNfit <- fitKDSN (y=y, X=X, levels=levels, Dim=Dim, sigma=sigma, lambda=lambda, 
                      alpha=rep(0, levels), info=FALSE, seedW=seedW, standX=TRUE)
  
  # Compute fitted values
  fittedVal <- predict(KDSNfit, newx=X)
  #   if(any(fittedVal==1 | fittedVal==0)) {
  #     fittedVal [fittedVal==1] <- 1 - sqrt(.Machine$double.eps)
  #     fittedVal [fittedVal==0] <- sqrt(.Machine$double.eps)
  #   }
  
  # Calculate weight matrix of observations
  varFunc <- varMu (mu=fittedVal)
  linkDeriv <- gDerivMu (mu=fittedVal)
  Wdiag <- calcWdiag (varMu=varFunc, gDerivMu=linkDeriv)
  
  # Extract deviance of last level
  fitDev <- devStandard(preds=fittedVal, ytest=c(y))
  
  # Calculate influence (hat) matrix
  trAmat <- calcTrAFast (X=t(KDSNfit$Output$rftX [[KDSNfit$Input$levels]] $Z), w=Wdiag, lambda=lambda [levels])
  
  # Compute loss
  Loss <- lossGCV (n=dim(X) [1], Dev=fitDev, trA=trAmat, gammaPar=gammaPar)
  # Numerical check of loss
  CheckCondition <- is.infinite(Loss) | 
    is.null(Loss) | 
    is.nan(Loss) |
    is.na(Loss)
  if(any(CheckCondition)) {
    Loss[CheckCondition] <- max(Loss, na.rm=TRUE)*2
  }
  attributes(Loss) <- list(Fit=KDSNfit)
  return(Loss)
}

lossKDSNgivenModel <- function(KDSNfit, y, X, gammaPar=1) {
  # Compute fitted values
  fittedVal <- predict(KDSNfit, newx=X)
  #   if(any(fittedVal==1 | fittedVal==0)) {
  #     fittedVal [fittedVal==1] <- 1 - sqrt(.Machine$double.eps)
  #     fittedVal [fittedVal==0] <- sqrt(.Machine$double.eps)
  #   }
  
  # Calculate weight matrix of observations
  varFunc <- varMu (mu=fittedVal)
  linkDeriv <- gDerivMu (mu=fittedVal)
  Wdiag <- calcWdiag (varMu=varFunc, gDerivMu=linkDeriv)
  
  # Extract deviance of last level
  fitDev <- devStandard(preds=fittedVal, ytest=c(y))

  # Calculate influence (hat) matrix
  trAmat <- calcTrAFast (X=KDSNfit$Output$rftX [[KDSNfit$Input$levels]] $Z, w=Wdiag, lambda=KDSNfit$Input$lambda [KDSNfit$Input$levels])
  
  # Compute loss
  Loss <- lossGCV (n=dim(X) [1], Dev=fitDev, trA=trAmat, gammaPar=gammaPar)
  attributes(Loss) <- list(Fit=KDSNfit)
  return(Loss)
}

# GCV loss
lossGCV <- function (n, Dev, trA, gammaPar=1) {
    stopifnot((n - gammaPar * trA) >= 0)
    n * Dev / (n - gammaPar * trA)^2
}

# Calculate W matrix in P-IRLS
calcWdiag <- function (varMu, gDerivMu) {
  c((varMu * gDerivMu^2) ^ (-1/2))
}

# Variance function of glm
varMu <- function (mu) {
    return(rep(1, length(mu)))
}

# Derivative of link function of glm 
gDerivMu <- function (mu) {
    return(rep(1, length(mu)))
}

# Calculate influence (hat) Matrix
calcTrA <- function (W, X, lambda=0) {
  
  # Matrix calculation
  B_mat <- W %*% X
  
  # QR decomposition
  R_mat <- qr.R(qr(B_mat)) # p x p
  E_mat <- diag(rep(sqrt(lambda), dim(X) [2])) # p x p
  R_E <- rbind(R_mat, E_mat)
  
  # Singular value decomposition
  U_mat_1 <- svd(R_E)$u [dim(R_mat) [1], ]
  
  # Compute trace of hat matrix A
  trA <- sum(diag(U_mat_1 %*% t(U_mat_1)))
  return(trA)
}

calcTrAFast <- function(X, w, lambda=0) {
  Xtilde <- X * w
  XtildeT <- crossprodRcpp(Xtilde)[[1]]
  eigenVal <- getEigenValuesRcpp(XtildeT)
  sum(eigenVal/(eigenVal + lambda))
}

##################################################
# cross validation loss

devStandard <- function (preds, ytest) {
    return(sum((ytest-preds)^2))
}

lossCvKDSN <- function (parOpt, y, X, cvIndex, seedW, lossFunc=devStandard) {
  # Checks
  lenParOpt <- length(parOpt)
  stopifnot((lenParOpt %% 3) ==0)
  levels <- lenParOpt / 3
  stopifnot(length(seedW)==levels)
  
  # Extract parameters
  Dim <- round (parOpt [seq(1, lenParOpt, 3)])
  sigma <- parOpt [seq(2, lenParOpt, 3)]
  lambda <- parOpt [seq(3, lenParOpt, 3)]
  
  # Compute deviance on test data
  cvSamples <- length(cvIndex)
  vecLoss <- vector("numeric", cvSamples)
  for(j in 1:cvSamples) {
    
    # Fit model
    fitKDSN <- fitKDSN (y=y [cvIndex [[j]] ], X=as.matrix(X [cvIndex [[j]], ]), 
                        levels=levels, Dim=Dim, sigma=sigma, lambda=lambda, 
                        alpha=rep(0, levels), info=FALSE, seedW=seedW, standX=TRUE)
    
    # Compute predicted values
    predVal <- predict(fitKDSN, newx=as.matrix(X [-cvIndex [[j]], ]))

    # Compute loss
    vecLoss [j] <- lossFunc(preds=predVal, ytest=y [-cvIndex [[j]] ])
  }
  
  # Calculate average loss and output
  Loss <- mean(vecLoss)
  attributes(Loss) <- list(Fit=fitKDSN)
  return(Loss)
}

# Help function
# Computes the predictive log probability for model based optimization
predLogProb <- function(predMean, predSigma, y, X){
  sum(-log(predSigma)/2 - 
        (y-predMean)^2 / (2 * predSigma) - 
        log(2*pi) / 2)
}
