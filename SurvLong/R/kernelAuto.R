#******************************************************************************#
# Auto Tune algorithm for full and half kernel methods.                        #
#******************************************************************************#
#                                                                              #
# Inputs                                                                       #
#                                                                              #
#  Z              an object of class data.frame.                               #
#                 The structure of the data.frame must be                      #
#                 \{patient ID, time of measurement, measurement(s)\}.         #
#                 Patient IDs must be of class integer or be able to be        #
#                 coerced to class integer without loss of information.        #
#                 Missing values must be indicated as NA.                      #
#                                                                              #
#  X              an object of class data.frame.                               #
#                 The structure of the data.frame must be                      #
#                 \{patient ID, event time, event indicator\}.                 #
#                 Patient IDs must be of class integer or be able to be        #
#                 coerced to class integer without loss of information.        #
#                 Missing values must be indicated as NA.                      #
#                                                                              #
#  tau            an object of class numeric.                                  #
#                 The desired time point.                                      #
#                                                                              #
#  kType          an object of class character indicating the type of          #
#                 smoothing kernel to use in the estimating equation.          #
#                 Must be one of \{"epan", "uniform", "gauss"\}, where         #
#                 "epan" is the Epanechnikov kernel and "gauss" is the         #
#                 Gaussian kernel.                                             #
#                                                                              #
#  tol            an object of class numeric.                                  #
#                 maximum allowed change in parameter estimates, beyond which  #
#                 the parameter estimates are deemed to have converged.        #
#                                                                              #
#  maxiter        an object of class numeric.                                  #
#                 maximum number of iterations allowed to attain convergence   #
#                                                                              #
#  scoreFunction  an object of class character.                                #
#                 the name of the function to be used to calculate the score   #
#                                                                              #
#  Outputs                                                                     #
#                                                                              #
#  Returns a list                                                              #
#                                                                              #
# betaHat The estimated model coefficients.                                    #
# stdErr  The standard error for each coefficient.                             #
# zValue  The estimated z-value for each coefficient.                          #
# pValue  The p-value for each coefficient.                                    #
#                                                                              #
# If the bandwidth is determined automatically, two additional list            #
# elements are returned:                                                       #
#                                                                              #
#  optBW   The estimated optimal bandwidth.                                    #
#  minMSE  The mean squared error at the optimal bandwidth.                    #
#                                                                              #
#******************************************************************************#
kernelAuto <- function(X, 
                       Z, 
                       tau, 
                       kType,
                       tol,
                       maxiter,
                       scoreFunction){

  #--------------------------------------------------------------------------#
  # Process and verify input datasets                                        #
  #--------------------------------------------------------------------------#
  pre <- preprocessInputs(data.x = X, data.z = Z)

  X <- pre$data.x
  Z <- pre$data.z
  nCov <- ncol(Z) - 2L

  rm(pre)

  #--------------------------------------------------------------------------#
  # Determine search space for bandwidth                                     #
  #--------------------------------------------------------------------------#
  effn <- sum(X[,3L])
  limit <- 2.0*(quantile(Z[,2L], 0.75) - quantile(Z[,2L], 0.25))
  bw <- seq(from = limit * (effn)^(-0.7), 
            to = limit * (effn)^(-0.3),
            length.out = 50)
  lbd <- length(bw)

  #--------------------------------------------------------------------------#
  # initialize matrix for parameter estimates                                #
  #--------------------------------------------------------------------------#
  betaHat0 <- matrix(data = 0.0, nrow = lbd, ncol = nCov)

  #--------------------------------------------------------------------------#
  # initialize matrices for variance estimate                                #
  #--------------------------------------------------------------------------#
  hatV <- matrix(data = 0.0, nrow = lbd, ncol = nCov)

  #--------------------------------------------------------------------------#
  # Randomly assign each patient to group 1 or 2                             #
  #--------------------------------------------------------------------------#
  cvlabel <- sample(x = rep( c(1L:2L), length.out = nrow(X) ))

  #--------------------------------------------------------------------------#
  # Estimate parameters at each bandwidth                                    #
  #--------------------------------------------------------------------------#

  guess0 <- NULL
  guess1 <- NULL
  guess2 <- NULL

  for( bd in 1L:lbd ) {

    betaHat0[bd,] <- betaEst(Z = Z,
                             X = X, 
                             tau = tau, 
                             h = bw[bd],
                             kType = kType,
                             betaGuess = guess0,
                             tol = tol,
                             maxiter = maxiter,
                             scoreFunction = scoreFunction)

    guess0 <- betaHat0[bd,]

    tst <- cvlabel == 1L

    betaHat1 <- betaEst(Z = Z,
                        X = X[tst,,drop=FALSE], 
                        tau = tau, 
                        h = bw[bd],
                        kType = kType,
                        betaGuess = guess1,
                        tol = tol,
                        maxiter = maxiter,
                        scoreFunction = scoreFunction)

    guess1 <- betaHat1

    tst <- cvlabel == 2L

    betaHat2 <- betaEst(Z = Z,
                        X = X[tst,,drop=FALSE], 
                        tau = tau, 
                        h = bw[bd],
                        kType = kType,
                        betaGuess = guess2,
                        tol = tol,
                        maxiter = maxiter,
                        scoreFunction = scoreFunction)

    guess2 <- betaHat2

    betaDiff <- betaHat2 - betaHat1

    hatV[bd, ] <- nrow(X) * bw[bd] * ( betaDiff * betaDiff ) * 0.25

  }

  #--------------------------------------------------------------------------#
  # Estimate slop of bias expression; Calculate MSE                          #
  #--------------------------------------------------------------------------#
  hatC <- array(data = 0.0, dim = nCov)
  MSE <- array(data = 0.0, dim = lbd)

  if( scoreFunction == "scoreHalf" ) {

    for( p in 1L:nCov ) {
      hatC[p] <- lm( betaHat0[,p] ~ bw )$coef[2]
    }

    #----------------------------------------------------------------------#
    # MSE = hat(C) * hat(C) * bw^2 + hat(V)                                #
    #----------------------------------------------------------------------#
    MSE <- bw^2 * as.vector(hatC %*% hatC) + rowSums(hatV)

  } else if( scoreFunction == "scoreFull" ) {

    for( p in 1L:nCov ) {
      hatC[p] <- lm( betaHat0[,p] ~ bw^2 )$coef[2]
    }

    #----------------------------------------------------------------------#
    # MSE = hat(C) * hat(C) * bw^4 + hat(V)                                #
    #----------------------------------------------------------------------#
    MSE <- bw^4 * as.vector(hatC %*% hatC) + rowSums(hatV)

  }

  #--------------------------------------------------------------------------#
  # Identify minimum MSE                                                     #
  #--------------------------------------------------------------------------#
  tst <- MSE > 0.0
  if( !any(tst) ) stop("no positive values")
  MSE[!tst] <- NA

  opt_h <- which.min(MSE)

  if( length(opt_h) > 1L ) {
    warning("Multiple minimums. Smallest bandwidth used.")
    opt_h <- opt_h[1]
  }

  if( isTRUE(all.equal(opt_h, bw[1L])) || isTRUE(all.equal(opt_h, bw[lbd])) ) {
    warning("Minimum is at bandwidth boundary.")
  }

  minMSE <-  MSE[opt_h]

  #--------------------------------------------------------------------------#
  # Estimate parameters and standard error at optimal bandwidth              #
  #--------------------------------------------------------------------------#
  bHat <- betaEst(Z = Z,
                  X = X, 
                  tau = tau, 
                  h = bw[opt_h],
                  kType = kType,
                  betaGuess = guess0,
                  tol = tol,
                  maxiter = maxiter,
                  scoreFunction = scoreFunction)

  names(bHat) <- colnames(Z)[3L:(2L+nCov)]

  argList <- list("beta" = bHat,
                  "Z" = Z,
                  "X" = X, 
                  "tau" = tau, 
                  "h" = bw[opt_h],
                  "kType" = kType)

  score <- do.call(what = scoreFunction,
                   args = argList)

  invdU <- try(solve(score$dUdBeta), silent = TRUE)

  if( class(invdU) == 'try-error' ) {
    cat("Unable to invert derivative of estimating equation.\n")
    stop(attr(invdU,"condition"))
  }

  sig <- invdU %*% (score$mMatrix) %*% invdU

  sdVec <- sqrt(diag(sig)) 

  names(sdVec) <- names(bHat)

  #--------------------------------------------------------------------------#
  # Generate results matrix                                                  #
  #--------------------------------------------------------------------------#
  results <- matrix(data = 0.0,
                    nrow = nCov,
                    ncol = 4L,
                    dimnames = list(paste("beta",0L:{nCov-1L},sep=""),
                                    c("estimate","stdErr","z-value","p-value")))

  results[,1L] <- bHat
  results[,2L] <- sdVec
  results[,3L] <- bHat/sdVec
  results[,4L] <- 2.0*pnorm(-abs(results[,3L]))

  cat("\nBandwidth search range: ", bw[1], " - ", bw[lbd], "\n", sep="")
  cat("Optimal bandwidth: ", bw[opt_h], "\n", sep="")
  cat("MSE: ", minMSE, "\n\n", sep="")

  print(results)

  zv <- bHat/sdVec
  pv <- 2.0*pnorm(-abs(zv))
  MSE <- matrix(MSE, ncol = 1L)
  rownames(MSE) <- round(bw,4L)

  return( list( "betaHat" = bHat,
                "stdErr"  = sdVec,
                "zValue"  = zv,
                "pValue"  = pv,
                "optBW"   = bw[opt_h],
                "minMSE" =  minMSE,
                "MSE" = MSE ) )

}
