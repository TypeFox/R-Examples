#******************************************************************************#
# Fixed bandwidth algorithm for full and half kernel methods.                  #
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
#  bandwidth      an object of class numeric.                                  #
#                 bandwidth value(s) at which parameters are to be estimated.  #
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
#  Returns a list of matrices                                                  #
#                                                                              #
# betaHat The estimated model coefficients.                                    #
# stdErr  The standard error for each coefficient.                             #
# zValue  The estimated z-value for each coefficient.                          #
# pValue  The p-value for each coefficient.                                    #
#                                                                              #
#******************************************************************************#
kernelFixed <- function(X, 
                        Z, 
                        tau, 
                        bandwidth,
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
  # Determine number of bandwidths provided by user                          #
  #--------------------------------------------------------------------------#
  lbd <- length(bandwidth)

  #--------------------------------------------------------------------------#
  # initialize matrices for parameter estimates and standard deviations      #
  #--------------------------------------------------------------------------#
  bHat <- matrix(data = 0.0, 
                 nrow = lbd, 
                 ncol = nCov,
                 dimnames = list(NULL,colnames(Z)[3L:(nCov+2L)]))

  sdVec <- bHat

  guess <- NULL

  results <- matrix(data = 0.0,
                    nrow = nCov,
                    ncol = 4L,
                    dimnames = list(paste("beta",0L:{nCov-1L},sep=""),
                                    c("estimate","stdErr","z-value","p-value")))
  for( bd in 1L:lbd ) {

    cat("\nBandwidth: ", bandwidth[bd], "\n", sep="")

    bHat[bd,] <- betaEst(Z = Z,
                         X = X, 
                         tau = tau, 
                         h = bandwidth[bd],
                         kType = kType,
                         betaGuess = guess,
                         tol = tol,
                         maxiter = maxiter,
                         scoreFunction = scoreFunction)

    guess <- bHat[bd,]

    argList <- list("beta" = bHat[bd,],
                    "Z" = Z,
                    "X" = X, 
                    "tau" = tau, 
                    "h" = bandwidth[bd],
                    "kType" = kType)

    score <- do.call(what = scoreFunction,
                     args = argList)

    invdU <- try(solve(score$dUdBeta), silent = TRUE)

    if( class(invdU) == 'try-error' ) {
      cat("Unable to invert derivative of estimating equation.\n")
      stop(attr(invdU,"condition"))
    }

    sig <- invdU %*% (score$mMatrix) %*% invdU

    sdVec[bd,] <- sqrt(diag(sig)) 

    results[,1L] <- bHat[bd,]
    results[,2L] <- sdVec[bd,]
    results[,3L] <- bHat[bd,]/sdVec[bd,]
    results[,4L] <- 2.0*pnorm(-abs(results[,3L]))

    print(results)
    cat("\n")

  }

  zv <- bHat/sdVec
  pv <- 2.0*pnorm(-abs(zv))

  return( list( "betaHat" = bHat,
                "stdErr"  = sdVec,
                "zValue" = zv,
                "pValue" = pv ) )

}
