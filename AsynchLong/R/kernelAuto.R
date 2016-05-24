#----------------------------------------------------------------------#
# kernelAuto : autotune bandwidths for all methods                     #
#----------------------------------------------------------------------#
#                                                                      #
# data.x    : Matrix of covariates. It is assumed that the first column#
#             contains integer patient IDs, the second column contains #
#             the time of measurement, and the remaining columns       #
#             contain the values of the covariates.                    #
#                                                                      #
# data.y    : Matrix of responses. It is assumed that the first column #
#             contains integer patient IDs, the second column contains #
#             the time of measurement, and the third column contains   #
#             the value of the measurement.                            #
#                                                                      #
# kType     : a character. One of "epan", "uniform", or "gauss".       #
#             Specifies the form of the kernel function.               #
#                                                                      #
# lType     : a character. One of "identity", "log", "logistic".       #
#             Specifies the form of the link function.                 #
#                                                                      #
# times     : a numeric or vector object containing time points        #
#                                                                      #
# distanceFunction : a character object giving the distance function   #
#             to be used.                                              #
#                                                                      #
# nCores    : a numeric. Number of cores to use for auto-tune          #
#             implementation                                           #
#                                                                      #
#----------------------------------------------------------------------#
#                                                                      #
# Returns a list. Values are for optimal bandwidth                     #
#                                                                      #
#   betaHat : matrix, row k contains parameter estimates for kth time  #
#             point                                                    #
#   sd      : matrix, row k contains standard deviations for kth time  #
#             point                                                    #
#   zValue  : matrix, row k contains z-values for kth time point       #
#   pValue  : matrix, row k contains p-values for kth time point       #
#   optBW   : matrix, row k contains optimal bw for kth time point     #
#   minMSE  : matrix, row k contains minimum standard error for kth    #
#             time point                                               #
#                                                                      #
#----------------------------------------------------------------------#
kernelAuto <- function(data.x, 
                       data.y,
                       kType, 
                       lType,
                       time,
                       distanceFunction,
                       nCores, ...){

  #------------------------------------------------------------------#
  # Number of covariates (assumes column of ids and measurement times#
  #------------------------------------------------------------------#
  nCov <- ncol(data.x) - 2L
  patientIDs <- sort(unique(data.y[,1L]))
  nPatients <- length(patientIDs)

  #------------------------------------------------------------------#
  # For each response, determine which covariates are for the same   #
  # sample and the number of measurements.                           #
  #------------------------------------------------------------------#
  xIs <- list()
  for( i in 1L:nrow(data.y) ) {
    xIs[[i]] <- list()
    xIs[[i]]$v <- which( data.x[,1L] == data.y[i,1L] )
    xIs[[i]]$n <- as.integer(round(length(xIs[[i]]$v),0L))
  }

  #------------------------------------------------------------------#
  # For each response, determine which response measurements are for #
  # the same sample and the number of measurements.                  #
  #------------------------------------------------------------------#
  yIs <- list()
  for( i in 1L:nPatients ) {
    yIs[[i]] <- list()
    yIs[[i]]$v <- which( data.y[,1L] == patientIDs[i] )
    yIs[[i]]$n <- as.integer(round(length(yIs[[i]]$v),0L))
  }

  #------------------------------------------------------------------#
  # Calculate range of bandwidths to be considered.                  #
  #------------------------------------------------------------------#
  range <- c(data.x[, 2L], data.y[,2L])

  bw <- seq(from = 2*(stats::quantile(x = range, probs = 0.75) - 
                      stats::quantile(x = range, probs = 0.25)) * 
                   nPatients^(-0.7),
            to =   2*(stats::quantile(x = range, probs = 0.75) - 
                      stats::quantile(x = range, probs = 0.25)) * 
                   nPatients^(-0.3), 
            length = 50)

  lbd <- length(bw)

  #------------------------------------------------------------------#
  # Assign each patient to a cross validation group                  #
  #------------------------------------------------------------------#
  cvlabel <- sample(x = rep( c(1L:2L), length.out = nPatients ))

  #------------------------------------------------------------------#
  # initialize matrices for parameter estimates & standard deviations#
  #------------------------------------------------------------------#
  betaHat0 <- matrix(data = 0.0, nrow = lbd, ncol = nCov)

  #------------------------------------------------------------------#
  # initialize matrices for variance estimate                        #
  #------------------------------------------------------------------#
  hatV <- matrix(data = 0.0, nrow = lbd, ncol = nCov)

  guess0 <- NULL
  guess1 <- NULL
  guess2 <- NULL

  tempFunc <- function(x, 
                       data.y,
                       data.x,
                       bw,
                       kType,
                       lType,
                       nPatients,
                       xIs,
                       yIs,
                       tt,
                       guess0,
                       guess1,
                       guess2,
                       distanceFunction,
                       cvlabel){

    betaHat0 <- betaHat(data.y = data.y,
                        data.x = data.x,
                        bandwidth = bw[x],
                        kType = kType,
                        lType = lType,
                        nPatients = nPatients,
                        xIs = xIs,
                        yIs = yIs,
                        tt = time,
                        guess = guess0,
                        distanceFunction = distanceFunction)

    tst <- cvlabel == 1L

    betaHat1 <- betaHat(data.y = data.y,
                        data.x = data.x,
                        bandwidth = bw[x],
                        kType = kType,
                        lType = lType,
                        nPatients = sum(tst),
                        xIs = xIs,
                        yIs = yIs[tst],
                        tt = time,
                        guess = guess1,
                        distanceFunction = distanceFunction)

    tst <- cvlabel == 2L

    betaHat2 <- betaHat(data.y = data.y,
                        data.x = data.x,
                        bandwidth = bw[x],
                        kType = kType,
                        lType = lType,
                        nPatients = sum(tst),
                        xIs = xIs,
                        yIs = yIs[tst],
                        tt = time,
                        guess = guess2,
                        distanceFunction = distanceFunction)

    betaDiff <- betaHat2 - betaHat1

    hatV <- nPatients * bw[x] * ( betaDiff * betaDiff ) * 0.25


    return( list( "beta0" = betaHat0,
                  "beta1" = betaHat1,
                  "beta2" = betaHat2,
                  "hatV" = hatV) )

  }

  if( nCores > 1.1 ) {
    cl <- makeCluster(nCores)
    res <- parLapply(cl, 1L:lbd, tempFunc, data.y = data.y,
                                           data.x = data.x,
                                           bw = bw,
                                           kType = kType,
                                           lType = lType,
                                           nPatients = nPatients,
                                           xIs = xIs,
                                           yIs = yIs,
                                           tt = time,
                                           guess0 = NULL,
                                           guess1 = NULL,
                                           guess2 = NULL,
                                           distanceFunction = distanceFunction,
                                           cvlabel = cvlabel)

    stopCluster(cl)

    for( bd in 1L:lbd ) {

      betaHat0[bd,] <- res[[bd]]$beta0
      hatV[bd,] <- res[[bd]]$hatV

    }

  } else {

    for( bd in 1L:lbd ) {

      res <- tempFunc(data.y = data.y,
                      data.x = data.x,
                      bw = bw[bd],
                      kType = kType,
                      lType = lType,
                      nPatients = nPatients,
                      xIs = xIs,
                      yIs = yIs,
                      tt = time,
                      guess0 = guess0,
                      guess1 = guess1,
                      guess2 = guess2,
                      distanceFunction = distanceFunction,
                      cvlabel = cvlabel)

      betaHat0[bd,] <- res$beta0
      guess0 <- res$beta0
      guess1 <- res$beta1
      guess2 <- res$beta2
      hatV[bd,] <- res$hatV

    }
  }

  hatC <- array(data = 0.0, dim = nCov)

  for( p in 1L:nCov ) {
    hatC[p] <- lm( betaHat0[,p] ~ bw^2 )$coef[2]
  }

  #------------------------------------------------------------------#
  # MSE = hat(C) * hat(C) * bw^2 + hat(V)                            #
  #------------------------------------------------------------------#
  MSE <- bw^4 %o% ( hatC * hatC ) + hatV 

  #------------------------------------------------------------------#
  # Identify minimum MSE                                             #
  #------------------------------------------------------------------#
  mseFunc <- function(x) {

    tst <- x > 0.0
    if( sum(tst) == 0L ) {
      stop("no positive MSE values", call. = FALSE)
    }
    x[!tst] <- NA

    opt_h <- which.min(x)

    if( length(opt_h) > 1L ) {
      warning("Multiple minimums. Smallest bandwidth used.",
              call. = FALSE)
      opt_h <- opt_h[1]
    }
    return(opt_h)
  }

  opt_h <- apply(X = MSE, MARGIN = 2L, FUN = mseFunc)

  if( any(opt_h == 1L) || any(opt_h == lbd) ) {
    warning("At least 1 minimum is at bandwidth boundary.",
            call. = FALSE)
  }

  tminMSE <- array(data = 0.0, dim = nCov)
  for( i in 1L:nCov ) {
    tminMSE[i] <- MSE[opt_h[i],i]
  }

  optH <- bw[opt_h]

  minMSE <- tminMSE

  bHat <- betaHat(data.y = data.y,
                  data.x = data.x,
                  bandwidth = bw[opt_h],
                  kType = kType,
                  lType = lType,
                  tt = time,
                  guess = guess0,
                  xIs = xIs,
                  yIs = yIs,
                  nPatients = nPatients,
                  distanceFunction = distanceFunction)

  sdVec <- SD(data.y = data.y,
              data.x = data.x,
              bandwidth = bw[opt_h],
              kType = kType,
              lType = lType,
              bHat = bHat,
              tt = time,
              xIs = xIs,
              yIs = yIs,
              nPatients = nPatients,
              distanceFunction = distanceFunction)

  results <- matrix(data = 0.0,
                    nrow = nCov,
                    ncol = 6L,
                    dimnames = list(paste("beta",0L:{nCov-1L},sep=""),
                                    c("estimate","stdErr","z-value",
                                      "p-value", "optBW", "minMSE")))

  results[,1L] <- bHat
  results[,2L] <- sdVec
  results[,3L] <- bHat/sdVec
  results[,4L] <- 2.0*pnorm(-abs(results[,3L]))
  results[,5L] <- bw[opt_h]
  results[,6L] <- tminMSE

  print(results)
  cat("\n")

  return( list( "betaHat" = results[,1L],
                "stdErr"  = results[,2L],
                "zValue"  = results[,3L],
                "pValue"  = results[,4L],
                "optBW"   = results[,5L],
                "minMSE" =  results[,6L] ) )

}
