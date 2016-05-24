#----------------------------------------------------------------------#
# kernelFixed : fixed bandwidths for all methods                       #
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
# bandwidth : vector of bandwidth values                               #
#                                                                      #
# kType     : a character. One of "epan", "uniform", or "gauss".       #
#             Specifies the form of the kernel function.               #
#                                                                      #
# lType     : a character. One of "identity", "log", "logistic".       #
#             Specifies the form of the link function.                 #
#                                                                      #
# time      : a numeric or vector object containing time points        #
#                                                                      #
# distanceFunction : a character object giving the distance function   #
#             to be used.                                              #
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
#                                                                      #
#----------------------------------------------------------------------#
kernelFixed <- function(data.x, 
                        data.y,
                        bandwidth,
                        kType, 
                        lType,
                        time,
                        distanceFunction, ...){

  #------------------------------------------------------------------#
  # bandwidths must be positive.                                     #
  #------------------------------------------------------------------#
  if( any(bandwidth < 1.5e-8) ) {
    stop("All bandwidths must be > 0", call. = FALSE)
  }

  #------------------------------------------------------------------#
  # nCov assumes a column for patient ids and measurement times      #
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
    xIs[[i]]$n <- length(xIs[[i]]$v) 
  }

  #------------------------------------------------------------------#
  # For each response, determine which response measurements are for #
  # the same sample and the number of measurements.                  #
  #------------------------------------------------------------------#
  yIs <- list()
  for( i in 1L:nPatients ) {
    yIs[[i]] <- list()
    yIs[[i]]$v <- which( data.y[,1L] == patientIDs[i] )
    yIs[[i]]$n <- length(yIs[[i]]$v)
  }

  lbd <- length(bandwidth)

  bHat <- matrix(data = 0.0, 
                 nrow = lbd, 
                 ncol = nCov,
                 dimnames = list(NULL,colnames(data.x)[3L:(nCov+2L)]))

  sdVec <- bHat

  results <- matrix(data = 0.0,
                    nrow = nCov,
                    ncol = 4L,
                    dimnames = list(paste("beta",0L:{nCov-1L},sep=""),
                                    c("estimate","stdErr","z-value",
                                      "p-value")))

  guess <- NULL

  for( bd in 1L:lbd ) {

    if(distanceFunction != "distanceLV") {
      cat("Bandwidth: ", bandwidth[bd], "\n")
    }

    bHat[bd,] <- betaHat(data.y = data.y,
                         data.x = data.x,
                         bandwidth = bandwidth[bd],
                         kType = kType,
                         lType = lType,
                         tt = time,
                         xIs = xIs,
                         yIs = yIs,
                         nPatients = nPatients,
                         guess = guess,
                         distanceFunction = distanceFunction)

    results[,1L] <- bHat[bd,]

    guess <- bHat[bd,]

    sdVec[bd,] <- SD(data.y = data.y,
                     data.x = data.x,
                     bandwidth = bandwidth[bd],
                     kType = kType,
                     lType = lType,
                     bHat = bHat[bd,],
                     tt = time,
                     xIs = xIs,
                     yIs = yIs,
                     nPatients = nPatients,
                     distanceFunction = distanceFunction)

    results[,2L] <- sdVec[bd,]
    results[,3L] <- bHat[bd,]/sdVec[bd,]
    results[,4L] <- 2.0*stats::pnorm(-abs(results[,3L]))

    print(results)
    cat("\n")
  }

  zv <- bHat/sdVec
  pv <- 2.0*pnorm(-abs(results[,3L]))

  return( list( "betaHat" = bHat,
                "stdErr"  = sdVec,
                "zValue"  = zv,
                "pValue"  = pv ) )

}
