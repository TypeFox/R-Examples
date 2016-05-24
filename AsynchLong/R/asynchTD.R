#----------------------------------------------------------------------#
# asynchTD : Time-dependent coefficients                               #
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
# times     : A vector of time points.                                 #
#                                                                      #
# kType     : a character. One of "epan", "uniform", or "gauss".       #
#             Specifies the form of the kernel function.               #
#                                                                      #
# lType     : a character. One of "identity", "log", "logistic".       #
#             Specifies the form of the link function.                 #
#                                                                      #
# bw        : a numeric or NULL. Kernel bandwidth.                     #
#             If NULL, autotune is used to determine optimal bandwidth.#
#                                                                      #
# nCores    : a numeric. Number of cores to use for auto-tune          #
#             implementation                                           #
#                                                                      #
#----------------------------------------------------------------------#
asynchTD <- function(data.x, 
                     data.y,
                     times, 
                     kType = "epan", 
                     lType = "identity",
                     bw = NULL, 
                     nCores = 1, ...) {

  #------------------------------------------------------------------#
  # Process and verify input datasets                                #
  #------------------------------------------------------------------#
  data.x <- preprocessX(data.x = data.x)
  data.y <- preprocessY(data.y = data.y)

  #------------------------------------------------------------------#
  # Scale times.                                                     #
  #------------------------------------------------------------------#
  rge <- range(c(data.y[,2L],data.x[,2L]))
  data.y[,2L] <- (data.y[,2L] - rge[1L])/(rge[2L] - rge[1L])
  data.x[,2L] <- (data.x[,2L] - rge[1L])/(rge[2L] - rge[1L])

  nTimes <- length(times)
  #------------------------------------------------------------------#
  # Assumes patient id and measurement time columns                  #
  #------------------------------------------------------------------#
  nCov <- ncol(data.x) - 2L

  bHat <- matrix(data = 0.0, 
                 nrow = nTimes, 
                 ncol = nCov,
                 dimnames = list(NULL,colnames(data.x)[3L:(nCov+2L)]))

  sdVec <- bHat
  optH <- bHat
  minMSE <- bHat

  if( !is(bw,"NULL") ) {
    if( length(bw) > 1L ) {
      stop("Only a single bandwidth can be provided.")
    }
  }

  for( var in 1L:nTimes ) {

    cat("Time point: ", times[var], "\n")

    if( is(bw,"NULL") ) {

      result <- kernelAuto(data.x = data.x,
                           data.y = data.y,
                           kType = kType,
                           lType = lType,
                           time = times[var],
                           distanceFunction = "distanceTD",
                           nCores = nCores, ...)

      bHat[var,] <- result$betaHat
      sdVec[var,] <- result$stdErr
      optH[var,] <- result$optBW
      minMSE[var,] <- result$minMSE

    } else {

      result <- kernelFixed(data.y = data.y,
                            data.x = data.x,
                            bandwidth = bw,
                            kType = kType,
                            lType = lType,
                            time = times[var],
                            distanceFunction = "distanceTD", ...)

      bHat[var,] <- result$betaHat
      sdVec[var,] <- result$stdErr
    }
  }

  plotTD(bHat, sdVec, times)

  zv <- bHat/sdVec
  pv <- 2.0*pnorm(-abs(zv))

  if( is.null(bw) ) {
    return( list( "betaHat" = bHat,
                  "stdErr"  = sdVec,
                  "zValue"  = zv,
                  "pValue"  = pv,
                  "optBW"   = optH,
                  "minMSE"  = minMSE ) )
  } else {
    return( list( "betaHat" = bHat,
                  "stdErr"  = sdVec,
                  "zValue"  = zv,
                  "pValue"  = pv ) )

  }

}
