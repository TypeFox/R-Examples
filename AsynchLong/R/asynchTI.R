#----------------------------------------------------------------------#
# asynchTI : Time-invariant coefficients                               #
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
# bw        : a vector or numeric object of bandwidths                 #
#                                                                      #
# nCores    : a numeric. Number of cores to use for auto-tune          #
#             implementation                                           #
#                                                                      #
#----------------------------------------------------------------------#
asynchTI <- function(data.x, 
                     data.y,
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
  # Scale times                                                      #
  #------------------------------------------------------------------#
  rge <- range(c(data.y[,2L],data.x[,2L]))
  data.y[,2L] <- (data.y[,2L] - rge[1L])/(rge[2L] - rge[1L])
  data.x[,2L] <- (data.x[,2L] - rge[1L])/(rge[2L] - rge[1L])

  if( is(bw, "NULL") ) {

    result <- kernelAuto(data.x = data.x,
                         data.y = data.y,
                         kType = kType,
                         lType = lType,
                         time = NULL,
                         distanceFunction = "distanceTI", 
                         nCores = nCores, ...)

  } else {

    result <- kernelFixed(data.y = data.y,
                          data.x = data.x,
                          bandwidth = bw,
                          kType = kType,
                          lType = lType,
                          time = NULL,
                          distanceFunction = "distanceTI", ...)

  }

  return( result )

}
