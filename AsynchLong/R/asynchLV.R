#----------------------------------------------------------------------#
# asynchLV : Last Value Carried Forward                                #
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
# lType     : a character. One of "identity", "log", "logistic".       #
#             Specifies the form of the link function.                 #
#                                                                      #
#----------------------------------------------------------------------#
asynchLV <- function(data.x, 
                     data.y, 
                     lType = "identity", ...) {

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

  result <- kernelFixed(data.y = data.y,
                        data.x = data.x,
                        bandwidth = 0.01,
                        kType = NULL,
                        lType = lType,
                        time = NULL,
                        distanceFunction = "distanceLV", ...)

  return( result )

}
