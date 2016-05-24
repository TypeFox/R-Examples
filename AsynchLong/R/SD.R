#----------------------------------------------------------------------#
# SD : Estimate the standard deviation for parameter estimates         #
#----------------------------------------------------------------------#
#                                                                      #
# data.y    : Matrix of responses. It is assumed that the first column #
#             contains integer patient IDs, the second column contains #
#             the time of measurement, and the third column contains   #
#             the value of the measurement.                            #
#                                                                      #
# data.x    : Matrix of covariates. It is assumed that the first column#
#             contains integer patient IDs, the second column contains #
#             the time of measurement, and the remaining columns       #
#             contain the values of the covariates.                    #
#                                                                      #
# bandwidth : a vector or numeric object of bandwidths                 #
#                                                                      #
# kType     : a character. One of "epan", "uniform", or "gauss".       #
#             Specifies the form of the kernel function.               #
#                                                                      #
# lType     : a character. One of "identity", "log", "logistic".       #
#             Specifies the form of the link function.                 #
#                                                                      #
# bHat      : an object of class numeric                               #
#             parameter estimates.                                     #
#                                                                      #
# xIs       : an object of class list.                                 #
#             v - which elements of x match y[i]'s id                  #
#             n - number of elements that match.                       #
#                                                                      #
# yIs       : an object of class list.                                 #
#             v - which elements of y match patientID[i]               #
#             n - number of elements that match.                       #
#                                                                      #
# nPatients : an object of class numeric.                              #
#             the number of patients in dataset.                       #
#                                                                      #
# distanceFunction : an object of class character.                     #
#             name of the distance function to use for calculation     #
#                                                                      #
# tt        : If provided, a vector of times at which to evaluate the  #
#             kernel                                                   #
#                                                                      #
#----------------------------------------------------------------------#
#                                                                      #
# Returns a vector of standard deviations.                             #
#                                                                      #
#----------------------------------------------------------------------#
SD <- function(bHat,
               data.y, 
               data.x,
               bandwidth, 
               kType,
               lType,
               nPatients,
               xIs,
               yIs,
               distanceFunction,
               tt = NULL) {

  #------------------------------------------------------------------#
  # Minimize u-function to estimate parameters                       #
  #------------------------------------------------------------------#
  argList <- list("xIs" = xIs,
                  "data.x" = data.x,
                  "data.y" = data.y,
                  "bandwidth" = bandwidth,
                  "kType" = kType,
                  "tt" = tt)

  dis <- do.call(what = distanceFunction, args = argList)

  xIs <- dis$xIs
  dis <- dis$dis

  res <- cmp_duFunc(pars = bHat,
                    data.y = data.y, 
                    data.x = data.matrix(data.x[, 3L:ncol(data.x), drop=FALSE]),
                    kernel = dis, 
                    lType = lType,
                    xIs = xIs,
                    yIs = yIs,
                    nPatients = nPatients)

  dU <- res$du
  var <- res$var

  var <- t(var) %*% var

  invdU <- try(solve(dU), silent = TRUE)

  if( class(invdU) == 'try-error' ) {
    cat("Unable to invert derivative of estimating equation.\n")
    stop(attr(invdU,"condition"), call. = FALSE)
  }

  var <- invdU %*% var %*% t(invdU)

  return(sqrt( diag(var) ))

}

