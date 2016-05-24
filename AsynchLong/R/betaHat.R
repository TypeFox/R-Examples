#----------------------------------------------------------------------#
# betaHat : uses R's optim to estimate parameters                      #
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
# nPatients : an object of class numeric.                              #
#             the number of patients in dataset.                       #
#                                                                      #
# xIs       : an object of class list.                                 #
#             v - which elements of x match y[i]'s id                  #
#             n - number of elements that match.                       #
#                                                                      #
# yIs       : an object of class list.                                 #
#             v - which elements of y match patientID[i]               #
#             n - number of elements that match.                       #
#                                                                      #
# distanceFunction : an object of class character.                     #
#             name of the distance function to use for calculation     #
#                                                                      #
# tt        : If provided, a vector of times at which to evaluate the  #
#             kernel                                                   #
#                                                                      #
#                                                                      #
# guess     : If provided, initial guess for beta                      #
#----------------------------------------------------------------------#
#                                                                      #
# Returns a vector of parameter estimates.                             #
#                                                                      #
#----------------------------------------------------------------------#
betaHat <- function(data.y, 
                    data.x,
		    bandwidth, 
                    kType,
                    lType,
                    nPatients,
                    xIs,
                    yIs,
                    distanceFunction, 
		    tt,
                    guess = NULL) {

  #------------------------------------------------------------------#
  # Determine the number of covariates & initialize parameter        #
  # estimates. Note that an intercept term is assumed here.          #
  #------------------------------------------------------------------#
  nCov <- ncol(data.x) - 2L

  if( is(guess, "NULL") ) {
    beta <- array(runif(nCov, min=0.0, max=0.5))
  } else {
    beta <- guess
  }

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

  if( lType == 'identity' ) {
    par <- uFuncIden(data.y = data.y,
                     data.x = data.matrix(data.x[,3L:ncol(data.x), drop=FALSE]),
                     kernel = dis,
                     xIs = xIs,
                     yIs = yIs,
                     nPatients = nPatients)
  } else {

    opt <- stats::optim(par = beta,
                        fn = optFunc,
                        gr = doptFunc,
                        method = "Nelder-Mead",
                        data.y = data.y,
                        data.x = data.matrix(data.x[,3L:ncol(data.x), drop=FALSE]),
                        kernel = dis,
                        lType = lType,
                        xIs = xIs,
                        yIs = yIs,
                        nPatients = nPatients)

    if(opt$convergence != 0) {
      warning(paste("optim did not converge. Code: ", 
                    opt$convergence,
                    "\nMessage: ", opt$message,
                    "\nvalue: ", opt$value), call. = FALSE)
    }

    par <- opt$par
  }

  return(par)

}

