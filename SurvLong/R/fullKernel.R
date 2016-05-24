#******************************************************************************#
# Public calling function for full kernel method                               #
#******************************************************************************#
#                                                                              #
# Inputs                                                                       #
#                                                                              #
#  X              an object of class data.frame.                               #
#                 The structure of the data.frame must be                      #
#                 \{patient ID, event time, event indicator\}.                 #
#                 Patient IDs must be of class integer or be able to be        #
#                 coerced to class integer without loss of information.        #
#                 Missing values must be indicated as NA.                      #
#                                                                              #
#  Z              an object of class data.frame.                               #
#                 The structure of the data.frame must be                      #
#                 \{patient ID, time of measurement, measurement(s)\}.         #
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
#  bw             an object of class numeric or NULL.                          #
#                 If numeric, parameter estimates will be obtained at each     #
#                 value. If null, auto tune method will be used.               #
#                                                                              #
#  tol            an object of class numeric.                                  #
#                 maximum allowed change in parameter estimates, beyond which  #
#                 the parameter estimates are deemed to have converged.        #
#                                                                              #
#  maxiter        an object of class numeric.                                  #
#                 maximum number of iterations allowed to attain convergence   #
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
fullKernel <- function(X, 
                       Z, 
                       tau,
                       kType = "epan", 
                       bw = NULL,
                       tol = 0.001,
                       maxiter = 100){

  if( is(bw,"NULL") ) {
    result <- kernelAuto(X = X,
                         Z = Z,
                         tau = tau,
                         kType = kType,
                         tol = tol,
                         maxiter = maxiter,
                         scoreFunction = "scoreFull")
  } else {
    result <- kernelFixed(X = X,
                          Z = Z,
                          tau = tau,
                          bandwidth = bw,
                          kType = kType,
                          tol = tol,
                          maxiter = maxiter,
                          scoreFunction = "scoreFull")
  }


  return(result)

}
