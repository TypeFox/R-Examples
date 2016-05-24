#******************************************************************************#
# Newton Raphson algorithm for all methods.                                    #
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
#  betaGuess      an object of class numeric or NULL.                          #
#                 If numeric, beta will be initialized to the values           #
#                 provided.                                                    #
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
#  Returns an object of class numeric containing the parameter estimates.      #
#                                                                              #
#******************************************************************************#
betaEst <- function(Z,  
                    X,  
                    tau,  
                    h,
                    kType,
                    betaGuess,
                    tol,
                    maxiter,
                    scoreFunction){


  #--------------------------------------------------------------------------#
  # If a starting value for Newton-Raphson provided, use. Else, initialize   #
  # to small positive value (0.01).                                          #
  #--------------------------------------------------------------------------#
  if( is(betaGuess, "NULL") ) {
    beta <- numeric(ncol(Z) - 2L)
    beta[] <- 0.01
  } else {
    beta <- betaGuess
  }

  iter <- 0L

  while( TRUE ) {

    #----------------------------------------------------------------------#
    # Calculate Score Function                                             #
    #----------------------------------------------------------------------#
    argList <- list("beta" = beta, 
                    "Z" = Z, 
                    "X" = X, 
                    "tau" = tau, 
                    "h" = h, 
                    "kType" = kType)
    
    Lvec <- do.call(scoreFunction, argList)

    #----------------------------------------------------------------------#
    # Verify derivative can be inverted.                                   #
    #----------------------------------------------------------------------#
    Ldet <- det(Lvec$dUdBeta)
    if( (Ldet < 1.5e-8 && Ldet > -1.5e-8) ) {
      stop("Singular matrix encountered in Newton-Raphson.")
    }

    #----------------------------------------------------------------------#
    # Calculate new parameter values                                       #
    #----------------------------------------------------------------------#
    change <- solve(Lvec$dUdBeta) %*% Lvec$U

    beta.hat <- beta - change

    #----------------------------------------------------------------------#
    # Determine if parameter estimates have converged.                     #
    #----------------------------------------------------------------------#
    test <- TRUE
    for( i in 1L:length(beta) ) {
      if( abs(beta[i]) > 0.001 ) {
        if( abs(change[i])/abs(beta[i]) > tol ) test <- FALSE
      } else {
        if( abs(change[i]) > tol ) test <- FALSE
      }
    }

    if( test ) break

    beta <- beta.hat

    #----------------------------------------------------------------------#
    # Increment iterations and verify that maximum not yet reached         #
    #----------------------------------------------------------------------#
    iter <- iter + 1L
    if(iter >= maxiter) {
      warning(paste("Parameter estimates did not converge within ", 
                    maxiter, "iterations.", sep=""))
      break
    }
  }

  #--------------------------------------------------------------------------#
  # Return parameter estimates                                               #
  #--------------------------------------------------------------------------#
  return(beta)
}


