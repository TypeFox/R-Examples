ramsvm <- function(x, 
                   y, 
                   lambda, 
                   gamma = 0.5, 
                   weight = NULL, 
                   kernel = "linear", 
                   kparam = NULL, 
                   large = FALSE, 
                   epsilon = NULL, 
                   warm = NULL, 
                   nb.core = NULL) {

  #------------------------------------------------------------------#
  # Verify that kernel is one of linear, polynomial, gaussian        #
  #------------------------------------------------------------------#
  kernel <- tolower(kernel)
  if( !(kernel %in% c("linear","polynomial", "gaussian")) ) {
    stop("kernel must be one of ('linear', 'polynomial', 'gaussian')")
  }

  if( kernel %in% c("polynomial", "gaussian") ) {
    if( is(kparam, "NULL") ) kparam = 1.0
  }

  #------------------------------------------------------------------#
  # Verify that covariates are provided in matrix form with no NA/NAN#
  #------------------------------------------------------------------#
  if( !is(x, "matrix") ) stop("The covariates should be a matrix.")  
  if( any(is.na(x)) || any(is.nan(x)) ) {
    stop("There should be no NA/NaN in the covariates.")
  }

  nobs <- nrow(x)

  #------------------------------------------------------------------#
  # Verify that the number of labels matches the number of obs.      #
  #------------------------------------------------------------------#
  if( length(y) != nobs ) {
    stop(paste("The dimension of covariates should match the",
               "length of the label."))
  }

  #------------------------------------------------------------------#
  # determine the number of classes; must be > 1.                    #
  #------------------------------------------------------------------#
  k <- length(levels(as.factor(y)))
  if( k < 2L) stop("There must be at least two classes.")
  
  if( is(epsilon, "NULL") ) {
    #--------------------------------------------------------------#
    # If not provided, calculate epsilon.                          #
    #--------------------------------------------------------------#
    epsilon <- 0.0001 * nobs * k
  } else {
    #--------------------------------------------------------------#
    # If provided, verify that epsilon is a positive scalar.       #
    #--------------------------------------------------------------#
    if( is.na(epsilon) || is.nan(epsilon) ) {
      stop("epsilon should not be NA/NaN.")
    }
    if( !is(epsilon, "numeric") ) stop("Epsilon should be numeric.")
    if( epsilon <= 0.0 ) stop("Epsilon should be strictly positive.")
    if( length(epsilon) > 1.5 ) stop("Epsilon should be a scalar.")
  }
    
  if( is(weight, "NULL") ) {
    #--------------------------------------------------------------#
    # If not provided, set weight to 1.0                           #
    #--------------------------------------------------------------#
    weight <- numeric(nobs) + 1.0
  } else {
    #--------------------------------------------------------------#
    # If provided, verify that weight is a positive, numeric vector#
    # of length equivalent to the number of observations.          #
    #--------------------------------------------------------------#
    if( !is(weight, "numeric") ) {
      stop("The weight vector should be numeric.")
    }
    if( any(is.na(weight)) || any(is.nan(weight)) ) {
      stop("There should be no NA/NaN in the weight vector.")
    }
    if( min(weight) < 0.0 ) {
      stop("The weight vector should be non-negative.")
    }
    if( length(weight) != nobs ) {
      stop(paste("The length of the weight should agree with",
                 "the number of observations."))
    }
  }

  if( !is(warm, "NULL") ) {
    #--------------------------------------------------------------#
    # If provided, verify dims of starter matrix are appropriate.  #
    #--------------------------------------------------------------#
    if( {nrow(warm) != nobs} || {ncol(warm) != k} ) {
      stop("A dimension of the warm start matrix is incorrect.")
    }
  } else {
    #--------------------------------------------------------------#
    # If not provided, default starter matrix to zero              #
    #--------------------------------------------------------------#
    warm <- matrix(data = 0.0, nrow = nobs, ncol = k)
  }

  #------------------------------------------------------------------#
  # Verify that all lambdas are positive.                            #
  #------------------------------------------------------------------#
  if( any(lambda <= 0) ) stop("All lambdas must be positive.")

  #------------------------------------------------------------------#
  # If the lambdas are not sorted warn user that they are reorder.   #
  #------------------------------------------------------------------#
  if( any( order(lambda) != order(sort(lambda,TRUE))) ) {
    warning("The order of lambda has been changed.")
  }

  #------------------------------------------------------------------#
  # Sort lambda high to low.                                         #
  #------------------------------------------------------------------#
  lambda = as.double(rev(sort(lambda)))

  if( is(nb.core, "NULL") && large ) {
    #--------------------------------------------------------------#
    # If this is a large calculation and the number of cores are   #
    # not specified, detect the number of available cores.         #
    #--------------------------------------------------------------#
    nb.core <- parallel::detectCores()
  } else if( !large ) {
    nb.core <- 0L
  }

  #------------------------------------------------------------------#
  # Call solver method.                                              #
  #------------------------------------------------------------------#
  fit <- RAMSVM_solve(x = x,
                      y = y,
                      gamma = gamma,
                      lambda = lambda,
                      kernel = kernel,
                      kparam = kparam,
                      weight = weight,
                      epsilon = epsilon,
                      warm = warm,
                      nb.core = nb.core)

  fit@call <- match.call()
  return(fit)
}
