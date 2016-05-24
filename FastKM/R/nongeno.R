nongeno <- function(mat,
                    kernel = 'linear',
                    weights = NULL) {

  #--------------------------------------------------------------#
  # Internally, mat must be an object of class matrix.           #
  #--------------------------------------------------------------#
  if( is(mat, "data.frame") ) {
    mat <- data.matrix(mat)
  } else if( !is(mat, "matrix") ) {
    mat <- matrix(data = mat, ncol = 1L)
  }

  #--------------------------------------------------------------#
  # Shift to lower case to simply matching logic.                #
  #--------------------------------------------------------------#
  kernel <- tolower(kernel)

  #------------------------------------------------------------------#
  # Verify kernel.                                                   #
  #------------------------------------------------------------------#
  testKernel <- kernel %in% c("interactive", "linear", "quadratic")

  if( !testKernel ) {
    #--------------------------------------------------------------#
    # kernel must be one of interactive, linear, or quadratic.     #
    #--------------------------------------------------------------#
    stop(paste("For non-genotype data, kernel must be one of",
               " {interactive, linear, quadratic}", sep=""),
         call. = FALSE)
  }

  if( ncol(mat) == 1L ) {
    #--------------------------------------------------------------#
    # If only 1 variable is given, kernel must be linear.          #
    #--------------------------------------------------------------#
    if( kernel != "linear" ) {
      kernel <- "linear"
      cat("mat kernel changed to linear.\n", sep="")
    }
  }

  #------------------------------------------------------------------#
  #                          Verify weights                          #
  #------------------------------------------------------------------#
  if( {kernel == "interactive"} && !is(weights,"NULL") ) {
    cat("weights ignored for interactive kernel.\n")
    weights <- NULL
  }

  weights <- processWgts(mat = mat, weights = weights)

  obj <- methods::new(Class = "nongeno",
                      mat = mat,
                      kernel = kernel,
                      weights = weights)

  return(obj)

}
