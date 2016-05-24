processWgts <- function(mat, weights){

  if( is(weights, "NULL") ) {

    #--------------------------------------------------------------#
    # If weights are not provided for mat set weights to 1.        #
    #--------------------------------------------------------------#
    weights <- matrix(data = 1.0, nrow = ncol(mat), ncol = 1L) 

  } else {

    #--------------------------------------------------------------#
    # If weights are provided for mat, verify dimensions.          #
    #--------------------------------------------------------------#
    if( is(weights, "matrix") ) {
      dimOK <- (nrow(weights) == ncol(weights)) &&
               (nrow(weights) == ncol(mat) )

      if( !dimOK ) {
        stop(paste("weights matrix provided for mat ",
                   "is not appropriately dimensioned.", sep=""),
             call. = FALSE)
      }
    } else {
      dimOK <- length(weights) == ncol(mat)

      if( !dimOK ) {
        stop(paste("weights vector provided for mat",
                   " is not appropriately dimensioned.", sep=""),
             call. = FALSE)
      }
      weights <- matrix(weights, ncol = 1L)
    }
  }

  return(weights)

}
