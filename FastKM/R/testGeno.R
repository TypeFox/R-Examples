#----------------------------------------------------------------------#
# If a matrix is indicated as genotype, determine if data are integer  #
# and verify that kernel and inherit are appropriately defined.        #
#----------------------------------------------------------------------#
testGeno <- function(mat, inherit, kernel) {

  #------------------------------------------------------------------#
  # kernel must be one of ibs, linear, or quadratic                  #
  #------------------------------------------------------------------#
  testKernel <- kernel %in% c("ibs", "linear", "quadratic")

  if( !testKernel ) {
    stop(paste("For genotype data, ",
               "kernel must be one of {ibs, linear, quadratic}",
               sep = ""), call. = FALSE)
  }

  #------------------------------------------------------------------#
  # Test to determine if all elements are integer.                   #
  #------------------------------------------------------------------#
  isInteger <- isTRUE(all.equal(round(mat,0),mat))

  if( !isInteger ) {

    cat("Data identified to be non-integer.\n", sep="")

    #--------------------------------------------------------------#
    # For non-integer data, kernel must be linear/quad             #
    #--------------------------------------------------------------#
    if( kernel == "ibs" ) {
      cat("     kernel reset to 'linear.'\n", sep="")
      kernel <- "linear"
    } 

    #--------------------------------------------------------------#
    # For non-integer data, mode of inheritance must be add or NA  #
    #--------------------------------------------------------------#
    if( !is.na(inherit) && {inherit != "add"} ) {
      cat("     inheritMode reset to 'add.'\n", sep="")
      inherit <- "add"
    }
  }

  if( is.na(inherit) ) {

    cat("Data in one-column-per-marker format.\n", sep="")

  } else {
    #--------------------------------------------------------------#
    # For two-column-per-marker data, verify even # of columns.    #
    #--------------------------------------------------------------#
    if( {ncol(mat) %% 2} > 0.5 ) {
      stop(paste("InheritMode != NA indicates ",
                 "two-column-per-marker format. ",
                 "Matrix must have an even ",
                 "number of columns.",sep=""), call. = FALSE)
    }
  }

  isOne <- {ncol(mat) == 2L && !is.na(inherit)} ||
           {ncol(mat) == 1L &&  is.na(inherit)}

  if( isOne ) {
    #--------------------------------------------------------------#
    # If only 1 loci given, kernel must be ibs or linear based on  #
    # type of data (integer/non-integer).                          #
    #--------------------------------------------------------------#

    if( isInteger && {kernel != "ibs"}) {
      kernel <- "ibs"
      cat("kernel changed to ibs.\n", sep="")
    } else if( !isInteger && {kernel != "linear"}) {
      kernel <- "linear"
      cat("kernel changed to linear.\n", sep="")
    }
  }

  return(list( "kernel" = kernel,
              "inherit" = inherit))

}

