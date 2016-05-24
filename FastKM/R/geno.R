#----------------------------------------------------------------------#
# Define object of class geno.                                         #
#----------------------------------------------------------------------#
geno <- function(mat,
                 kernel = 'linear',
                 weights = NULL,
                 inheritMode = NA) {

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
  inheritMode <- tolower(inheritMode)

  #--------------------------------------------------------------#
  # Send mat to testGeno to verify kernel and inheritance        #
  #--------------------------------------------------------------#
  res <- testGeno(mat = mat, 
                  inherit = inheritMode, 
                  kernel = kernel)

  kernel <- res$kernel
  inheritMode <- res$inherit

  #------------------------------------------------------------#
  # Identify 9's, which indicate missing data.                 #
  #------------------------------------------------------------#
  incomplete <- apply(X = mat, 
                      MARGIN = 1, 
                      FUN = function(x){any( x > 8.5 )})

  if( any(incomplete) ) {
    #----------------------------------------------------------#
    # Change incomplete data to NA.                            #
    #----------------------------------------------------------#
    mat[incomplete,] <- NA
  }

  #------------------------------------------------------------------#
  #                          Verify weights                          #
  #------------------------------------------------------------------#
  matA <- toSNP(mat = mat, 
                inherit = inheritMode, 
                weights = weights)

  weights <- matA$weights
  mat <- matA$snp

  if( {ncol(weights) != 1L} && {kernel == 'ibs'} ) {
    stop("Must use vector weights for ibs", call. = FALSE)
  }

  obj <- methods::new(Class = "geno",
                      mat = mat,
                      kernel = kernel,
                      weights = weights)

  return(obj)

}
