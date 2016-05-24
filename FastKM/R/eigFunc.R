#----------------------------------------------------------------------#
# Given a matrix kmat and the number of variable to consider in matrix,#
# find the low-rank approximation of kmat.                             #
#----------------------------------------------------------------------#
 
eigFunc <- function(kmat, p){

  errorFlag <- FALSE

  #------------------------------------------------------------------#
  # the maximum number of non-zero eigenvalues is the rank, assumed  #
  # to be min(nrow(kmat),p)                                          #
  #------------------------------------------------------------------#
  if( nrow(kmat) <= p ) {
    #--------------------------------------------------------------#
    # if( x < p ) do full decomposition.                           #
    #--------------------------------------------------------------#
    egDecomp <- try(expr = eigen(x = kmat, 
                                 symmetric = TRUE, 
                                 only.values = TRUE),
                   silent = TRUE)

    if( is(egDecomp, "try-error") ) {
      stop("Unable to obtain eigenvalue decomposition of kernel matrix.",
           call. = FALSE)
    }

  } else {

    #--------------------------------------------------------------#
    # Calculate only the p largest eigenvalues.                    #
    #--------------------------------------------------------------#
    egDecomp <- try(expr = rARPACK::eigs_sym(A = kmat, 
                                             k = p, 
                                             opts = list(retvec=FALSE)),
                    silent = TRUE)

    if( is(egDecomp, "try-error") ) {

      egDecomp <- try(expr = eigen(x = kmat, 
                                   symmetric = TRUE, 
                                   only.values = TRUE),
                      silent = TRUE)

      if( is(egDecomp, "try-error") ) {
        stop("Unable to obtain eigenvalue decomposition of kernel matrix.",
             call. = FALSE)
      }

      errorFlag <- TRUE

    }

  }

  #------------------------------------------------------------------#
  # Remove all zero valued eigenvalues                               #
  #------------------------------------------------------------------#
  if( any(egDecomp$values < -1.5e-8) ) {
    stop("Negative eigenvalues encountered.", call.=FALSE)
  }

  keep <- egDecomp$values > 1.5e-8
  if( all(keep) ) {
    #--------------------------------------------------------------#
    # If no zeros identified, send back NULL indicating a larger p #
    # should be used.                                              #
    #--------------------------------------------------------------#
    if( p < ncol(kmat) ) return(NULL)
  }

  egDecomp$values <- egDecomp$values[keep]
  gSum <- cumsum(egDecomp$values)
  gSum <- gSum/gSum[length(gSum)]
  nEV <- length(egDecomp$values)

  #------------------------------------------------------------------#
  # Calculate the nEV largest eigenvectors                           #
  #------------------------------------------------------------------#
  if( nEV < nrow(kmat) && !errorFlag ) {
    Zg <- rARPACK::eigs_sym(A = kmat, k = nEV)
  } else {
    Zg <- eigen(x = kmat, symmetric = TRUE)
    Zg$values <- Zg$values[1:nEV]
    Zg$vectors <- Zg$vectors[,1:nEV]
  }

  return(list(   "Zg" = Zg, 
              "propV" = gSum))

}

