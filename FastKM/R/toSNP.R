#----------------------------------------------------------------------#
# If appropriate, convert two-column-per-marker genotype data to       #
# one-column-per-marker using specified inheritance mode.              #
# Determine weights if no weights are provided by the user.            #
#----------------------------------------------------------------------#
toSNP <- function(mat, 
                  inherit,  
                  weights){

  nSamples <- nrow(mat)

  #------------------------------------------------------------------#
  # inherit == NA indicates one-column-per-marker format.            #
  #------------------------------------------------------------------#
  if( !is.na(inherit) ) {
    #--------------------------------------------------------------#
    # For two-column-per-marker format, use additive mode of       #
    # inheritance to obtain allele frequencies.                    #
    #--------------------------------------------------------------#
    odds <- seq(from = 1L, to = ncol(mat), by = 2L)
    evens <- odds + 1L

    snp <- mat[,odds,drop=FALSE] + mat[,evens,drop=FALSE]

    #--------------------------------------------------------------#
    # Calculate mean value of the columns.                         #
    #--------------------------------------------------------------#
    alleleFrequency <- colSums(snp)/(2.0*nSamples)

    #--------------------------------------------------------------#
    # Verify that minor allele frequencies are given.              #
    #--------------------------------------------------------------#
    if( any(alleleFrequency > 0.5) ) {
      #----------------------------------------------------------#
      # If some means are > 0.5, convert to minor allele freq.   #
      #----------------------------------------------------------#
      mat <- 1.0 - mat
 
      #----------------------------------------------------------#
      # Recalculate mean value of the columns.                   #
      #----------------------------------------------------------#
      alleleFrequency <- 1.0 - alleleFrequency

      if( any(alleleFrequency > 0.5) ) {

        #------------------------------------------------------#
        # If some means still > 0.5, warn user. Revert data.   #
        #------------------------------------------------------#
        warning("Encountered allele frequency > 0.5. \n", call.=FALSE)
        mat <- 1.0 - mat
        alleleFrequency <- 1.0 - alleleFrequency

      } else {

        #------------------------------------------------------#
        # Warn user that data interpreted to be major allele   #
        # frequencies & converted to minor allele frequencies. #
        #------------------------------------------------------#
        cat("Genotype data identified to be coded for major allele ",
            "frequencies. ", 
            "Converted to minor allele frequency (1-mat).\n",
            sep="")
      }
    }

    #--------------------------------------------------------------#
    # Combine genotype data to obtain one-column-per-marker        #
    # format using inheritance mode input.                         #
    #--------------------------------------------------------------#
    if( inherit == "add" ) {
      snp <- mat[,odds,drop=FALSE] + mat[,evens,drop=FALSE]
    } else if( inherit == "rec" ) {
      snp <- 1.0*({mat[,odds] >= 0.5} & {mat[,evens] >= 0.5})
    } else if( inherit == "dom" ) {
      snp <- 1.0*({mat[,odds] >= 0.5} | {mat[,evens] >= 0.5})
    } else {
      stop("Error: Mode of inheritance must be one of: add, dom, rec.",
           call. = FALSE)
    }

  } else {
    #--------------------------------------------------------------#
    # For one-column-per-marker, do nothing to matrix.             #
    #--------------------------------------------------------------#
    snp <- mat

    #--------------------------------------------------------------#
    # Calculate mean value of the columns.                         #
    #--------------------------------------------------------------#
    alleleFrequency <- colSums(snp)/(2.0*nSamples)

    if( any(alleleFrequency > 0.5) ) {
      #----------------------------------------------------------#
      # If some means are > 0.5, warn user.                      #
      #----------------------------------------------------------#
      warning("Some allele frequencies > 0.5. Verify data.",
              call. = FALSE)
    }
  }

  if( is(weights, "NULL") ) {
    #--------------------------------------------------------------#
    # If user did not provide a weight vector, use minor allele    #
    # frequencies to define default.                               #
    #--------------------------------------------------------------#
    weights <- matrix((1.0 - alleleFrequency)^{24}, ncol = 1L)

    if( is.na(inherit) ) {
      #----------------------------------------------------------#
      # If data in one-column-per-marker format, check for 2's.  #
      #----------------------------------------------------------#
      if( !any(snp > 1.5) ) {
        #------------------------------------------------------#
        # If no 2's, warn.                                     #
        #------------------------------------------------------#
        warning(paste("Additive model assumed in obtaining weights. ",
                      "Note only 0/1 values in genotype matrix.",
                       sep=""),
                call. = FALSE)
      }
    }
  } else {
    #--------------------------------------------------------------#
    # If user provided weights, verify number provided.            #
    #--------------------------------------------------------------#
    if( is(weights, "vector") ) {

      weights <- matrix(weights, ncol = 1L)

      numOK <- {nrow(weights) == ncol(snp)}

    } else {

      numOK <- {ncol(weights) == ncol(snp)} && 
               {nrow(weights) == ncol(snp)}

    }

    if( !numOK ) {
      stop(paste("Dimension(s) of weights ", 
                 "must be equivalent to the # of loci.", sep=""),
           call. = FALSE)
    }
  }

  return(list(    "snp" = snp,
              "weights" = weights))
}
