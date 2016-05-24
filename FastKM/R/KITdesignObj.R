#**********************************************************************#
#  KITdesignObj : A Fast Multiple-Kernel Method Based on a Low-Rank    #
#                 Approximation with Design Matrix Inputs defined by   #
#                 classes geno and nongeno                             #
#**********************************************************************#
#                                                                      #
#  Inputs :                                                            #
#                                                                      #
#           y : a vector, matrix, or data.frame of traits. Must be     #
#               a continuous, dichotomous, or survival trait. The code #
#               will deduce the trait type using the following logic:  #
#               If input argument delta is not NULL, assume survival.  #
#               If y is integer with only 2 unique values, dichotomous.#
#               If y is not integer or has more than two unique vales, #
#               continuous.                                            #
#               Missing values must be coded as NA.                    #
#                                                                      #
#        matA : an object of class geno or nongeno.                    #
#                                                                      #
#        matB : an object of class geno or nongeno.                    #
#                                                                      #
#        matC : an object of class geno, nongeno, or NULL.             #
#               If NULL, the product definition of Akernel x Bkernel   #
#               will be used to calculate the interaction kernel.      #
#                                                                      #
#           x : a vector, matrix, or data.frame. The design matrix of  #
#               covariates that are not included in either matA or matB#
#               Missing values should be coded as NA.                  #
#                                                                      #
#       delta : the status indicator in survival analyses.             #
#               Usually, 0=alive, 1=dead; TRUE/FALSE (TRUE=death); or  #
#               1/2 (2 = death).                                       #
#                                                                      #
# standardize : TRUE/FALSE indicating if non-genotype data should be   #
#               centered and scaled.                                   #
#                                                                      #
#    ellipsis : arguments to be passed to SKAT, survival, or coxKM     #
#                                                                      #
#  Outputs :                                                           #
#                                                                      #
# A list is returned containing:                                       #
#   probA : The proportion of variability maintained in matA           #
#   probB : The proportion of variability maintained in matB           #
#   pValues : A vector of p-values.                                    #
#                                                                      #
#**********************************************************************#
KITdesign <- function(y,
                      matA,
                      matB,
                      matC = NULL,
                      x = NULL, 
                      trait = NULL,
                      delta = NULL, 
                      standardize = TRUE, ...){

  #------------------------------------------------------------------#
  #                          PROCESS INPUTS                          #
  #------------------------------------------------------------------#
  givenC <- !is(matC, "NULL")
  if( !givenC ) matC <- methods::new("nongeno")

  givenX <- !is(x, "NULL")

  if( !is(matA, "geno") && !is(matA, "nongeno")) {
    stop("matA must be an object of class geno or nongeno.")
  }

  if( !is(matB, "geno") && !is(matB, "nongeno")) {
    stop("matB must be an object of class geno or nongeno.")
  }

  if( givenC && !is(matC, "nongeno")) {
    stop("matC must be an object of class nongeno.")
  }

  #------------------------------------------------------------------#
  # y and x must be matrices. If provided as data.frames or vectors  #
  # convert accordingly.                                             #
  #------------------------------------------------------------------#
  if( is(y, "data.frame") ) {
    y <- data.matrix(y)
  } else if( !is(y, "matrix") ) {
    y <- matrix(data = y, ncol = 1L)
  }

  if( givenX ) {
    if( is(x, "data.frame") ) {
      x <- data.matrix(x)
    } else if( !is(x, "matrix") ) {
      x <- matrix(data = x, ncol = 1L)
    }
  }

  #------------------------------------------------------------------#
  # Verify dimensions of matrices.                                   #
  #------------------------------------------------------------------#
  nSamples <- nrow(y)

  nms <- c("matA","matB","matC","x")
  tst <- c(nrow(matA@mat), 
           nrow(matB@mat), 
           min(nrow(matC@mat), nSamples), 
           min(nrow(x), nSamples))

  if( any(tst != nSamples) ) {
    stop(paste("Dimensions of ", nms[tst != nSamples],
               " do not agree with y.", sep=""), call. = FALSE)
  }

  #------------------------------------------------------------------#
  # Identify any incomplete cases (NA) and remove from all matrices  #
  #------------------------------------------------------------------#
  complete <- stats::complete.cases(y, matA@mat, matB@mat, matC@mat, x)
  if( any(!complete) ) {
    y <- y[complete,,drop=FALSE]
    matA@mat <- matA@mat[complete,,drop=FALSE]
    matB@mat <- matB@mat[complete,,drop=FALSE]
    if( givenC ) matC@mat <- matC@mat[complete,,drop=FALSE]
    if( givenX ) x <- x[complete,,drop=FALSE]
    #--------------------------------------------------------------#
    # Warn user that samples were removed.                         #
    #--------------------------------------------------------------#
    cat("Removed", sum(!complete), "incomplete case(s).\n", sep=" ")
  }

  #------------------------------------------------------------------#
  #                Standardize non-genotype matrices                 #
  #------------------------------------------------------------------#
  if( is(matA,"nongeno") && standardize ) {
    matA@mat <- scale(x = matA@mat, center = TRUE, scale = TRUE)
  }

  if( is(matB,"nongeno") && standardize  ) {
    matB@mat <- scale(x = matB@mat, center = TRUE, scale = TRUE)
  }

  if( givenC && is(matC,"nongeno") && standardize ) {
    matC@mat <- scale(x = matC@mat, center = TRUE, scale = TRUE)
  }

  if( givenX && standardize ) {
    x <- scale(x = x, center = TRUE, scale = TRUE)
  }

  #------------------------------------------------------------------#
  #                     Calculate kernel matrices                    #
  #------------------------------------------------------------------#
  aKernel <- calcKernObj(matA)
  bKernel <- calcKernObj(matB)

  if( !givenC ) {
    #--------------------------------------------------------------#
    # If matC not provided, calculate AxB kernel                   #
    #--------------------------------------------------------------#
    AhasOne <- 1.0*(matA@kernel %in% c("quadratic", "interactive"))
    BhasOne <- 1.0*(matB@kernel %in% c("quadratic", "interactive"))

    cKernel <- aKernel - AhasOne
    cKernel <- cKernel * (bKernel - BhasOne)

  } else {

    cKernel <- calcKernObj(matC)

  }

  if( any(is.nan(aKernel)) ) stop("Encountered NaN in A kernel.", 
                                  call. = FALSE)
  if( any(is.nan(bKernel)) ) stop("Encountered NaN in B kernel.",  
                                  call. = FALSE)
  if( any(is.nan(cKernel)) ) stop("Encountered NaN in C kernel.",  
                                  call. = FALSE)

  if( any(is.na(aKernel)) ) stop("Encountered NA in A kernel.", 
                                  call. = FALSE)
  if( any(is.na(bKernel)) ) stop("Encountered NA in B kernel.",  
                                  call. = FALSE)
  if( any(is.na(cKernel)) ) stop("Encountered NA in C kernel.",  
                                  call. = FALSE)

  res <- KITkernel(y = y,
                   kmatA = aKernel,
                   kmatB = bKernel,
                   kmatC = cKernel,
                   x = x,
                   AkernelC = 0.0,
                   BkernelC = 0.0,
                   delta = delta, 
                   trait = trait, 
                   standardize = FALSE, ...)


  return(res)


}
