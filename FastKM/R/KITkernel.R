#**********************************************************************#
#  KITkernel : A Fast Multiple-Kernel Method Based on a Low-Rank       #
#              Approximation with kernel matrix inputs.                #
#**********************************************************************#
#                                                                      #
#  Inputs :                                                            #
#                                                                      #
#           y : a vector, matrix, or data.frame of traits. Must be     #
#               a continuous, dichotomous, or survival trait. The code #
#               will deduce the trait type using the following logic:  #
#               If input argument delta is not NULL, assumed survival. #
#               If y is integer with only 2 unique values, dichotomous.#
#               If y is not integer or has more than two unique vales, #
#               continuous.                                            #
#               Data must be complete.                                 #
#                                                                      #
#       kmatA : a matrix or data.frame. A kernel matrix.               #
#               Matrix must be complete.                               #
#                                                                      #
#       kmatB : a matrix or data.frame. A kernel matrix.               #
#               Matrix must be complete.                               #
#                                                                      #
#       kmatC : a matrix or data.frame. The kernel matrix for the      #
#               interaction. If provided, matrix must be complete. If  #
#               not provided, AxB kernel will be calculated from kmatA #
#               and kmatB.                                             #
#                                                                      #
#           x : a vector, matrix, or data.frame. The design matrix of  #
#               covariates that are not included in either kmatA or    #
#               kmatB. Data must be complete.                          #
#               Matrix x will be standardized.                         #
#                                                                      #
#    AkernelC : If kmatC is not provided and kmatA is polynomial or    #
#               interactive, the value of the constant term.           #
#               For example if kmatA = (1+X^T X), c = 1.0              #
#                                                                      #
#    BkernelC : If kmatC is not provided and kmatB is polynomial or    #
#               interactive, the value of the constant term.           #
#               For example if kmatB = X^T X, c = 0.0                  #
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
#   probA : The proportion of variability maintained in kmatA          #
#   probB : The proportion of variability maintained in kmatB          #
#   pValues : A vector of p-values.                                    #
#                                                                      #
#**********************************************************************#
KITkernel <- function(y,
                      kmatA,
                      kmatB,
                      kmatC = NULL,
                      x = NULL, 
                      AkernelC = 0.0,
                      BkernelC = 0.0,
                      trait = NULL, 
                      delta = NULL, 
                      standardize = TRUE, ...){

  #------------------------------------------------------------------#
  #                          PROCESS INPUTS                          #
  #------------------------------------------------------------------#
  # Identify if default values were used in formals.                 #
  #------------------------------------------------------------------#
  givenC <- !is(kmatC, "NULL")
  givenX <- !is(x, "NULL")
  givenDelta <- !is(delta, "NULL")

  #------------------------------------------------------------------#
  # y, kmatA, kmatB, kmatC, and x must be matrices.                  #
  #------------------------------------------------------------------#
  testMatrix <- function(mat) {
                  if( is(mat, "data.frame") ) {
                    mat <- data.matrix(mat)
                  } else if( !is(mat, "matrix") ) {
                    mat <- matrix(data = mat, ncol = 1L)
                  }
                  return(mat)
                }

  y <- testMatrix(y)
  kmatA <- testMatrix(kmatA)
  kmatB <- testMatrix(kmatB)
  if( givenC ) kmatC <- testMatrix(kmatC)
  if( givenX ) x <- testMatrix(x)

  #------------------------------------------------------------------#
  # Verify dimensions of matrices.                                   #
  #------------------------------------------------------------------#
  nSamples <- nrow(y)

  nms <- c("kmatA","kmatB","kmatC","x")
  tst1 <- c(nrow(kmatA), 
            nrow(kmatB), 
            min(nrow(kmatC), nSamples), 
            min(nrow(x), nSamples))
  tst2 <- c(ncol(kmatA), 
            ncol(kmatB), 
            min(ncol(kmatC), nSamples), 
            nSamples)

  if( any({tst1 != nSamples} | {tst2 != nSamples}) ) {
    stop(paste("Dimensions of ",
               nms[{tst1 != nSamples} | {tst2 != nSamples}],
               " do not agree with y.",sep=""), call. = FALSE)
  }

  #------------------------------------------------------------------#
  # Stop if any incomplete cases (NA) are identified                 #
  #------------------------------------------------------------------#
  complete <- stats::complete.cases(y, kmatA, kmatB, kmatC, x)
  if( any(!complete) ) stop("Data must be complete.", call. = FALSE)

  #------------------------------------------------------------------#
  #                           Standardize x                          #
  #------------------------------------------------------------------#
  if( givenX && standardize ) {
    x <- scale(x = x, center = TRUE, scale = TRUE)
  }

  #------------------------------------------------------------------#
  #                     Calculate kernel matrix                      #
  #------------------------------------------------------------------#
  if( !givenC ) {
    kmatC <- kmatA - AkernelC
    kmatC <- kmatC * (kmatB - BkernelC)
  }

  #------------------------------------------------------------------#
  #                      Obtain low-rank approximations              #
  #------------------------------------------------------------------#
  # Compute eigen decomposition of kmatA.                            #
  #------------------------------------------------------------------#
  dn <- max(floor(ncol(kmatA)*0.10),1L)
  pA <- 0L
  while( TRUE ) {
    pA <- min(pA + dn, ncol(kmatA))
    Za <- eigFunc(kmatA, pA)
    if( !is(Za, "NULL") ) break
  }

  nZa <- length(Za$Zg$values)
  aEV <- Za$propV[1L:nZa]
  Za <- Za$Zg$vectors %*% diag(x = sqrt(Za$Zg$values), 
                               nrow = nZa, 
                               ncol = nZa)

  #------------------------------------------------------------------#
  # Compute eigen decomposition of kmatB.                            #
  #------------------------------------------------------------------#
  dn <- max(floor(ncol(kmatB)*0.10),1L)
  pB <- 0L
  while( TRUE ) {
    pB <- min(pB + dn, ncol(kmatB))
    Zb <- eigFunc(kmatB, pB)
    if( !is(Zb, "NULL") ) break
  }
  nZb <- length(Zb$Zg$values)
  bEV <- Zb$propV[1L:nZb]
  Zb <- Zb$Zg$vectors %*% diag(x = sqrt(Zb$Zg$values), 
                               nrow = nZb,  
                               ncol = nZb)

  #------------------------------------------------------------------#
  #                 Create augmented covariate matrix                #
  #------------------------------------------------------------------#
  A <- cbind(x, Za, Zb)

  #------------------------------------------------------------------#
  #                    Deduce the type of trait.                     #
  #------------------------------------------------------------------#
  trait <- deduceTrait(givenDelta = givenDelta, 
                       y = y,  
                       trait = trait)

  #------------------------------------------------------------------#
  #                   Perform Kernel Machine Method                  #
  #------------------------------------------------------------------#
  if( trait == "c" || trait == "d" ) {
    #--------------------------------------------------------------#
    # Run KM analysis with SKAT.                                   #
    #--------------------------------------------------------------#
    result <- skatKM(A = A,
                     cKernel = kmatC,
                     y = y,
                     aEV = aEV,
                     bEV = bEV,
                     opt = toupper(trait), ...)

  } else if( trait == "s" ) {
    #--------------------------------------------------------------#
    # Run KM analysis with coxKM.                                  #
    #--------------------------------------------------------------#
    result <- survivalKM(A = A,
                         cKernel = kmatC,
                         y = y,
                         delta = delta,
                         aEV = aEV,
                         bEV = bEV, ...)
    
  }

  #------------------------------------------------------------------#
  #   Determine proportion of variability kept for kmatA and kmatB   #
  #------------------------------------------------------------------#
  if( !is(result$pv, "NULL") ) {

    #--------------------------------------------------------------#
    # If successful, retrieve proportions kept.                    #
    #--------------------------------------------------------------#
    propA <- aEV[result$nZa]
    propB <- bEV[result$nZb]

  } else {

    #--------------------------------------------------------------#
    # If unsuccessful, set proportions of a and b kept as NA.      #
    #--------------------------------------------------------------#
    propA <- NA
    propB <- NA

  }

  res <- list( "propA" = propA,
               "propB" = propB,
              "pValue" = result$pv)

  show(res)

  return(res)


}
