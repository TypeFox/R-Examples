# Constructs the overall matrix of the model.
# pars:      Paramter values to be used.
# sigmaType: Structure of the matrices in the model.
#             0 - Diagonal
#             1 - Exchangeable
#             2 - AR(1)
# kK:        Number of random effects. This is the dimension of the overall matrix model.
# kR:        Number of variance components.
# kLh:       Number of subvariance components in each variance component. Each subvariance component has a covariance matrix.
# kLhi:      Number of random effects in each subvariance component. This are the dimensions of the subvariance component matrices.

constructSigma <- function(pars, sigmaType, kK, kR, kLh, kLhi) {
  # We call ovSigma the overall covariance matrix
  # Dimension of overall sigma.
  ovSigma <- matrix(0, kK, kK)
  
  # counterPars has the position of the paramters that are in each sigma matrix.
  counterPars <- 1
  counterSubVar <- 1
  # For each variance component we need to paste its variance matrix into ovSigma.
  # sigmaDim has the dimension of the individual variance matrix per subvariance component.
  counterDim <- 1
  for (i in 1:kR) {
    # Sigma is a diagonal matrix (one parameter.)
    if (sigmaType[i] == 0) {
      par1 <- pars[counterPars]
      counterPars <- counterPars + 1
      # Number of subvariance matrices: kLh[i]. These share one paramter.
      for (j in 1:kLh[i]) {
        dim0 <- kLhi[counterSubVar] # Dimension of current subvariance matrix
        counterSubVar <- counterSubVar + 1 # Index of subvariance matrices
        tmp_mat <- par1 * diag(dim0)
        # print(counterDim:(counterDim + dim0 - 1))
        # print(tmp_mat)
        ovSigma[counterDim:(counterDim + dim0 - 1), counterDim:(counterDim + dim0 - 1)] <- tmp_mat
        counterDim <- counterDim + dim0
      }
    }
    
    # Exchangeable matrix. Two parameters, diagonal and off-diagonal.
    if (sigmaType[i] == 1) {
      # Number of subvariance matrices: kLh[i]. These share the two paramters but may be of different sizes.
      par1 <- pars[counterPars]
      par2 <- pars[counterPars + 1]
      counterPars <- counterPars + 2
      for (j in 1:kLh[i]) {
        dim0 <- kLhi[counterSubVar] # Dimension of current subvariance matrix
        counterSubVar <- counterSubVar + 1 # Index of subvariance matrices
        tmp_mat <- par1 * diag(dim0)
        tmp_mat[lower.tri(tmp_mat)] <- par2
        tmp_mat[upper.tri(tmp_mat)] <- par2
        # print(counterDim:(counterDim + dim0 - 1))
        # print(tmp_mat)
        ovSigma[counterDim:(counterDim + dim0 - 1), counterDim:(counterDim + dim0 - 1)] <- tmp_mat
        counterDim <- counterDim + dim0
      }
    }
    
    # AR(1) matrix. Two parameters, sigma^2 (1) and pho (2).
    if (sigmaType[i] == 2) {
      # Number of subvariance matrices: kLh[i]. These share the two paramters but may be of different sizes.
      par1 <- pars[counterPars]
      par2 <- pars[counterPars + 1]
      counterPars <- counterPars + 2
      for (j in 1:kLh[i]) {
        dim0 <- kLhi[counterSubVar] # Dimension of current subvariance matrix
        counterSubVar <- counterSubVar + 1 # Index of subvariance matrices
        d0 <- abs(outer(1:dim0, 1:dim0, "-"))
        tmp_mat <- par1 * par2^d0
        # print(counterDim:(counterDim + dim0 - 1))
        # print(tmp_mat)
        ovSigma[counterDim:(counterDim + dim0 - 1), counterDim:(counterDim + dim0 - 1)] <- tmp_mat
        counterDim <- counterDim + dim0
      }
    }
  }
  return(ovSigma)
}
