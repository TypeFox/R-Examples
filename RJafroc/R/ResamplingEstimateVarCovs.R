#' @importFrom stats cov
ResamplingEstimateVarCovs <- function(resampleMatrix) {
  I <- dim(resampleMatrix)[1]
  J <- dim(resampleMatrix)[2]
  covariances <- array(dim = c(I, I, J, J))
  
  for (i in 1:I) {
    for (ip in 1:I) {
      for (j in 1:J) {
        for (jp in 1:J) {
          covariances[i, ip, j, jp] <- cov(resampleMatrix[i, j, ], resampleMatrix[ip, jp, ])
        }
      }
    }
  }
  
  ret <- VarCovs(covariances)
  return(list(var = ret$var, cov1 = ret$cov1, cov2 = ret$cov2, cov3 = ret$cov3))
}

