validate.target.cormat.BinOrdNN <-
function (plist, skew.vec, kurto.vec, no.bin, no.ord, no.NN, CorrMat) 
{
  if (is.positive.definite(CorrMat) == FALSE) {
    stop("Specified correlation matrix is not positive definite! \n")
  }
  if (isSymmetric(CorrMat) == FALSE) {
    stop("Specified correlation matrix is not symmetric! \n")
  }
  if (ncol(CorrMat) != (no.bin + no.ord + no.NN)) {
    stop("Dimension of correlation matrix does not match the number of variables!\n")
  }
  if (sum(diag(CorrMat)) != ncol(CorrMat)) {
    stop("Diagonal elements of correlation matrix do not all equals to 1!\n")
  }
  Limits <- valid.limits.BinOrdNN(plist, skew.vec, kurto.vec, no.bin, no.ord, no.NN)
  minmat <- Limits$lower
  maxmat <- Limits$upper
  rangemat <- (minmat <= CorrMat & CorrMat <= maxmat)
  if (sum(!rangemat) > 0) {
    d <- ncol(rangemat)
    cat("Target matrix is not valid. The following values are invalid.\n")
    for (i in 1:d) {
      for (j in i:d) {
        if (rangemat[i, j] == FALSE) {
          cat("CorrMat[", i, ",", j, "] must be between", round(minmat[i, j], 3), "and", round(maxmat[i, j], 3), "\n")
        }
      }
    }
    stop("Range violation occurred in the target correlation matrix.\n")
  }
  return(TRUE)
}
