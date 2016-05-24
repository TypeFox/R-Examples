cmat.star.BinOrdNN <-
function(plist, skew.vec, kurto.vec, no.bin, no.ord, no.NN, CorrMat) 
{
  no.binord <- no.bin + no.ord
  
  if (no.NN == 0) {
    Sigma <- IntermediateOO(plist, CorrMat)
  }
  if (no.binord == 0) {
    Sigma <- IntermediateNonNor(skew.vec, kurto.vec, CorrMat)
  }
  
  if (no.NN > 0 & no.binord > 0) {
    if (validate.target.cormat.BinOrdNN(plist, skew.vec, kurto.vec, no.bin, no.ord, no.NN, CorrMat)) {
      k <- length(plist)
      tot <- nrow(CorrMat)
      if (no.binord > 1) OO <- IntermediateOO(plist, CorrMat[1:k, 1:k])
      else OO <- 1
      ONN <- IntermediateONN(plist, skew.vec, kurto.vec, CorrMat[(k + 1):tot,1:k])
      if (no.NN > 1) NN <- IntermediateNonNor(skew.vec, kurto.vec, CorrMat[(k + 1):tot, (k + 1):tot])
      else NN <- 1
      Sigma <- cbind(rbind(OO, ONN), rbind(t(ONN), NN))
      if (!is.positive.definite(Sigma)) {
        warning("Intermediate correlation matrix is not positive definite. A nearPD function is applied.")
        Sigma <- as.matrix(nearPD(Sigma, corr = TRUE, keepDiag = TRUE)$mat)
      }
      Sigma <- (Sigma + t(Sigma))/2
    }
  }
  return(Sigma)
}
