normalizeAdj <-
function(Amat, alpha = 1) {
  if (is.null(dim(Amat))) {
    ncond = length(Amat)
    p = nrow(Amat[[1]])
  } else {
    ncond <- 1
    p <- dim(Amat)[1]
    Amat <- list(Amat)
  }
  
  LapMat = Amat 
  Lmat = Amat
  normA = Amat
  InfMat = Amat
  
  Ip = diag(rep(1, p), p, p)
  
  for (i in 1:ncond) {
    LapResults = graphlaplacian(Amat[[i]])
    LapMat[[i]] = LapResults$Lnorm
    Lmat[[i]] = LapResults$Lunnorm
    normA[[i]] = alpha * LapMat[[i]]
    InfMat[[i]] = t(chol(pseudoinverse(Ip - normA[[i]]))) 
  }
  
  return(list(normA = normA, Lmat = Lmat, LapMat = LapMat, InfMat = InfMat))
}
