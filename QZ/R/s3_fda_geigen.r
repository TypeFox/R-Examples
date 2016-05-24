### This function is equivalent to fda::geigen(Amat, Bmat, Cmat)

fda.geigen <- function(Amat, Bmat, Cmat){
  #  solve the generalized eigenanalysis problem
  #
  #    max {tr L'AM / sqrt[tr L'BL tr M'CM] w.r.t. L and M
  #
  #  Arguments:
  #  Amat ... p by q matrix
  #  Bmat ... order p symmetric positive definite matrix
  #  Cmat ... order q symmetric positive definite matrix
  #  Returns:
  #  VALUES ... vector of length s = min(p,q) of eigenvalues
  #  Lmat   ... p by s matrix L
  #  Mmat   ... q by s matrix M

  ### Check dimensions.
  Adim <- dim(Amat)
  Bdim <- dim(Bmat)
  Cdim <- dim(Cmat)
  if(Adim[1] != Bdim[1]){
    stop('dim(Amat)[1] != dim(Bmat)[1]')
  }
  if(Adim[2] != Cdim[1]){
    stop('dim(Amat)[2] != dim(Cmat)[1]')
  }
  if(Bdim[1] != Bdim[2]){
    stop('Bmat is not square')
  }
  if(Cdim[1] != Cdim[2]){
    stop('Cmat is not square')
  }

  ### Check symmetric.
  if(! isSymmetric.matrix(Bmat)){
     stop('Bmat is not symmetric')
  }
  if(! isSymmetric.matrix(Cmat)){
    stop('Cmat is not symmetric')
  }
  Bmat <- (Bmat + t(Bmat)) / 2
  Cmat <- (Cmat + t(Cmat)) / 2

  ### Rebuild (A,B).
  p <- Bdim[1]
  q <- Cdim[1]
  s <- min(c(p, q))
  A <- rbind(cbind(matrix(0, nrow = p, ncol = p), Amat),
             cbind(t(Amat), matrix(0, nrow = q, ncol = q)))
  B <- rbind(cbind(Bmat, matrix(0, nrow = p, ncol = q)),
             cbind(matrix(0, nrow = q, ncol = p), Cmat))
  ret.qz <- qz.dggev(A, B)

  ### Summary and scaling.
  values <- ret.qz$ALPHA / ret.qz$BETA
  id <- order(values, decreasing = TRUE)[1:s]
  Lmat <- ret.qz$U[1:p, id]
  Mmat <- ret.qz$U[-(1:p), id]
  scale.L <- sqrt(diag(t(Lmat) %*% Bmat %*% Lmat))
  scale.M <- sqrt(diag(t(Mmat) %*% Cmat %*% Mmat))
  Lmat <- sweep(Lmat, 2, scale.L, "/")
  Mmat <- sweep(Mmat, 2, scale.M, "/")

  ### For returns
  ret <- list(values = values[id], Lmat = Lmat, Mmat = Mmat)
  ret
} # End of fda.geigen().

