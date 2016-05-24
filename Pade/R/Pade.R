Pade <- function(L, M, A){
  if (floor(L) != L || floor(M) != M) stop ("Polynomial orders need to be integers.")
  Lp1 <- L + 1
  MatSize <- Lp1 + M
  if(length(A) < MatSize) stop ("Not enough Taylor series coefficients provided.")
  PQ <- matrix(0, ncol = MatSize, nrow = MatSize)
  PQ[1:Lp1, 1:Lp1] <- -diag(Lp1)
  for (i in seq_len(M)){
    PQ[, Lp1 + i] <- c(rep.int(0, i), head(A, (MatSize - i)))
  }
  PadeCoeff <- solve(PQ, -head(A, MatSize))
  Numer <- head(PadeCoeff, Lp1)
  Denom <- c(1, tail(PadeCoeff, M))
  return(list(Px = Numer, Qx = Denom))
}
