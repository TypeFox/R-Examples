## decomposes an asymmetric proximity matrix into a symmetric and a skew-symmetric part (P = M + N)

symdecomp <- function(P) {
  ## P ... input proximity matrix
  P <- as.matrix(P)
  if (ncol(P) != nrow(P)) stop("Input matrix needs to be square!")
  M <- (P + t(P))/2
  N <- (P - t(P))/2
  return(list(M = M, N = N))
}

