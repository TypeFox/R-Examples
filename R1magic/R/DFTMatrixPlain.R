DFTMatrixPlain <- function(N) {
  # Part of R1Magic by mehmet.suzen@physics.org
  omega <- exp(-2*pi*1i/N)
  sq <- rep(0:(N-1),N)*rep(0:(N-1),each=N)
  AA <- sapply(sq,oo,omega) 
  A <- matrix(AA,N,N)
 return(A);
}

