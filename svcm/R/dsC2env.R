dsC2env <- function(A, Ccode = TRUE) {
  ##A: dsCMatrix in upper triangle form (requested by Ccode)
  
  N <- nrow(A)
  ##calculate lower bandwidth of each row
  if (A@uplo == "U") {
    bw <- 0:(N - 1) - A@i[A@p + 1][-(N + 1)]
  } else if (A@uplo == "L") {
    bw <- (A@i[A@p[-1]] - 0:(N - 1))[N:1]
  }
  Nenv <- sum(bw)
  ##xenv[i] tells us where the first entry of row i is stored
  xenv <- cumsum(c(1, bw)) - 1 ##for 0-based indices to be used in C
  env <- numeric(Nenv)
  diagvec <- numeric(N)
  if (Ccode) {
    .C("dsC2env", N, A@x, as.integer(A@i), as.integer(A@p),
       diagvec = diagvec, env = env, xenv = as.integer(xenv), DUP = FALSE)
  } else {
    for(i in 1:N) {
      if (bw[i] != 0) {
        env[(xenv[i] + 1):xenv[i + 1]] <- A[i, (i - bw[i]) : (i - 1)]
      }
    }
    diagvec <- diag(as(A, "dgCMatrix"))
  }
  
  return(list(diagvec = diagvec, env = env, xenv = as.integer(xenv),
              N = as.integer(N), Nenv = as.integer(Nenv),
              Nplus1 = as.integer(N + 1)))
}
