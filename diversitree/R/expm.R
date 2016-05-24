## This is significantly less efficient than proper checkpointing
## within ddexpmv -> ddexpmvi.  However, the interface should be OK.
expmv.expokit.dense <- function(Q, t, v) {
  n <- as.integer(length(v))
  nt <- length(t)
  if ( nt != 1 ) {
    ret <- matrix(0, length(v), nt)
    for ( i in seq_len(nt) ) {
      res <- .Fortran("ddexpmv", Q, n, v, t[i], out=numeric(n),
                      iflag=numeric(1), package="diversitree")
      if ( res$iflag != 0 )
        stop("Expokit failed with code ", res$iflag)
      ret[,i] <- res$out
    }
  } else {
    res <- .Fortran("ddexpmv", Q, n, v, t, out=numeric(n),
                    iflag=numeric(1), PACKAGE="diversitree")
    if ( res$iflag != 0 )
      stop("Expokit failed with code ", res$iflag)
    ret <- res$out
  }
  ret
}

expmv.expokit.sparse <- function(Q.sparse, t, v, tol=1e-8) {
  if ( is.matrix(Q.sparse) )
    Q.sparse <- expm.expokit.sparse.pars(Q.sparse)
  Q     <- Q.sparse$Q
  iq    <- Q.sparse$iq
  jq    <- Q.sparse$jq
  qnorm <- Q.sparse$qnorm

  n <- as.integer(length(v))
  lt <- as.integer(length(t))

  res <- .Fortran("dsexpmvi", Q, n, iq, jq, length(iq),
                  qnorm, v, t, lt, tol, out=numeric(n*lt),
                  iflag=integer(1), PACKAGE="diversitree")
  if ( res$iflag != 0 )
    stop("Expokit failed with code ", res$iflag)

  matrix(res$out, n, lt)
}

## Convert a sparse Q matrix into COO format in preparation for
## running expmv(sparse) on it.
expm.expokit.sparse.pars <- function(Q) {
  idx <- which(Q != 0, TRUE)
  nz <- as.integer(nrow(idx))
  qq <- Q[idx]
  list(Q=qq, iq=idx[,1], jq=idx[,2], nz=nrow(idx), qnorm=max(abs(qq)))
}

expm.dense <- function(Q, t) {
  if ( !is.matrix(Q) )
    stop("Q must be a matrix")
  n <- as.integer(nrow(Q))
  if ( n != ncol(Q) )
    stop("Q must be a square matrix")
  if ( n < 1 )
    stop("Q must have positive dimensions")
  t <- as.numeric(t)
  if ( length(t) != 1 )
    stop("t must be a scalar")
  if ( !is.finite(t) )
    stop("t must be finite")
  
  matrix(.Fortran("dexpmf", Q, n, t, out=numeric(n*n),
                  iflag=numeric(1), PACKAGE="diversitree")$out, n, n)
}
