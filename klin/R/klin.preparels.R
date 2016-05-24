`klin.preparels` <-
function(A) {
  ## dimensions of matrices
  m <- sapply(A,nrow)
  n <- sapply(A,ncol)
  K <- length(A)
  ## check consistency
  if (!is.list(A))
    stop("argument should be a list")
  if (any(n>m))
    stop("incompatible dimensions")
  ## pick square matrices
  square <- n==m
  ## calculate left right hand side
  lhs <- rhs <- list(K)
  for (i in seq_len(K)) {
    if (square[i]) {
      lhs[[i]] <- A[[i]]
      rhs[[i]] <- Diagonal(n[[i]])
    } else {
      lhs[[i]] <- crossprod(A[[i]])
      rhs[[i]] <- A[[i]]
    }
  }
  ## return
  structure(list(lhs=lhs,rhs=rhs),class="klin.prepls")
}

