svdFunc <- function(x, y, lambda) {

  k <- length(lambda)

  Ds <- try(svd(x = x), silent = TRUE)
  if( is(Ds, "try-error") ) {
    stop("Unable to obtain singular value decomposition.", 
         call. = FALSE)
  }

  ld <- length(Ds$d)

  div <- Ds$d*Ds$d + rep(lambda, rep(ld, k))

  rhs <- t(Ds$u) %*% y
  tmp <- drop(Ds$d * rhs) / div
  dim(tmp) <- c(ld, k)

  coef <- Ds$v %*% tmp

  return( coef )

}
