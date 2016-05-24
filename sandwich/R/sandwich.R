sandwich <- function(x, bread. = bread, meat. = meat, ...)
{
  if(is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"
  if(is.function(bread.)) bread. <- bread.(x)
  if(is.function(meat.)) meat. <- meat.(x, ...)
  n <- NROW(estfun(x))
  return(1/n * (bread. %*% meat. %*% bread.))
}

meat <- function(x, adjust = FALSE, ...)
{
  if(is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"
  psi <- estfun(x, ...)
  k <- NCOL(psi)
  n <- NROW(psi)
  rval <- crossprod(as.matrix(psi))/n
  if(adjust) rval <- n/(n-k) * rval
  rownames(rval) <- colnames(rval) <- colnames(psi)
  return(rval)
}

vcovOPG <- function(x, adjust = FALSE, ...)
{
  if(is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"
  psi <- estfun(x, ...)
  k <- NCOL(psi)
  n <- NROW(psi)
  rval <- chol2inv(qr.R(qr(psi)))
  if(adjust) rval <- n/(n-k) * rval
  rownames(rval) <- colnames(rval) <- colnames(psi)
  return(rval)
}
