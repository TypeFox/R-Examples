dropvar <-
function(x, tol=1e-7, LAPACK=FALSE, silent=FALSE)
### works if ncol(X) >= 0 and nrow(X) >= 0
{
  ## test and match arguments:
  stopifnot(is.matrix(x))
  silent <- as.logical(silent)[1]
  ## perform the qr-decomposition of X using LINPACK methods:
  qr.X <- qr(x, tol=tol, LAPACK=LAPACK)
  if(qr.X$rank == NCOL(x))
    return(x) ## return x if x has full column rank
  if(!silent){ ## message the no. of dropped columns:
    cat("regressor-matrix is column rank deficient, so dropping", NCOL(x) - qr.X$rank, "regressors\n")
    cat("\n")
  }
#OLD:
#    message(gettextf("regressor-matrix is column rank deficient, so dropping %d regressors",
#                     NCOL(x) - qr.X$rank))
  ## return the columns correponding to the first qr.x$rank pivot
  ## elements of x:
  newX <- x[, qr.X$pivot[1:qr.X$rank], drop = FALSE]
  ## did we succeed? stop-if-not:
  if(qr.X$rank != qr(newX)$rank)
    stop(gettextf("determination of full column rank design matrix failed"),
         call. = FALSE)
  return(newX)
}
