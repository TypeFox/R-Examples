#' apply() for big.matrix objects
#' @description \link{apply} for \code{\link[bigmemory]{big.matrix}} objects.
#' Note that the performance may be degraded (compared to \code{apply} with 
#' regular \R matrices) because of S4 overhead associated with extracting data 
#' from \code{big.matrix} objects.  This sort of limitation is unavoidable and 
#' would be the case (or even worse) with other "custom" data structures.  Of 
#' course, this would only be partically significant if you are applying over 
#' lengthy rows or columns.
#' @rdname apply-methods
#' @param X a big.matrix object.
#' @param MARGIN the margin. May be may only be 1 or 2, but otherwise 
#' conforming to what you would expect from \code{apply()}.
#' @param FUN the function to apply.
#' @param \dots other parameters to pass to the FUN parameter.
#' @docType methods
#' @export
#' @examples
#' library(bigmemory)
#' options(bigmemory.typecast.warning=FALSE)
#' x <- big.matrix(5, 2, type="integer", init=0,
#'                 dimnames=list(NULL, c("alpha", "beta")))
#' x[,] <- round(rnorm(10))
#' biganalytics::apply(x, 1, mean)
setMethod('apply', signature(X="big.matrix"),
  function(X, MARGIN, FUN, ...) return(bmapply(X, MARGIN, FUN, ...)))

bmapply <- function(X, MARGIN, FUN, ...)
{
  if (length(MARGIN)>1) 
    stop("MARGIN > 1 not supported with big.matrix objects.\n")
  FUN <- match.fun(FUN)
  dn.ans <- dimnames(X)[MARGIN]
  if (MARGIN==1) {
    d2 <- nrow(X)
    ans <- vector("list", nrow(X))
    for (i in 1:d2) {
      tmp <- FUN(X[i,], ...)
      if (!is.null(tmp)) ans[[i]] <- tmp
    }
  } else {
    if (MARGIN==2) {
      d2 <- ncol(X)
      ans <- vector("list", ncol(X))
      for (i in 1:d2) {
        tmp <- FUN(X[,i], ...)
        if (!is.null(tmp)) ans[[i]] <- tmp
      }
    } else {
      stop("Only MARGIN equal to 1 or 2 is supported for a big.matrix.\n")
    }
  }
  ans.list <- is.recursive(ans[[1]])
  l.ans <- length(ans[[1]])
  ans.names <- names(ans[[1]])
  if (!ans.list) 
    ans.list <- any(unlist(lapply(ans, length)) != l.ans)
  if (!ans.list && length(ans.names)) {
    all.same <- sapply(ans, function(x) identical(names(x), ans.names))
    if (!all(all.same)) ans.names <- NULL
  }
  if (ans.list) len.a <- d2
  else len.a <- length(ans <- unlist(ans, recursive = FALSE))
  if (len.a == d2) {
    if (length(dn.ans[[1]])) names(ans) <- dn.ans[[1]]
    return(ans)
  }
  if (len.a > 0 && len.a%%d2 == 0) {
    if (is.null(dn.ans)) 
      dn.ans <- vector(mode = "list", length(d2))
    dn.ans <- c(list(ans.names), dn.ans)
    return(array(ans, c(len.a%/%d2, d2), if (!all(sapply(dn.ans, 
                 is.null))) dn.ans))
  }
  return(ans)
}


