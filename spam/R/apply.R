# This is file ../spam/R/apply.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     


















# primitive apply function.

apply.spam <- function(X, MARGIN=NULL, FUN, ...){
  if (!is.spam(X)) stop("'X' must be of type 'spam'")
  if (!is.null(dimnames(X))) warning("dimnames are stripped")
  FUN <- match.fun(FUN)
  d <- dim(X)
  d.ans <- d[MARGIN]
  dn.ans <- NULL
  if (is.null(MARGIN)|| length(MARGIN)==2){
    ans <- FUN(X@entries,...)
    if (length( ans)!=length( X@entries))
      stop("'FUN' does not return an appropriate vector")
    if (any(!is.finite(ans))) {
      warning("'NA/NaN/Inf' coerced to zero")
      ans[!is.finite(ans)] <- 0
    }

    X@entries <- ans
    return(X)
  }
  ans <- vector("list",d.ans)
  if (MARGIN==1) {
    for (i in 1:d[1])
      ans[[i]] <- FUN(X[i,,drop=F]@entries,...)
  } else   if (MARGIN==2) {  
    for (i in 1:d[2])
      ans[[i]] <- FUN(X[,i,drop=F]@entries,...)
  } else stop("'MARGIN' must be 1, 2 or c(1,2)")


  # Block very similar to 'apply'
  d2 <- prod(d.ans)
  ans.list <- is.recursive(ans[[1]])
  l.ans <- length(ans[[1]])
  ans.names <- names(ans[[1]])

  if (!ans.list){
    ans.list <- any(unlist(lapply(ans, length)) != l.ans)
  }
  if (!ans.list && length(ans.names)) {
    all.same <- sapply(ans, function(x) identical(names(x), ans.names))
    if (!all(all.same))
      ans.names <- NULL
  }
  len.a <- if (ans.list)     d2   else length(ans <- unlist(ans, recursive = FALSE))
  if (length(MARGIN) == 1 && len.a == d2) 
    return(ans)
  if (len.a == d2)
    return(array(ans, d.ans))
  if (len.a > 0 && len.a%%d2 == 0) {
    dn.ans <- vector(mode = "list", length(d.ans))
    dn.ans <- c(list(ans.names), dn.ans)
    return(array(ans, c(len.a%/%d2, d.ans), if (!all(sapply(dn.ans,
                                                            is.null))) dn.ans))
  }
  return(ans)

    
}








