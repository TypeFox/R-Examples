outerList <- function(X, Y, FUN, ...) {
    nX <- length(X)
    nY <- length(Y)
    
    x <- rep(X, nY)
    y <- rep(Y, each = nX)
    
    res <- mapply(FUN, x, y, ...)
    
    matrix(res, nrow = nX, dimnames = list(names(X), names(Y)))
}

zipLists <- function(...)
    mapply(function(...) list(...), ..., SIMPLIFY = FALSE)
