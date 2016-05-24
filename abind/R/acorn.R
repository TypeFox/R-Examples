acorn <- function(x, n=6, m=5, r=1, ...) UseMethod('acorn')

acorn.default <- function(x, n=6, m=5, r=1, ..., addrownums = TRUE) {
    stopifnot(length(n) == 1L)
    stopifnot(length(m) == 1L)
    stopifnot(length(r) == 1L)
    if (is.null(dim(x)))
        if (n >= 0) return(head(x, n)) else return(tail(x, -n))
    n <- sign(n)*min(abs(n), nrow(x))
    ns <- if (n >= 0) seq(len=n) else seq(len=-n, to=nrow(x))
    as <- list(ns)
    if (length(dim(x))>1) {
        m <- sign(m)*min(abs(m), ncol(x))
        ms <- if (m >= 0) seq(len=m) else seq(len=-m, to=ncol(x))
        as <- c(as, list(ms))
    }
    if (length(dim(x))>2) {
        r <- sign(r)*min(abs(r), dim(x)[3])
        rs <- if (r >= 0) seq(len=r) else seq(len=-r, to=dim(x)[3])
        as <- c(as, list(rs))
    }
    if (length(dim(x))>3) {
        args <- list(...)
        for (i in seq(4, length(dim(x)))) {
            if (length(args) >= i-3)
                a <- args[[i-3]]
            else
                a <- 1
            a <- sign(a)*min(abs(a), dim(x)[i])
            as <- c(as, list(if (a >= 0) seq(len=a) else seq(len=-a, to=dim(x)[i])))
        }
    }
    y <- asub(x, as, drop=FALSE)
    if (addrownums) {
        if (is.null(dimnames(y)))
            dimnames(y) <- vector('list', length(dim(x)))
        for (i in seq(len=length(dim(x))))
            if (is.null(dimnames(y)[[i]]))
                dimnames(y)[[i]] <- paste0('[', as[[i]], ']')
    }
    y
}
