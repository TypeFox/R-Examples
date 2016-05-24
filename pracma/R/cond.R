##
##  c o n d . R  Matrix Condition
##


cond <- function(M, p = 2) {
    if (length(M) == 0)
        return(0)
    if (!is.numeric(M))
        stop("Argument 'M' must be a numeric matrix.")
    if (is.vector(M))
        M <- matrix(c(M), nrow = length(M), ncol = 1)

    if (length(M) == 0) return(c())
    if (ncol(M) != nrow(M) && p != 2)
        stop("Matrix 'M' must be square if p = 2.")

    if (p == 2) {
        s <- svd(M)$d
        cnd <- if (any(s == 0)) Inf else max(s) / min(s)
    } else {
        stop("At the moment, p-norms other than p = 2 are not implemented.")
        #cnd <- norm(M, p) * norm(inv(M), p)
    }
    return(cnd)
}
