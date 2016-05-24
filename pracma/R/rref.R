##
##  r r e f . R  Reduced Row Echelon Form
##


rref <- function(A) {
    stopifnot(is.numeric(A))
    if (!is.matrix(A))
        stop("Input parameter 'A' must be a matrix.")

    nr <- nrow(A); nc <- ncol(A)
    tol <- eps() * max(nr, nc) * max(abs(A))

    r <- 1
    for (i in 1:nc) {
        pivot <- which.max(abs(A[r:nr, i]))
        pivot <- r + pivot - 1
        m <- abs(A[pivot, i])
        if (m <= tol) {
            A[r:nr, i] <- 0  # zeros(nr-r+1, 1)
        } else {
            A[c(pivot, r), i:nc] <- A[c(r, pivot), i:nc]
            A[r, i:nc] <- A[r, i:nc] / A[r, i]
            if (r == 1) {
                ridx <- c((r+1):nr)
            } else if (r == nr) {
                ridx <- c(1:(r-1))
            } else {
                ridx <- c(1:(r-1), (r+1):nr)
            }
            A[ridx, i:nc] <- A[ridx, i:nc] - 
                             A[ridx, i, drop=FALSE] %*% A[r, i:nc, drop=FALSE]
            if (r == nr) break
            r <- r+1
        }
    }
    A[abs(A) < tol] <- 0
    return(A)
}

