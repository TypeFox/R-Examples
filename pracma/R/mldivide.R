##
##  m l d i v i d e . R  Matlab backslash Operator
##


mldivide <- function(A, B, pinv = TRUE) {
    stopifnot(is.numeric(A) || is.complex(A),
              is.numeric(B) || is.complex(B))
    if (is.vector(A)) A <- as.matrix(A)
    if (is.vector(B)) B <- as.matrix(B)
    if (nrow(A) != nrow(B))
        stop("Matrices 'A' and 'B' must have the same number of rows.")
    if (pinv) {
        pinv(t(A) %*% A) %*% t(A) %*% B
    } else {
        qr.solve(A, B)
    }
}

mrdivide <- function(A, B, pinv = TRUE) {
    stopifnot(is.numeric(A) || is.complex(A),
              is.numeric(B) || is.complex(B))
    if (is.vector(A)) A <- t(A)
    if (is.vector(B)) B <- t(B)
    if (ncol(A) != ncol(B))
        stop("Matrices 'A' and 'B' must have the same number of columns.")

    t(mldivide(t(B), t(A), pinv = pinv))
}
