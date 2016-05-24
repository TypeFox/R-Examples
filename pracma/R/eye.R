##
##  e y e . R  Generate basic Matrices
##


eye <- function(n, m = n) {
    stopifnot(is.numeric(n), length(n) == 1,
              is.numeric(m), length(m) == 1)
    n <- floor(n)
    m <- floor(m)
    if (n <= 0 || m <= 0) return(matrix(NA, 0, 0))
    else                  return(base::diag(1, n, m))
}

ones <- function(n, m = n) {
        stopifnot(is.numeric(n), length(n) == 1,
                  is.numeric(m), length(m) == 1)
        n <- floor(n)
        m <- floor(m)
        if (n <= 0 || m <= 0) return(matrix(1, 0, 0))
        else                  return(matrix(1, n, m))
}

zeros <- function(n, m = n) {
    stopifnot(is.numeric(n), length(n) == 1,
              is.numeric(m), length(m) == 1)
    n <- floor(n)
    m <- floor(m)
    if (n <= 0 || m <= 0) return(matrix(0, 0, 0))
    else                  return(matrix(0, n, m))
}
