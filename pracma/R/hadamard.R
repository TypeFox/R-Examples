##
##  h a d a m a r d . R  Hadamard Matrix
##


hadamard <- function(n) {
    if (!is.numeric(n) || length(n) != 1 ||
        floor(n) != ceiling(n) || n <= 1)
        stop("Argument 'n' must be a positiv integer.")

    x <- log2(c(n, n/12, n/20))
    k <- which(floor(x) == ceiling(x))

    if (length(k) == 0)
        stop("Argument 'n' is not of the form 2^e, 12*2^e , or 20*2^e.")

    e <- x[k]
    if (k == 1) {
        H <- c(1)
    } else if (k == 2) {
        H <- ones(12, 12)
        H[2:12, 2:12] <- Toeplitz(c(-1, -1, 1, -1, -1, -1, 1, 1, 1, -1, 1),
                                  c(-1, 1, -1, 1, 1, 1, -1, -1, -1, 1, -1))
    } else if (k == 3) {
        H <- ones(20, 20)
        H[2:20, 2:20] <-
  hankel(c(-1, -1, 1, 1, -1, -1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1, -1, -1, 1),
         c(1, -1, -1, 1, 1, -1, -1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1, -1, -1))
    }
    if (e >= 1) {
        for ( i in 1:e) {
            H <- matrix(c(1, 1, 1, -1), 2, 2) %x% H  # Kronecker product
        }
    }
    return(H)
}


Toeplitz <- function(a, b = a) {
    if (!is.vector(a) || !is.vector(b))
        stop("Arguments 'a' and 'b' must be vectors.")
    if (a[1] != b[1])
        warning("First elements of vectors 'a', 'b' are not equal.")

    n <- length(a)
    m <- length(b)
    T <- matrix(nrow = n, ncol = m)
    T[1, ] <- b
    T[, 1] <- a
    for (i in 2:n) {
        T[i, 2:m] <- T[i-1, 1:(m-1)]
    }
    return(T)
}
