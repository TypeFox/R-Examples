###
### m a g i c . R -- Create a Magic Square
###


magic <- function(n) {
    if (!is.numeric(n)) {
        stop("Argument 'n' must be numeric.")
    } else if (!(length(n) == 1)) {
        stop("Argument 'n' must be of length 1.")
    }

    oddOrder <- function(n) {
        J <- matrix(rep(1:n, each = n),  n, n)
        I <- matrix(rep(1:n, times = n), n, n)
        A <- (I + J - (n + 3) / 2) %% n
        B <- (I + 2 * J - 2) %% n
        M <- n * A + B + 1;

        return(M)
    }

    doublyEvenOrder <- function(n) {
        J <- matrix(rep(1:n, each = n),  n, n)
        I <- matrix(rep(1:n, times = n), n, n)
        K <- trunc((I %% 4) / 2) == trunc((J %% 4) / 2)
        M <- t(matrix(1:(n*n), n, n))
        # M <- t(pracma::Reshape(as.matrix(1:(n * n)), n, n))
        M[K] = n * n + 1 - M[K]

        return(M)
    }

    singlyEvenOrder <- function(n) {
        p <- n / 2
        M <- magic(p)
        M <- rbind(cbind(M, M + 2 * p ^ 2),
                   cbind(M + 3 * p ^ 2, M + p ^ 2))
        if (!(n == 2)) {
            i <- t(1:p)
            k <- (n - 2) / 4
            j <- c(1:k, if ((n - k + 2) <= n) (n - k + 2):n)
            M[cbind(i, i + p), j] <- M[cbind(i + p, i), j]
            i <- k + 1
            j <- c(1, i)
            M[cbind(i, i + p), j] <- M[cbind(i + p, i), j]
        }

        return(M)
    }

    n <- floor(n)
    M <- if (n <= 0) {
             matrix(numeric(0), 0, 0)       # degenerate
         } else if (n == 1) {
             matrix(as.numeric(1), 1, 1)    # degenerate
         } else {
             if (pracma::mod(n, 2) == 1) {
                 oddOrder(n)
             } else if (pracma::mod(n, 4) == 0) {
                 doublyEvenOrder(n)
             } else {
                 singlyEvenOrder(n)
             }
         }

    if (n == 2)                             # impossible
        warning("There is no magic square of order 2.")
    return(M)
}

