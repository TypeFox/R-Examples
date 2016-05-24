###
### $Id: magic.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Create a magic square.
###


##-----------------------------------------------------------------------------
magic <- function(n) {
    if (!is.numeric(n)) {
        stop(sprintf("argument %s must be numeric", sQuote("n")))
    } else if (!(length(n) == 1)) {
        stop(sprintf("argument %s must be of length 1", sQuote("n")))
    }

    oddOrder <- function(n) {
        ans <- matlab::meshgrid(1:n)
        J <- ans$x
        I <- ans$y
        A <- matlab::mod(I + J - (n + 3) / 2, n)
        B <- matlab::mod(I + 2 * J - 2, n)
        M <- n * A + B + 1;

        return(M)
    }

    doublyEvenOrder <- function(n) {
        ans <- matlab::meshgrid(1:n)
        J <- ans$x
        I <- ans$y
        K <- matlab::fix(matlab::mod(I, 4) / 2) ==
             matlab::fix(matlab::mod(J, 4) / 2)
        M <- t(matlab::reshape(as.matrix(1:(n * n)), n, n))
        M[K] = n * n + 1 - M[K]

        return(M)
    }

    singlyEvenOrder <- function(n) {
        p <- n / 2
        fun <- sys.function(sys.parent())
        M <- fun(p)             # same as matlab::magic(p)
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
             if (matlab::mod(n, 2) == 1) {
                 oddOrder(n)
             } else if (matlab::mod(n, 4) == 0) {
                 doublyEvenOrder(n)
             } else {
                 singlyEvenOrder(n)
             }
         }

    return(M)
}

