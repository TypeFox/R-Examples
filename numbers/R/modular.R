##
##  m o d u l a r . R  Modular functions
##

mod <- function(n, m) {
    stopifnot(is.numeric(n), is.numeric(m))
    if (length(n) == 1) {
        n <- rep(n, length(m))
    } else if (length(m) == 1) {
        m <- rep(m, length(n))
    }
    ln <- length(n); lm <- length(m)
    if (ln != lm)
        stop("Arguments 'n', 'm' must be scalars or have the same length.")
    if (any(floor(n) != ceiling(n)) || any(floor(m) != ceiling(m)))
        stop("Arguments 'n', 'm' must be integers or vectors of integers.")

    k <- ifelse(m != 0, n %% m, n)
    return(k)
}


rem <- function(n, m) {
    stopifnot(is.numeric(n), is.numeric(m))
    if (length(n) == 1) {
        n <- rep(n, length(m))
    } else if (length(m) == 1) {
        m <- rep(m, length(n))
    }
    ln <- length(n); lm <- length(m)
    if (ln != lm)
        stop("Arguments 'n', 'm' must be scalars or have the same length.")
    if (any(floor(n) != ceiling(n)) || any(floor(m) != ceiling(m)))
        stop("Arguments 'n', 'm' must be integers or vectors of integers.")

    k <- ifelse(m != 0, n %% m, NaN)
    k <- ifelse(m != 0 & sign(n) != sign(m) & k != 0, k - m, k)
    return(k)
}


div <- function(n, m) {
    stopifnot(is.numeric(n), is.numeric(m))
    if (length(n) == 1) {
        n <- rep(n, length(m))
    } else if (length(m) == 1) {
        m <- rep(m, length(n))
    }
    ln <- length(n); lm <- length(m)
    if (ln != lm)
        stop("Arguments 'n', 'm' must be scalars or have the same length.")
    if (any(floor(n) != ceiling(n)) || any(floor(m) != ceiling(m)))
        stop("Arguments 'n', 'm' must be integers or vectors of integers.")

    k <- n %/% m
    return(k)
}


modinv <- function(n, m) {
    stopifnot(is.numeric(n), is.numeric(m))
    v <- extGCD(n, m)
    if (v[1] == 0 || v[1] > 1) return(NA)
    if (v[2] >= 0) v[2] else v[2] + m
}

modlin <- function(a, b, n) {
    stopifnot(is.numeric(a), is.numeric(b), is.numeric(n))
    if (length(a) != 1 || length(b) != 1 || length(n) != 1 ||
        floor(a) != ceiling(a) || floor(b) != ceiling(b) || floor(n) != ceiling(n) ||
        a < 1 || b < 1 || n < 1)
        stop("All inputs 'a', 'b', 'n' must be integers.")

    def <- extGCD(a, n)
    d <- def[1]; e <- def[2]; f <- def[3]

    x <- c()
    if (b %% d == 0) {
        x0 <- (e * (b/d)) %% n
        for (i in 0:(d-1)) {
            x <- c(x, (x0 + i*(n/d)) %% n)
        }
    }
    return(x)
}
