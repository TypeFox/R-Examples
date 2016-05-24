##
##  m o d u l a r . R  Modular functions
##

ceil <- function(n) ceiling(n)


Fix <- function(n) trunc(n)


mod <- function(n, m) {
    stopifnot(is.numeric(n), is.numeric(m))
    if (length(m) != 1 || floor(m) != ceiling(m))
        stop("Argument 'm' must be an integer scalar.")
    if (m == 0) return(n)
    else        return(n %% m)
}


rem <- function(n, m) {
    stopifnot(is.numeric(n), is.numeric(m))
    if (length(m) != 1 || floor(m) != ceiling(m))
        stop("Argument 'm' must be an integer scalar.")
    if (m == 0) return(NaN)
    k <- mod(n, m)
    if (sign(n) * sign(m) < 0) {
        k <- k - m
    }
    return(k)
}


idivide <- function(n, m, rounding = c("fix", "floor", "ceil", "round")) {
    stopifnot(is.numeric(n), is.numeric(m))
    rounding <- match.arg(rounding)
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

    k <- n / m
    if      (rounding == "fix")   k <- Fix(k)
    else if (rounding == "floor") k <- floor(k)
    else if (rounding == "ceil")  k <- ceil(k)
    else if (rounding == "round") k <- round(k)
    else
        stop("Rounding must be one of 'fix', 'floor', 'ceil', 'round'.")
    return(k)
}


gcd <- function(a, b, extended = FALSE) {
    stopifnot(is.numeric(a), is.numeric(b))
    if (length(a) == 1) {
        a <- rep(a, times=length(b))
    } else if (length(b) == 1) {
        b <- rep(b, times=length(a))
    }
    n <- length(a)

    e <- d <- g <- numeric(n)
    for (k in 1:n) {
        u <- c(1, 0, abs(a[k]))
        v <- c(0, 1, abs(b[k]))
        while (v[3] != 0) {
            q <- floor(u[3]/v[3])
            t <- u - v*q
            u <- v
            v <- t
        }
        e[k] <- u[1] * sign(a[k])
        d[k] <- u[2] * sign(a[k])
        g[k] <- u[3]
    }

    if (extended) {
        return(list(g = g, c = e, d = d))
    } else {
        return(g)
    }
}


Lcm <- function(a, b) {
    stopifnot(is.numeric(a), is.numeric(b))
    if (length(a) == 1) {
        a <- rep(a, times=length(b))
    } else if (length(b) == 1) {
        b <- rep(b, times=length(a))
    }

    g <- gcd(a, b, extended = FALSE)
    return( a / g * b )
}
