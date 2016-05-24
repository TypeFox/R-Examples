##
##  g a u s s L e g e n d r e . R  Gauss-Legendre et al. Quadrature Formula
##


gaussLegendre <- function(n, a, b) {
    stopifnot(is.numeric(a), length(a) == 1, is.numeric(b), length(b) == 1,
              is.numeric(n), length(n) == 1, n >= 2)

    i <- seq(1, n-1, by = 1)
    d <- i / sqrt(4*i^2 - 1)

    E <- eigen(Diag(d, 1) + Diag(d, -1), symmetric = TRUE)
    L <- E$values
    V <- E$vectors

    inds <- order(L)
    x <- L[inds]
    V <- t(V[, inds])
    w <- 2 * V[, 1]^2

    x <- 0.5 * ((b-a)*x + a+b)
    w <- -0.5 * (a-b)*w

    return(list(x = x, w = w))
}


gaussHermite <- function(n) {
    stopifnot(is.numeric(n), length(n) == 1, n >= 2)

    i <- seq(1, n-1, by = 1)
    d <- sqrt(i/2)

    E <- eigen(Diag(d, 1) + Diag(d, -1), symmetric = TRUE)
    L <- E$values
    V <- E$vectors

    inds <- order(L)
    x <- L[inds]
    V <- t(V[, inds])
    w <- sqrt(pi) * V[, 1]^2

    return(list(x = x, w = w))
}


gaussLaguerre <- function(n, a = 0) {
    stopifnot(is.numeric(n), length(n) == 1, n >= 2)
    stopifnot(is.numeric(a), length(a) == 1, a >= 0)

    i <- 1:n
    d <- (2*i - 1) + a
    b <- sqrt( i[1:(n-1)] * ((1:(n-1)) + a) )

    E <- eigen(Diag(d) + Diag(b, 1) + Diag(b, -1))
    L <- E$values
    V <- E$vectors

    inds <- order(L)
    x <- L[inds]
    V <- t(V[, inds])

    w <- gamma(a + 1) * V[, 1]^2

    return(list(x = x, w = w))
}
