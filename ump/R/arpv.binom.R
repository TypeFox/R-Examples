arpv.binom <- function(x, n, p, plot = TRUE, ...) {

    if (! is.numeric(x)) stop("x not numeric")
    if (! is.numeric(n)) stop("n not numeric")
    if (! is.numeric(p)) stop("p not numeric")

    if (length(x) != 1) stop("x not scalar")
    if (length(n) != 1) stop("n not scalar")
    if (length(p) != 1) stop("p not scalar")

    if (as.integer(x) != x) stop("x not integer")
    if (as.integer(n) != n) stop("n not integer")

    if (! (n > 0)) stop("n not positive")
    if (! (0 <= x & x <= n)) stop("x not in 0, ..., n")
    if (! (0 < p & p < 1)) stop("p not in (0, 1)")

    mu <- n * p

    foo <- sign(x - mu)

    if (foo == 0) {

        p <- dbinom(x, n, p)
        alpha <- c(1 - p, 1)
        phi <- c(0, 1)

    } else {

    if (foo > 0) {
        c1 <- seq(0, floor(mu))
        c2 <- x
    } else {
        c1 <- x
        c2 <- seq(n, ceiling(mu))
    }

    p1 <- dbinom(c1, n, p)
    p2 <- dbinom(c2, n, p)
    P1 <- pbinom(c1 - 1, n, p)
    P2 <- pbinom(c2, n, p, lower.tail = FALSE)
    M1 <- mu * pbinom(c1 - 2, n - 1, p)
    M2 <- mu * pbinom(c2 - 1, n - 1, p, lower.tail = FALSE)

    alpha.max.1 <- (p1 * (c2 - c1) - M1 - M2 + c2 * P1 + c2 * P2) / (c2 - mu)
    alpha.max.2 <- (p2 * (c1 - c2) - M1 - M2 + c1 * P2 + c1 * P1) / (c1 - mu)
    alpha.min.1 <- (- M1 - M2 + c2 * P1 + c2 * P2) / (c2 - mu)
    alpha.min.2 <- (- M1 - M2 + c1 * P1 + c1 * P2) / (c1 - mu)
    alpha.max <- pmin(alpha.max.1, alpha.max.2)
    alpha.min <- pmax(alpha.min.1, alpha.min.2)
    alpha.max[c2 - c1 == 1] <- 1
    alpha.min[c2 - c1 == n] <- 0

    gamma.max.1 <- (alpha.max * (c2 - mu) + (M1 - c2 * P1) + (M2 - c2 * P2)) /
        (p1 * (c2 - c1))
    gamma.max.2 <- (alpha.max * (c1 - mu) + (M2 - c1 * P2) + (M1 - c1 * P1)) /
        (p2 * (c1 - c2))
    gamma.min.1 <- (alpha.min * (c2 - mu) + (M1 - c2 * P1) + (M2 - c2 * P2)) /
        (p1 * (c2 - c1))
    gamma.min.2 <- (alpha.min * (c1 - mu) + (M2 - c1 * P2) + (M1 - c1 * P1)) /
        (p2 * (c1 - c2))

    inies <- alpha.max > alpha.min

    alpha.min <- alpha.min[inies]
    alpha.max <- alpha.max[inies]

    alpha.min[alpha.min < 0] <- 0
    alpha.max[alpha.max > 1] <- 1

    if (foo > 0) {
        ##### c2 == x #####
        phi.min <- gamma.min.2[inies]
        phi.max <- gamma.max.2[inies]
    } else {
        ##### c1 == x #####
        phi.min <- gamma.min.1[inies]
        phi.max <- gamma.max.1[inies]
    }

    phi.min[phi.min < 0] <- 0
    phi.max[phi.max > 1] <- 1

    alpha.foo <- cbind(c(alpha.min, NA), c(NA, alpha.max))
    phi.foo <- cbind(c(phi.min, NA), c(NA, phi.max))

    alpha <- apply(alpha.foo, 1, mean, na.rm = TRUE)
    phi <- apply(phi.foo, 1, mean, na.rm = TRUE)

    }

    if (plot)
        arpv.plot(alpha, phi, ...)

    return(invisible(list(alpha = alpha, phi = phi)))

}
