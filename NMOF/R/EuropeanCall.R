EuropeanCall <- function(S0, X, r, tau, sigma, M = 101) {

    ## compute constants
    f7 <- 1
    dt <- tau / M
    v <- exp(-r * dt)
    u <- exp(sigma * sqrt(dt))
    d <- 1 /u
    p <- (exp(r * dt) - d) / (u - d)

    ## initialise asset prices at maturity (period M)
    S <- numeric(M + 1)
    S[f7 + 0] <- S0 * d^M
    for (j in 1:M) S[f7 + j] <- S[f7 + j - 1] * u / d

    ## initialise option values at maturity (period M)
    C <- pmax(S - X, 0)

    ## step back through the tree
    for (i in seq(M - 1, 0, by = -1))
        C <- v * (p * C[(1+f7):(i+1+f7)] + (1-p) * C[(0+f7):(i+f7)])
    C
}

