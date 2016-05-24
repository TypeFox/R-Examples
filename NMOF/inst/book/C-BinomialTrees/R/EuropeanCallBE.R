# EuropeanCallBE.R -- version 2010-12-28
EuropeanCallBE <- function(S0, X, r, T, sigma, M) {
    # compute constants
    dt <- T/M
    u <- exp(sigma * sqrt(dt))   
    d <- 1/u
    p <- (exp(r * dt) - d)/(u - d)
    
    # initialize asset prices at maturity (period M)
    C <- pmax(S0 * d^(M:0) * u^(0:M) - X, 0)
    
    # log/cumsum version
    csl <- cumsum(log(c(1,1:M)))
    tmp <- csl[M+1] - csl - csl[(M+1):1] + 
        log(p) * (0:M) + log(1-p) * (M:0)
    C0 <- exp(-r * T) * sum(exp(tmp) * C)
    C0
}