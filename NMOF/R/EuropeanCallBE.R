EuropeanCallBE <- function(S0, X, r, tau, sigma, M = 101) {

    ## compute constants
    dt <- tau / M
    u <- exp(sigma*sqrt(dt))   
    d <- 1 / u
    p <- (exp(r * dt) - d) / (u - d)
    
    ## initialise asset prices at maturity (period M)
    C <- pmax(S0 * d^(M:0) * u^(0:M) - X, 0)
    
    ## log/cumsum version; see Higham (2002) in References
    csl <- cumsum(log(c(1,1:M)))
    tmp <- csl[M+1] - csl - csl[(M+1):1] + log(p)*(0:M) + log(1-p)*(M:0)
    exp(-r*tau)*sum(exp(tmp)*C)
}

