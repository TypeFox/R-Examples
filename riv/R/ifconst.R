# 6) Some constants needed for the IF of S-est (Lopuhaa, 1989, p.1672)
#
# r=dim of data for IF of S-est.
IFconst <- function(r, c) {
    x.sample <- mvrnorm(10000, rep(0, r), diag(r))^2
    
    # norms of multivariate normal random vectors, dim(x)=1000
    x.norms <- sqrt(rowSums(x.sample))

    # nu
    fnu <- (1 - 1/r) * psi.bisquare(x.norms, c) +
           (1/r) * psi.prime(x.norms, c)
    nu <- mean(fnu)

    # b0
    b <- mean(rho.biweight(x.norms, c))

    # gamma3
    fg3 <- psi.bisquare(x.norms, c) * x.norms^2
    g3 <- mean(fg3)

    # gamma1
    fg1 <- psi.prime(x.norms, c) * x.norms^2 +
           (r + 1) * psi.bisquare(x.norms, c) * x.norms^2
    g1 <- mean(fg1)/(r + 2)

    # gamma 2
    fg2 <- 2 * psi.prime(x.norms, c) * x.norms^2 +
           r * psi.bisquare(x.norms, c) * x.norms^2
    g2 <- mean(fg2)/(2 * r * (r + 2))
    
    list(nu = nu, b0 = b, gamma1 = g1, gamma2 = g2, gamma3 = g3)
}
