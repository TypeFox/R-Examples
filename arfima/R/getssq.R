getsigma2muhat <- function(z, phi = numeric(0), phiseas = numeric(0), theta = numeric(0), 
    thetaseas = numeric(0), H = numeric(0), period = 0, dfrac = numeric(0), alpha = numeric(0), 
    dfs = numeric(0), Hs = numeric(0), alphas = numeric(0)) {
    
    n <- length(z)
    r <- tacvfARFIMA(phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, thetaseas = thetaseas, 
        dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, period = period, maxlag = n - 
            1, useCt = T)
    error <- numeric(n)
    sigmasq <- numeric(n)
    error[1] <- z[1]
    sigmasq[1] <- r[1]
    phee <- r[2]/r[1]
    error[2] <- z[2] - phee * z[1]
    sigmasqkm1 <- r[1] * (1 - phee^2)
    sigmasq[2] <- sigmasqkm1
    g <- r[1] * sigmasqkm1
    
    for (k in 2:(n - 1)) {
        phikk <- (r[k + 1] - phee %*% rev(r[2:k]))/sigmasqkm1
        sigmasqk <- sigmasqkm1 * (1 - phikk^2)
        phinew <- phee - phikk * rev(phee)
        phee <- c(phinew, phikk)
        sigmasqkm1 <- sigmasqk
        g <- g * sigmasqk
        error[k + 1] <- z[k + 1] - crossprod(phee, rev(z[1:k]))
        sigmasq[k + 1] <- sigmasqk
    }
    S <- sum((error * error)/sigmasq)
    
    sigma2 <- S/(n - length(phi) - length(theta) - length(thetaseas) - length(phiseas) - 
        length(H) - length(dfrac) - length(dfs) - length(alpha) - length(alphas))
    
    sigma2
} 
