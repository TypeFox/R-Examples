"lARFIMA" <- function(z, phi = numeric(0), theta = numeric(0), dfrac = numeric(0), phiseas = numeric(0), 
    thetaseas = numeric(0), dfs = numeric(0), H = numeric(0), Hs = numeric(0), alpha = numeric(0), 
    alphas = numeric(0), period = 0, useC = 3) {
    
    if (is.null(phi) || any(is.na(phi))) 
        phi <- numeric(0)
    if (is.null(theta) || any(is.na(theta))) 
        theta <- numeric(0)
    if (is.null(phiseas) || any(is.na(phiseas))) 
        phiseas <- numeric(0)
    if (is.null(thetaseas) || any(is.na(thetaseas))) 
        thetaseas <- numeric(0)
    
    if (round(useC) != useC) 
        stop("non-integer useC")
    n <- length(z)
    if (useC == 0) 
        useC <- useCt <- F else if (useC == 1) {
        useC <- T
        useCt <- F
    } else if (useC == 2) {
        useC <- F
        useCt <- T
    } else if (useC == 3) 
        useC <- useCt <- T else stop("invalid useC")
    r <- as.double(tacvfARFIMA(phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, 
        thetaseas = thetaseas, dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, 
        period = period, maxlag = (n - 1), useCt = useCt))
    
    if (useC) {
        logl <- tryCatch(DLLoglikelihood(r, z), error = function(err) NULL)
        if (is.null(logl)) 
            logl <- -1e+10
    } else {
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
        logl <- -n/2 * log(S/n) - log(g)/2
    }
    
    return(logl)
} 
