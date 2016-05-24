Boot.arfima <- function(obj, R = 1, pred = FALSE, seed = NA, ...) {
    
    if (length(seed) > 1 && length(seed) != R) {
        warning("set seed length is not equal to R: using only first seed")
        seed <- seed[1]
    }
    
    n <- length(obj$z)
    
    m <- length(obj$modes)
    zinit <- NULL
    dint <- obj$dint
    dseas <- obj$dseas
    res <- vector("list", m)
    
    for (i in 1:m) {
        
        cmode <- obj$modes[[i]]
        phi <- cmode$phi
        theta <- cmode$theta
        phiseas <- cmode$phiseas
        thetaseas <- cmode$thetaseas
        period <- cmode$period
        dfrac <- cmode$dfrac
        dfs <- cmode$dfs
        dint <- cmode$dint
        dseas <- cmode$dseas
        logl <- cmode$logl
        sigma2 <- cmode$sigma2
        muHat <- cmode$muHat
        H <- cmode$H
        Hs <- cmode$Hs
        alpha <- cmode$alpha
        alphas <- cmode$alphas
        
        if (length(muHat) == 0) 
            muHat <- 0
        H <- cmode$H
        Hs <- cmode$Hs
        
        if (R == 1) {
            if (length(seed) == 1) 
                set.seed(seed)
            if (pred) {
                indexes <- sample(1:n, replace = TRUE)
                innov <- cmode$residuals[indexes]
            } else innov <- NULL
            z <- arfima.sim(n, model = list(phi = phi, theta = theta, dint = dint, dfrac = dfrac, 
                H = H, alpha = alpha, seasonal = list(phi = phiseas, theta = thetaseas, 
                  dint = dseas, period = period, dfrac = dfs, H = Hs, alpha = alphas)), 
                sigma2 = sigma2, useC = 3, zinit = zinit, muHat = muHat, innov = innov)
            if (!is.null(obj$tsp)) 
                z <- ts(z, start = obj$tsp[1], frequency = obj$tsp[3])
        } else {
            z <- matrix(0, nrow = n, ncol = R)
            if (length(seed) == 1) 
                set.seed(seed)
            for (i in 1:R) {
                if (length(seed) == R) 
                  set.seed(seed[i])
                if (pred) {
                  indexes <- sample(1:n, replace = TRUE)
                  innov <- cmode$residuals[indexes]
                  
                } else innov <- NULL
                z[, i] <- arfima.sim(n, model = list(phi = phi, theta = theta, dint = dint, 
                  dfrac = dfrac, H = H, alpha = alpha, seasonal = list(phi = phiseas, theta = thetaseas, 
                    dint = dseas, period = period, dfrac = dfs, H = Hs, alpha = alphas)), 
                  sigma2 = sigma2, useC = 3, zinit = zinit, muHat = muHat, innov = innov)
            }
        }
        res[[i]] <- z
    }
    res
}



Boot.ARFIMA <- function(obj, dint, dseas, period, R = 1, n.ahead = NULL, n = 1, zinit = NULL, 
    lastpoint = 0, pred = FALSE, seed = NA, ...) {
    
    if (length(seed) > 1 && length(seed) != R) {
        warning("set seed length is not equal to R: using only first seed")
        seed <- seed[1]
    }
    if (is.null(n.ahead)) 
        n.ahead <- n
    phi <- obj$phi
    theta <- obj$theta
    phiseas <- obj$phiseas
    thetaseas <- obj$thetaseas
    dfrac <- obj$dfrac
    dfs <- obj$dfs
    logl <- obj$logl
    sigma2 <- obj$sigma2
    muHat <- if (!pred) 
        obj$muHat else 0
    H <- obj$H
    Hs <- obj$Hs
    H <- obj$H
    Hs <- obj$Hs
    
    if (R == 1) {
        if (length(seed) == 1 && !is.na(seed)) 
            set.seed(seed)
        if (pred) {
            indexes <- sample(1:n, replace = TRUE)
            innov <- obj$residuals[indexes]
        } else innov <- NULL
        if (!is.null(innov)) 
            sigma2 <- 1
        z <- arfima.sim(n.ahead, model = list(phi = phi, theta = theta, dint = dint, dfrac = dfrac, 
            H = H, seasonal = list(phi = phiseas, theta = thetaseas, dint = dseas, period = period, 
                dfrac = dfs, H = Hs)), sigma2 = sigma2, useC = 3, zinit = zinit, muHat = muHat, 
            innov = innov) + lastpoint
        
    } else {
        z <- matrix(0, nrow = n.ahead, ncol = R)
        if (length(seed) == 1 && !is.na(seed)) 
            set.seed(seed)
        for (i in 1:R) {
            if (length(seed) == R) 
                set.seed(seed[i])
            if (pred) {
                indexes <- sample(1:n, replace = TRUE)
                innov <- obj$residuals[indexes]
            } else innov <- NULL
            if (!is.null(innov)) 
                sigma2 <- 1
            z[, i] <- arfima.sim(n.ahead, model = list(phi = phi, theta = theta, dint = dint, 
                dfrac = dfrac, H = H, seasonal = list(phi = phiseas, theta = thetaseas, 
                  dint = dseas, period = period, dfrac = dfs, H = Hs)), sigma2 = sigma2, 
                useC = 3, zinit = zinit, muHat = muHat, innov = innov) + lastpoint
        }
    }
    
    z
} 
