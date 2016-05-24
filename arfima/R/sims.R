arfima.sim.work <- function(n, phi = numeric(0), theta = numeric(0), dfrac = numeric(0), 
    H = numeric(0), alpha = numeric(0), phiseas = numeric(0), thetaseas = numeric(0), dfs = numeric(0), 
    Hs = numeric(0), alphas = numeric(0), period = 0, useC = T, useCt = T, sigma2 = 1, rand.gen = rnorm, 
    innov = NULL, ...) {
    
    r <- tacvfARFIMA(phi = phi, theta = theta, dfrac = dfrac, H = H, alpha = alpha, phiseas = phiseas, 
        thetaseas = thetaseas, dfs = dfs, Hs = Hs, alphas = alphas, period = period, maxlag = n - 
            1, useCt = useCt, sigma2 = sigma2)
    
    if (length(innov) == 0) 
        z <- DLSimulate(n, r, rand.gen = rand.gen, ...) else z <- DLSimForNonPar(n, innov, r, dint = 0, dseas = 0, period = 0, muHat = 0, zinit = NULL)  #set to 0 and null because arfima.sim already adds mean and integrates.
    
    return(z)
}

DLSimForNonPar <- function(n, a, r, dint = 0, dseas = 0, period = 0, muHat = 0, zinit = NULL) {
    
    
    nr <- length(r)
    if (min(n, nr) <= 0) 
        stop("input error")
    if (nr < n) 
        r <- c(r, rep(0, n - nr))
    r <- r[1:n]
    EPS <- .Machine$double.eps  # 1+EPS==1, machine epsilon
    z <- numeric(n)
    out <- .C("durlevsim", z = as.double(z), as.double(a), as.integer(n), as.double(r), 
        as.double(EPS), fault = as.integer(1))
    fault <- out$fault
    if (fault == 1) 
        stop("error: sequence not p.d.")
    z <- out$z + muHat
    
    z
}


arfima.sim <- function(n, model = list(phi = numeric(0), theta = numeric(0), dint = 0, dfrac = numeric(0), 
    H = numeric(0), alpha = numeric(0), seasonal = list(phi = numeric(0), theta = numeric(0), 
        dint = 0, period = numeric(0), dfrac = numeric(0), H = numeric(0), alpha = numeric(0))), 
    useC = 3, sigma2 = 1, rand.gen = rnorm, muHat = 0, zinit = NULL, innov = NULL, ...) {
    
    H <- model$H
    dfrac <- model$dfrac
    dint <- model$dint
    alpha <- model$alpha
    
    if (length(H) == 0) 
        H <- numeric(0)
    if (length(dfrac) == 0) 
        dfrac <- numeric(0)
    if (length(alpha) == 0) 
        alpha <- numeric(0)
    if (length(H) + length(dfrac) + length(alpha) > 1) {
        warning("two or more of dfrac, H and alpha have been specified: using only dfrac")
        H <- numeric(0)
    }
    if (useC == 0) 
        useC <- useCt <- F else if (useC == 1) {
        useC <- T
        useCt <- F
    } else if (useC == 2) {
        useC <- F
        useCt <- T
    } else useC <- useCt <- T
    if (length(dint) == 0) 
        dint <- 0
    if (!is.null(dint) && (round(dint) != dint || dint < 0)) 
        stop("dint must be an integer >= 0")
    
    phi <- model$phi
    theta <- model$theta
    seasonal <- model$seasonal
    if (length(phi) == 0) 
        phi <- numeric(0)
    if (length(theta) == 0) 
        theta <- numeric(0)
    if (!is.null(seasonal)) {
        if ((length(seasonal$period) == 0) || (seasonal$period == 0)) {
            dseas <- 0
            period <- 0
            dfs <- numeric(0)
            phiseas <- numeric(0)
            thetaseas <- numeric(0)
            Hs <- numeric(0)
            alphas <- numeric(0)
        } else {
            dseas <- seasonal$dint
            period <- seasonal$period
            dfs <- seasonal$dfrac
            Hs <- seasonal$H
            alphas <- seasonal$alpha
            
            if (length(period) == 0 || period < 2 || period != round(period)) 
                stop("if using a periodic model, must have an integer period >=2")
            if (length(dseas) == 0) 
                dseas <- 0
            if (round(dseas) != dseas || dseas < 0) 
                stop("seasonal d must be an integer >= 0")
            if (length(dfs) == 0) 
                dfs <- numeric(0)
            if (length(Hs) == 0) 
                Hs <- numeric(0)
            if (length(alphas) == 0) 
                alphas <- numeric(0)
            if (length(Hs) + length(dfs) + length(alphas) > 1) {
                warning("two or more of seasonal dfrac, H and alpha have been specified: using only seasonal dfrac")
                Hs <- numeric(0)
            }
            phiseas <- seasonal$phi
            thetaseas <- seasonal$theta
            if (length(phiseas) == 0) 
                phiseas <- numeric(0)
            if (length(thetaseas) == 0) 
                thetaseas <- numeric(0)
        }
        
    } else {
        seasonal <- NULL
        dseas <- 0
        period <- 0
        dfs <- numeric(0)
        phiseas <- numeric(0)
        thetaseas <- numeric(0)
        Hs <- numeric(0)
        alphas <- numeric(0)
    }
    
    if (!IdentInvertQ(phi = phi, theta = theta, dfrac = dfrac, H = H, alpha = alpha, phiseas = phiseas, 
        thetaseas = thetaseas, dfs = dfs, Hs = Hs, alphas = alphas, period = period, ident = TRUE)) {
        stop("Model is non-causal, non-invertible or non-identifiable")
    }
    
    
    z <- arfima.sim.work(n, phi = phi, theta = theta, dfrac = dfrac, H = H, alpha = alpha, 
        phiseas = phiseas, thetaseas = thetaseas, dfs = dfs, Hs = Hs, alphas = alphas, period = period, 
        useC = useC, useCt = useCt, sigma2 = sigma2, rand.gen = rand.gen, innov = innov, 
        ...)
    z <- z + muHat
    if (dint + dseas > 0) {
        icap <- dint + dseas * period
        if (is.null(zinit)) 
            zinit <- rep(0, icap)
        z <- integ(z, zinit, dint, dseas, period)[-c(1:icap)]
    }
    as.ts(z)
} 
