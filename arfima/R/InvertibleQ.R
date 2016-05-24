InvertibleQ <- function(phi) {
    if (length(phi) == 0) 
        return(T)
    return(all(abs(ARToPacf(phi)) < 1))
}

InvertibleD <- function(d) {
    if (length(d) == 0) 
        return(T)
    return(d > -1 && d < 0.5)
}

InvertibleH <- function(H) {
    if (length(H) == 0) 
        return(T)
    return(H > 0 && H < 1)
}

InvertibleAlpha <- function(alpha) {
    if (length(alpha) == 0) 
        return(T)
    return(alpha > 0 && alpha < 3)
}

IdentifiableQ <- function(phi = numeric(0), theta = numeric(0)) {
    p <- length(phi)
    q <- length(theta)
    k <- p + q
    A <- matrix(0, nrow = k, ncol = k)
    if (p > 0) 
        for (j in 1:p) {
            for (i in 1:k) {
                if (i - j == 0) 
                  A[i, j] <- -1
                if (i - j <= q && i - j > 0) 
                  A[i, j] <- theta[i - j]
            }
        }
    if (q > 0) 
        for (j in 1:q) {
            for (i in 1:k) {
                if (i - j == 0) 
                  A[i, j + p] <- 1
                
                if (i - j <= p && i - j > 0) 
                  A[i, j + p] <- -phi[i - j]
            }
        }
    return(det(A) > 0)
}

IdentifiableQQ <- function(phi = numeric(0), theta = numeric(0), phiseas = numeric(0), thetaseas = numeric(0), 
    negphi = T) {
    p <- length(phi)
    q <- length(theta)
    ps <- length(phiseas)
    qs <- length(thetaseas)
    
    k <- p + ps + q + qs
    A <- matrix(0, nrow = k, ncol = k)
    if (p > 0) 
        for (j in 1:p) {
            for (i in 1:k) {
                if (i - j == 0) 
                  A[i, j] <- -1
                if (i - j <= q && i - j > 0) 
                  A[i, j] <- theta[i - j]
            }
        }
    if (q > 0) 
        for (j in 1:q) {
            for (i in 1:k) {
                if (i - j == 0) 
                  A[i, j + p] <- if (negphi) 
                    1 else -1
                if (i - j <= p && i - j > 0) 
                  A[i, j + p] <- if (negphi) 
                    -phi[i - j] else phi[i - j]
            }
        }
    if (ps > 0) 
        for (j in 1:ps) {
            for (i in 1:k) {
                if (i - j == 0) 
                  A[i + p + q, j + p + q] <- -1
                if (i - j <= qs && i - j > 0) 
                  A[i + p + q, j + p + q] <- thetaseas[i - j]
            }
        }
    if (qs > 0) 
        for (j in 1:qs) {
            for (i in 1:k) {
                if (i - j == 0) 
                  A[i + p + q, j + p + q + ps] <- if (negphi) 
                    1 else -1
                if (i - j <= ps && i - j > 0) 
                  A[i + p + q, j + p + q + ps] <- if (negphi) 
                    -phiseas[i - j] else phiseas[i - j]
            }
        }
    return(det(A) > 0)
}

IdentInvertQ <- function(phi = numeric(0), theta = numeric(0), phiseas = numeric(0), thetaseas = numeric(0), 
    dfrac = numeric(0), dfs = numeric(0), H = numeric(0), Hs = numeric(0), alpha = numeric(0), 
    alphas = numeric(0), delta = numeric(0), period = 0, debug = FALSE, ident = TRUE) {
    if (!(InvertibleQ(phi) && InvertibleQ(theta) && InvertibleQ(phiseas) && InvertibleQ(thetaseas) && 
        InvertibleD(dfrac) && InvertibleD(dfs) && InvertibleH(H) && InvertibleH(Hs) && InvertibleQ(delta) && 
        InvertibleAlpha(alpha) && InvertibleAlpha(alphas))) {
        if (debug) 
            warning("Model is non-stationary or non-invertible")
        return(FALSE)
    }
    if ((length(H) + length(dfrac) + length(alpha) > 1) || (length(Hs) + length(dfs) + length(alphas) > 
        1)) {
        if (debug) 
            warning("Model is contains fractional d or H in either seasonal or non-seasonal components")
        return(FALSE)
    }
    if (ident && (length(phi) > 0 || length(theta) > 0 || length(phiseas) > 0 || length(thetaseas) > 
        0 || length(dfrac) > 0 || length(dfs) > 0)) {
        if (length(dfrac) > 0) 
            d <- T else d <- F
        if (length(dfs) > 0) 
            dd <- T else dd <- F
        I <- iARFIMA(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, 
            dfrac = d, dfs = dd, period = period)
        if (!recdet(I)) {
            if (debug) 
                warning("Information matrix of SARMA or ARFIMA process is not positive definite")
            return(FALSE)
        }
    }
    
    TRUE
}


recdet = function(a) {
    p <- sqrt(length(a))
    if (p != round(p)) 
        stop("not a square matrix")
    a <- as.matrix(a, nrow = p)
    for (i in 1:p) {
        if (det(as.matrix(a[1:i, 1:i], nrow = i)) <= 0) 
            return(F)
    }
    return(T)
} 
