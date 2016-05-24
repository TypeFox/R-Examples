riv_sdest <- function(Y, Xend, Xex, Zinst, intercept) {
    if (is.null(Xex)) {
        X <- Xend
        Z <- cbind(Xend, Zinst, Y)
    } else {
        X <- cbind(Xend, Xex)
        Z <- cbind(Xend, Zinst, Xex, Y)
    }

    res1 <- CovSde(Z, tune = 0.95, prob = 0.99)
    L <- res1@center
    V <- res1@cov

    n <- length(Y)
    p <- ncol(X)
    k <- ncol(Zinst)
    r <- k + p + 1
    
    kend <- ncol(Xend)

    # Parameter Estimates
    Vm <- matrix(V[(kend + 1):(nrow(V) - 1), -((kend + 1):(kend + k))],
                 nrow = nrow(V) - kend - 1)
    Swx <- Vm[, -ncol(Vm)]
    Sxw <- t(Swx)
    Sww <- matrix(V[(kend + 1):(nrow(V) - 1), (kend + 1):(nrow(V) - 1)],
                  nrow = nrow(V) - kend - 1)
    Swy <- Vm[, ncol(Vm)]
    Lm <- L[-((kend + 1):(kend + k))]
    Mx <- Lm[1:(length(Lm) - 1)]
    My <- Lm[length(Lm)]
    part1 <- Sxw %*% solve(Sww) %*% Swx
    
    b1 <- solve(part1) %*% Sxw %*% solve(Sww) %*% Swy
    b0 <- My - sum(b1 * Mx)
    
    if (intercept) {
        beta.oiv <- matrix(rbind(b0, b1), ncol = 1, dimnames = NULL)
    } else {
        beta.oiv <- b1
    }
        
    # Summary Results
    tabRIV <- beta.oiv
    colnames(tabRIV) <- 'Coef'
    if (intercept)
        rownames(tabRIV) <- c('Intercept', colnames(X))
    else
        rownames(tabRIV) <- colnames(X)
    
    list(Summary.Table = tabRIV)
}    
