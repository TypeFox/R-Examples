riv_classical <- function(Y, Xend, Xex, Zinst, intercept) {
    if (is.null(Xex)) {
        X <- Xend
        Z <- cbind(Xend, Zinst, Y)
    } else {
        X <- cbind(Xend, Xex)
        Z <- cbind(Xend, Zinst, Xex, Y)
    }
    
    L <- apply(Z, 2, mean)
    V <- cov(Z)

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
    
    df.oiv <- n - p - intercept
    
    if (intercept) {
        beta.oiv <- matrix(rbind(b0, b1), ncol = 1, dimnames = NULL)
        resid <- Y - b0 - X %*% b1
        
        if (is.null(Xex)) 
            inst.w <- cbind(1, Zinst)
        else
            inst.w <- cbind(1, Zinst, Xex)
        x.oiv <- cbind(1, X)
    } else {
        beta.oiv <- b1
        resid <- Y - X %*% beta.oiv
        if (is.null(Xex)) 
            inst.w <- Zinst
        else
            inst.w <- cbind(Zinst, Xex)
        x.oiv <- cbind(X)
    }
    
    # Estimate Var-cov of OIV
    zinv <- solve(crossprod(inst.w))
    Pz <- inst.w %*% zinv %*% t(inst.w)
    XPz <- t(x.oiv) %*% Pz %*% x.oiv
    inv.XPz <- solve(XPz)
    
    sigma.hat <- as.numeric(crossprod(resid)/df.oiv)
    var.oiv <- sigma.hat * inv.XPz
    sd.oiv <- sqrt(diag(var.oiv))
    tval <- beta.oiv/sd.oiv
    pv <- 2 * (1 - pt(abs(tval), df.oiv))
    
    # Summary Results
    tabOIV <- cbind(beta.oiv, sd.oiv, tval, pv)
    colnames(tabOIV) <- c('Coef', 'Std.Err.', 't', 'p.values')
    if (intercept)
        rownames(tabOIV) <- colnames(var.oiv) <-
            rownames(var.oiv) <- c('Intercept', colnames(X))
    else
        rownames(tabOIV) <- colnames(var.oiv) <-
            rownames(var.oiv) <- colnames(X)
    
    list(Summary.Table = tabOIV,
         VC = var.oiv,
         MSE = sigma.hat)
}    
