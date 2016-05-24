riv_sest <- function(Y, Xend, Xex, Zinst, intercept) {
    if (is.null(Xex)) {
        X <- Xend
        Z <- cbind(Xend, Zinst, Y)
    } else {
        X <- cbind(Xend, Xex)
        Z <- cbind(Xend, Zinst, Xex, Y)
    }
    
    n <- length(Y)
    p <- ncol(X)
    k <- ncol(Zinst)
    r <- k + p + 1

    kend <- ncol(Xend)
    
    res1 <- CovSest(Z, method='bisquare')
    L <- res1@center
    names(L) <- colnames(Z)
    V <- res1@cov
    
    c <- Tbsc(.5, r)
    
    Econst <- IFconst(r, c)
    
    res2 <- IF.Sest(Z, L, V, c, Econst)
    IFL <- res2$IFL
    IFV <- res2$IFV

    IFriv <- IF.RIV(L, V, IFL, IFV, kend, k, n, intercept)
    
    AV.IFriv <- AVif(IFriv)

    # Standard errors and inference
    df.riv <- n - p - intercept
    sd.riv <- sqrt(diag(AV.IFriv))

    betaiv <- betaIV(L, V, X, Y, kend, k, intercept)
    biv <- betaiv$beta
    tval <- biv/sd.riv
    pv <- 2 * (1 - pt(abs(tval), df.riv))

    resid <- betaiv$resid
    sigma.hat1 <- as.numeric(crossprod(resid)/df.riv)
    sigma.hat3 <- (mad(resid))^2
    tabRIV <- cbind(biv, sd.riv, tval, pv)
    colnames(tabRIV) <- c('Coef', 'Std.Err.', 't', 'p.values')
    MD <- mahalanobis(Z, L, V)

    result <- list(Summary.Table = tabRIV,
                   VC = AV.IFriv,
                   MSE = c(sigma.hat1, sigma.hat3),
                   MD = MD)

    if (k == kend) {
      weight <- psi.bisquare(sqrt(MD), c)
      weight <- weight/sum(weight)
      result[['weight']] <- weight

      w.resid <- resid * weight
      sigma.hat2 <- as.numeric(crossprod(w.resid)/df.riv)
      result[['MSE']] <- c(sigma.hat1, sigma.hat2, sigma.hat3)
    }
    result
}
