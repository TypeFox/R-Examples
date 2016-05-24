`cusp.logist.fit` <-
function(Xa, Xb, Y, hessian=FALSE, use = 'new', nfits = 5, ndigit = 17, max.repeats = 4, 
    gtol = .Machine$double.eps^(1/3), ...)
{
    qxa <- qr(Xa)
    qxb <- qr(Xb)
    qy  <- if(is.qr(Y)) {Y} else {qr(Y)}
    
    rxa  <- qxa$rank; ixa <- 1:rxa
    rxb  <- qxb$rank; ixb <- 1:rxb
    ry   <- qy$rank;  iy  <- 1:ry
    Qxa <- qr.Q(qxa)[,ixa, drop = FALSE]
    Qxb <- qr.Q(qxb)[,ixb, drop = FALSE]
    Qy  <- qr.Q(qy)[,iy, drop = FALSE]
    Rxa <- qr.R(qxa)[ixa,ixa]
    Rxb <- qr.R(qxb)[ixb,ixb]
    Ry  <- qr.R(qy)[iy,iy]
    
    objf <- switch(use,old=cusp.logist.objf.old, new=cusp.logist.objf)
    fitq <- nlm(objf,rnorm(rxa+rxb-1), xa=Qxa, xb=Qxb, y=qr(Qy), ndigit=ndigit, ...)
    j=1
    repeat{
        for(i in 1:nfits){
            tmp <- nlm(objf,rnorm(rxa+rxb-1), xa=Qxa, xb=Qxb, y=Qy, ndigit=ndigit, ...)
            if(tmp$min<fitq$min) {fitq <- tmp}
        }
        if(max(abs(fitq$gradient*ifelse(abs(fitq$est)<1.0, 1.0, fitq$est))) < gtol || j>max.repeats) {
            break;
        }
        j <- j + 1
    }
    rss <- objf(fitq$est, xa=Qxa, xb=Qxb, y=Y)
    rsq <- attr(rss, 'RSq')
    ahat <- rep(NA, ncol(Xa))
    ahat[qxa$pivot[ixa]] <- backsolve(Rxa, fitq$est[ixa]);
    bhat <- rep(NA, ncol(Xb))
    bhat[qxb$pivot[ixb]] <- backsolve(Rxb, c(1,fitq$est[rxa+ixb[-1]-1]))
    if(bhat[1]){
        ahat <- ahat / bhat[1]^2
        bhat <- bhat / bhat[1]
    }
    # coefficients
    coef <- c(ahat, bhat)
    anames <- paste('a.', if(is.null(colnames(Xa))) {1:ncol(Xa)} else {colnames(Xa)}, sep='') 
    bnames <- paste('b.', if(is.null(colnames(Xb))) {1:ncol(Xb)} else {colnames(Xb)}, sep='')
    names(coef) <- c(anames, bnames)
    fitq$coefficients <- coef
    # predicted etc.
    .alpha <- Xa %*% ifelse(is.na(ahat), 0, ahat) # == attr(rss, 'alpha') which is calculated differently!
    .beta  <- Xb %*% ifelse(is.na(bhat), 0, bhat) # == attr(rss, 'beta')
    fitq$linear.predictors <- cbind(alpha=.alpha, beta=.beta)
    colnames(fitq$linear.predictors) <- c('alpha', 'beta')
    rownames(fitq$linear.predictors) <- rownames(Y)
    fitq$fitted.values <- drop((.alpha/.beta^2 > -50.0) / (1 + exp(-.alpha/.beta^2)))
    names(fitq$fitted.values) <- rownames(Y)
    res <- lm(fitq$fitted.values ~ Y-1)$resid
    fitq$residuals <- sqrt(rss) * res/sqrt(sum(res^2))
    names(fitq$residuals) <- rownames(Y)
    # deviance & likelihood
    nobs <- NROW(Y)
    fitq$rank <- qxa$rank + qxb$rank + qy$rank - 1 # one restriction
    fitq$deviance <- c(`Sum Sq. Err.` = rss) # drop attributes
    fitq$logLik <- -0.5*nobs*(1+log(2*pi*sum(fitq$residuals^2)/nobs))
    fitq$aic <- -2*fitq$logLik + 2*fitq$rank
    fitq$rsq <- rsq
    fitq$df.residual <- nobs - fitq$rank
    fitq$df.null <- nobs - 1
    # hessian
    fitq$rss <- c(rss)
    fitq$hessian <- if(hessian) {nlm(objf,fitq$est, xa=Qxa, xb=Qxb, y=Y, hessian=TRUE)$hessian} else {matrix(0,0,0)}
    fitq$Hessian <- "Not implemented"
    fitq$qr <- if(!is.null(fitq$hessian)) {qr(fitq$hessian)} else {"Not implemented"}
    # convergence
    fitq$converged <- fitq$code == 1 || max(abs(fitq$gradient*ifelse(abs(fitq$est)<1.0, 1.0, fitq$est))) < gtol
    fitq$converged <- fitq$converged  && (!hessian || all(eigen(fitq$hessian,,only.values=TRUE)$values>0))
    fitq$weights <- rep(1, NROW(Y))
    fitq$boundary <- rep(NA, length(coef))
    fitq
}

