`summary.cusp` <-
function (object, correlation = FALSE, symbolic.cor = FALSE, 
    logist=FALSE, ...) 
{
    est.disp <- FALSE
    df.r <- object$df.residual
#    dispersion <- if (df.r > 0) {
#        est.disp <- TRUE
#        if (any(object$weights == 0)) 
#            warning("observations with zero weight not used for calculating dispersion")
#        sum((object$weights * object$residuals^2)[object$weights > 
#            0])/df.r
#    }
#    else {
#        est.disp <- TRUE
#        NaN
#    }
    dispersion <- 1 # no least squares fit but full ML estimate
    aliased <- is.na(coef(object))
    p <- object$rank
    if (p > 0) {
        p1 <- 1:p
        Qr <- object$qr
        coef.p <- object$coefficients[!aliased][Qr$pivot[p1]]
        covmat.unscaled <- solve(qr.X(Qr)[p1, p1, drop=FALSE]) #chol2inv(Qr$qr[p1, p1, drop = FALSE])
        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
        nms = names(object$coefficients)[!aliased]
        covmat.unscaled <- covmat.unscaled[nms, nms]
        coef.p <- coef.p[nms]
        covmat <- dispersion * covmat.unscaled
        var.cf <- diag(covmat)
        s.err <- sqrt(var.cf)
        tvalue <- coef.p/s.err
        dn <- c("Estimate", "Std. Error")
        if (!est.disp) {
            pvalue <- 2 * pnorm(-abs(tvalue))
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p), c(dn, 
                "z value", "Pr(>|z|)"))
        }
        else if (df.r > 0) {
            pvalue <- 2 * pt(-abs(tvalue), df.r)
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p), c(dn, 
                "t value", "Pr(>|t|)"))
        }
        else {
            coef.table <- cbind(coef.p, NaN, NaN, NaN)
            dimnames(coef.table) <- list(names(coef.p), c(dn, 
                "t value", "Pr(>|t|)"))
        }
        df.f <- length(aliased)
    }
    else {
        coef.table <- matrix(, 0, 4)
        dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error", 
            "t value", "Pr(>|t|)"))
        covmat.unscaled <- covmat <- matrix(, 0, 0)
        df.f <- length(aliased)
    }
    fnobs <- sum(object$weights)

    # cusp criterion values
    wmeanerr  <- sum((object$weights * object$residuals)[object$weights > 0]) / fnobs
    wvarerr   <- sum((object$weights * object$residuals^2)[object$weights > 0]) / fnobs - wmeanerr^2
    r2cusp     <- 1 - wvarerr / var(object$y) #var(model.response(object$model))
    r2cusp.aic <- AIC(object)
    K <- attr(logLik(object),'df')
    r2cusp.logLik = logLik(object)[1];
    r2cusp.npar = K
    r2cusp.aicc <- r2cusp.aic + 2 * K * (K+1) / (fnobs - K - 1)
    r2cusp.bic <- r2cusp.aic - 2*K + K*log(fnobs) #AIC(object, k=log(fnobs))
    r2cusp.df  <- object$df.residual
    r2cusp.dev <- deviance(object)
#    mf <- object$call
#    m <- match(names(formals(lm)), names(mf), 0)
#    mf <- mf[c(1,m)]
#    mf[[1]] <- as.name("lm")
#    fit.lm <- eval(mf, parent.frame());
#    fit.lm.sum <- summary(fit.lm)
#    r2lin <- fit.lm.sum$r.squared
#    r2lin.df <- fit.lm.sum$df[2]
#    r2lin.aic <- AIC(fit.lm)
#    K <- attr(logLik(fit.lm),'df')
#    r2lin.aicc <- r2lin.aic + 2 * K * (K+1) / (fnobs - K - 1)
#    r2lin.bic <- AIC(fit.lm, k=log(NROW(object$model)))
#    r2lin.dev <- deviance(fit.lm)

    # linear model criterion values
#    browser()
    fit.lm <- cancor(object$x$X.y,cbind(object$x$X.alpha,object$x$X.beta),xcenter=FALSE,ycenter=FALSE)
    idx = which(zapsmall(fit.lm$cor)<1)[1]
    r2lin <- fit.lm$cor[idx]^2 # first canonical correlation smaller tnan 1
    K <- qr(cbind(object$x$X.y, object$x$X.alpha, object$x$X.beta))$rank #(NROW(fit.lm$xcoef) + NROW(fit.lm$ycoef) - 1) # one constraint is imposed (not 2, it's not cannonical correlation; cancor is just a way of computing)
    r2lin.df <- fnobs - K
    rss.lm <- cusp.subspacerss(dependents = object$x$X.y, predictors = cbind(object$x$X.alpha,object$x$X.beta))
    idx = which(zapsmall(rss.lm$cor)<1)[1]
#    r2lin.rss <- (1-r2lin) * sum((crossprod(object$x$X.y) %*% fit.lm$xcoef[,idx])^2) 
    r2lin.rss <- rss.lm$rss[idx]
    r2lin.logLik <- -0.5*fnobs*(log(2*pi*sum(r2lin.rss)/fnobs)+1)
    r2lin.npar <- K
    r2lin.aic <- -2*r2lin.logLik + 2 * K 
    r2lin.aicc <- r2lin.aic + 2 * K * (K+1) / (fnobs - K - 1)
    r2lin.bic <- -2*r2lin.logLik + K*log(fnobs)
    r2lin.dev <- r2lin.rss

    # logist model criterion values
#browser()
    mf <- object$call
    mf$start <- NULL # quick fix (better to select only arguments known to cusp.logist? extend logist with start?)
    mf[[1]] <- as.name("cusp.logist")
    mf$hessian <- FALSE
    fit.logist <- if(logist) {eval(mf, parent.frame())} else {list()} # glm(object$model[,1] ~ -1);
    r2log <- if(logist) {fit.logist$rsq} else {NA} #1 - var(fit.logist$residuals) / var(fit.logist$y) else  NA
    r2log.df <- if(logist) {fit.logist$df.resid} else{ NA}
    r2log.aic <- if(logist) {fit.logist$aic} else{ NA}
    K <- if(logist) {fnobs - fit.logist$df.resid } else {NA}
    r2log.logLik <- if(logist) {fit.logist$logLik} else {NA}
    r2log.npar <- if(logist) {K} else {NA}
    r2log.aicc <- if(logist){ r2log.aic + 2 * K * (K+1) / (fnobs - K - 1)} else { NA }
    r2log.bic <- if(logist) {r2log.aic - 2*K + K*log(fnobs)} else {NA}# AIC(fit.logist, k=log(NROW(fit.logist$model)))
    r2log.dev <- if(logist) {fit.logist$rss} else {NA} #deviance(fit.logist)
    
    dispers <- sum((object$weights * object$residuals^2)[object$weights > 0]) / df.r
    ans <- c(object[c("call", "terms", "family", "deviance", 
        "aic", "contrasts", "df.residual", "null.deviance", "df.null", 
        "iter", "na.action")], list(deviance.resid = residuals(object, 
        type = "deviance"), coefficients = coef.table, aliased = aliased, 
        dispersion = dispersion, df = c(object$rank, df.r, df.f), resid.name = colnames(object$residual)[1], 
        cov.unscaled = covmat.unscaled, cov.scaled = covmat),
        `r2lin`  = list(r.squared = r2lin,  dev=r2lin.dev,  df=r2lin.df,  logLik=r2lin.logLik,  npar=r2lin.npar,  aic=r2lin.aic,  aicc=r2lin.aicc, bic=r2lin.bic),
        `r2log`  = list(r.squared = r2log,  dev=r2log.dev,  df=r2log.df,  logLik=r2log.logLik,  npar=r2log.npar,  aic=r2log.aic,  aicc=r2log.aicc,  bic=r2log.bic),
        `r2cusp` = list(r.squared = r2cusp, dev=r2cusp.dev, df=r2cusp.df, logLik=r2cusp.logLik, npar=r2cusp.npar, aic=r2cusp.aic, aicc=r2cusp.aicc, bic=r2cusp.bic)
    )
    if (correlation && p > 0) {
        dd <- sqrt(diag(covmat.unscaled))
        ans$correlation <- covmat.unscaled/outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- c("summary.cusp","summary.glm")
    return(ans)
}

