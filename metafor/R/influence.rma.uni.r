influence.rma.uni <-
function (model, digits, ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    if (missing(digits)) 
        digits <- x$digits
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    tau2.del <- rep(NA_real_, x$k.f)
    delpred <- rep(NA_real_, x$k.f)
    vdelpred <- rep(NA_real_, x$k.f)
    QE.del <- rep(NA_real_, x$k.f)
    dffits <- rep(NA_real_, x$k.f)
    dfbetas <- matrix(NA_real_, nrow = x$k.f, ncol = length(x$b))
    cook.d <- rep(NA_real_, x$k.f)
    cov.r <- rep(NA_real_, x$k.f)
    weight <- rep(NA_real_, x$k.f)
    pred.full <- x$X.f %*% x$b
    if (x$weighted) {
        if (is.null(x$weights)) {
            wi <- 1/(x$vi + x$tau2)
            W <- diag(wi, nrow = x$k, ncol = x$k)
            svb <- crossprod(x$X, W) %*% x$X/x$s2w
        }
        else {
            svb <- chol2inv(chol(x$vb))
            A <- diag(x$weights, nrow = x$k, ncol = x$k)
            stXAX <- .invcalc(X = x$X, W = A, k = x$k)
            H <- x$X %*% stXAX %*% t(x$X) %*% A
        }
    }
    else {
        svb <- chol2inv(chol(x$vb))
        stXX <- .invcalc(X = x$X, W = diag(x$k), k = x$k)
        H <- x$X %*% stXX %*% t(x$X)
    }
    options(na.action = "na.pass")
    hii <- hatvalues(x)
    options(na.action = na.act)
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq_len(x$k.f)[x$not.na]) {
        res <- try(suppressWarnings(rma(x$yi.f[-i], x$vi.f[-i], 
            weights = x$weights.f[-i], mods = cbind(x$X.f[-i, 
                ]), method = x$method, weighted = x$weighted, 
            intercept = FALSE, knha = x$knha, control = x$control)), 
            silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        if (any(res$coef.na)) 
            next
        tau2.del[i] <- res$tau2
        QE.del[i] <- res$QE
        Xi <- matrix(x$X.f[i, ], nrow = 1)
        delpred[i] <- Xi %*% res$b
        vdelpred[i] <- Xi %*% tcrossprod(res$vb, Xi)
        if (x$weighted) {
            if (is.null(x$weights)) {
                dffits[i] <- (pred.full[i] - delpred[i])/sqrt(res$s2w * 
                  hii[i] * (tau2.del[i] + x$vi.f[i]))
            }
            else {
                dffits[i] <- (pred.full[i] - delpred[i])/sqrt(res$s2w * 
                  diag(H %*% diag(tau2.del[i] + x$vi, nrow = x$k, 
                    ncol = x$k) %*% t(H)))[i - sum(!x$not.na[1:i])]
            }
        }
        else {
            dffits[i] <- (pred.full[i] - delpred[i])/sqrt(res$s2w * 
                diag(H %*% diag(tau2.del[i] + x$vi, nrow = x$k, 
                  ncol = x$k) %*% t(H)))[i - sum(!x$not.na[1:i])]
        }
        dfbeta <- x$b - res$b
        if (x$weighted) {
            if (is.null(x$weights)) {
                vb.del <- .invcalc(X = x$X, W = diag(1/(x$vi + 
                  tau2.del[i]), nrow = x$k, ncol = x$k), k = x$k)
            }
            else {
                vb.del <- tcrossprod(stXAX, x$X) %*% A %*% diag(x$vi + 
                  tau2.del[i], nrow = x$k, ncol = x$k) %*% A %*% 
                  x$X %*% stXAX
            }
        }
        else {
            vb.del <- tcrossprod(stXX, x$X) %*% diag(x$vi + tau2.del[i], 
                nrow = x$k, ncol = x$k) %*% x$X %*% stXX
        }
        dfbetas[i, ] <- dfbeta/sqrt(res$s2w * diag(vb.del))
        cook.d[i] <- (crossprod(dfbeta, svb) %*% dfbeta)
        cov.r[i] <- det(res$vb)/det(x$vb)
    }
    delresid <- x$yi.f - delpred
    sedelresid <- sqrt(x$vi.f + vdelpred + tau2.del)
    standelres <- delresid/sedelresid
    options(na.action = "na.omit")
    weight[x$not.na] <- weights(x)
    options(na.action = na.act)
    inf <- cbind(standelres, dffits, cook.d, cov.r, tau2.del, 
        QE.del, hii, weight)
    dfb <- cbind(dfbetas)
    inf <- data.frame(inf)
    dfb <- data.frame(dfb)
    is.infl <- abs(inf$dffits) > 3 * sqrt(x$p/(x$k - x$p)) | 
        pchisq(inf$cook.d, df = x$p) > 0.5 | inf$hii > 3 * x$p/x$k | 
        apply(abs(dfb) > 1, 1, any)
    out <- list(inf = inf, dfb = dfb, tau2 = x$tau2, QE = x$QE, 
        ids = x$ids, not.na = x$not.na, is.infl = is.infl, k = x$k, 
        p = x$p, digits = digits)
    rownames(out$inf) <- x$slab
    rownames(out$dfb) <- x$slab
    colnames(out$dfb) <- rownames(x$b)
    colnames(out$inf) <- c("rstudent", "dffits", "cook.d", "cov.r", 
        "tau2.del", "QE.del", "hat", "weight")
    class(out) <- "infl.rma.uni"
    return(out)
}
