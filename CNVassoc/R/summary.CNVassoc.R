summary.CNVassoc <-
function (object, ref = 1, alpha = 0.05, ...)
{
    x <- object
    nCov <- attr(x, "nCov")
    z <- qnorm(1 - (alpha/2))
    F <- qr.solve(-x$hessian)
    family <- attr(x, "family")
    odds.ratio <- family!="gaussian"
    fun <- if (odds.ratio) exp else function(temp) temp
    if (attr(x, "model") == 1) {
        beta.ref <- x$coefficients[1, ref]
        beta <- c(0, x$coefficients[1, -ref] - beta.ref)
        cc <- NCOL(x$coefficients)
        se <- rep(NA, cc)
        for (i in 1:cc) {
            if (i != ref) se[i] <- sqrt(F[i, i] + F[ref, ref] - 2 * F[ref, i])
        }
        se <- se[!is.na(se)]
        se <- c(NA, se)
        stat <- beta/se
        pvalue <- 2 * (1 - pnorm(abs(stat)))
        lim.inf <- beta - z * se
        lim.sup <- beta + z * se
        coeff <- cbind(fun(cbind(beta, lim.inf, lim.sup)), se, stat, pvalue)
        rownames(coeff) <- paste("CNV", c(ref - 1, (1:cc)[-ref] - 1), sep = "")
    }
    else {
        beta <- x$coefficients[1:2, 1]
        se <- sqrt(diag(F)[1:2])
        stat <- beta/se
        pvalue <- 2 * (1 - pnorm(abs(stat)))
        lim.inf <- beta - z * se
        lim.sup <- beta + z * se
        coeff <- cbind(fun(cbind(beta, lim.inf, lim.sup)), se, stat, pvalue)
        rownames(coeff) <- c("(Intercept)", "trend")
        cc <- 2
    }
    if (nCov > 0) {
        se.cov <- rep(NA, nCov)
        ii <- 2 + as.integer(attr(x, "model") != 1)
        beta.cov <- x$coefficients[ii:NROW(x$coefficients)]
        cov.names <- rownames(x$coefficients)[ii:NROW(x$coefficients)]
        for (i in 1:nCov) {
            se.cov[i] <- sqrt(F[cc + i, cc + i])
        }
        lim.inf <- beta.cov - z * se.cov
        lim.sup <- beta.cov + z * se.cov
        stat <- beta.cov/se.cov
        pvalue <- 2 * (1 - pnorm(abs(stat)))
        coeff.cov <- cbind(fun(cbind(beta.cov, lim.inf, lim.sup)),
            se.cov, stat, pvalue)
    }
    else {
        coeff.cov <- NULL
        cov.names <- NULL
    }
    coeff <- rbind(coeff, coeff.cov)
    cn <- c("lower.lim", "upper.lim", "SE", "stat", "pvalue")
    if (family == "binomial") 
        cn <- c("OR", cn)
    if (family == "gaussian")
        cn <- c("beta", cn)
    if (family == "poisson")
        cn <- c("RR", cn)
    if (family == "weibull")
        cn <- c("HR", cn)
    rownames(coeff) <- c(rownames(coeff)[1:cc], cov.names)
    colnames(coeff) <- cn
    if (odds.ratio & attr(x, "model") != 1)
        coeff <- coeff[-1, , drop = FALSE]
    out <- x
    out$coefficients <- coeff
    loglike <- logLik(x)
    names(loglike) <- NULL
    out$deviance <- -2 * loglike[1]
    out$numparam <- loglike[2]
    class(out) <- c("summary.CNVassoc", "CNVassoc")
    out
}

