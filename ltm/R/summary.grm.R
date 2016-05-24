summary.grm <-
function (object, ...) {
    if (!inherits(object, "grm"))
        stop("Use only with 'grm' objects.\n")
    coefs <- if (object$IRT.param) IRT.parm(object, digits.abbrv = object$control$digits)$parms else object$coefficients
    coef.tab <- if (!is.null(object$hessian)) {
        p <- length(coefs)
        coefs. <- coefs
        ncatg <- sapply(coefs, length)
        if (object$constrained)
            coefs.[-p] <- lapply(coefs[-p], function (x) x[-length(x)])
        coefs. <- unlist(coefs., use.names = FALSE)
        if (all(!is.na(object$hessian) & is.finite(object$hessian))) {
            ev <- eigen(object$hessian, TRUE, TRUE)$values
            if (!all(ev >= -1e-06 * abs(ev[1]))) 
                warning("Hessian matrix at convergence is not positive definite; unstable solution.\n")
        } else 
            stop("Hessian matrix at convergence contains infinite or missing values; unstable solution.\n")
        se <- if (object$IRT.param) IRT.parm(object, TRUE, digits.abbrv = object$control$digits)$se else sqrt(diag(vcov(object)))
        z.vals <- coefs.
        for (i in seq(along = z.vals))
            z.vals[i] <- coefs.[i] / se[i]
        d <- data.frame(value = coefs., std.err = se, z.vals = z.vals)
        d <- if (object$constrained)
                split(d[c(seq(1, nrow(d) - 1), rep(nrow(d), p)), ], c(rep(1:p, ncatg - 1), 1:p)) 
            else
                split(d, rep(1:p, ncatg))
        names(d) <- names(coefs)
        for (i in seq(along = coefs))
            rownames(d[[i]]) <- names(coefs[[i]])
        d
        } else {
        coefs <- lapply(coefs, function (x) cbind(value = x))
    }
    out <- list(coefficients = coef.tab)
    out$logLik <- object$log.Lik
    df <- sapply(object$coef, length)
    df <- if (object$constrained) sum(df) - length(df) + 1 else sum(df)
    out$AIC <- -2 * object$log.Lik + 2 * df
    out$BIC <- -2 * object$log.Lik + log(nrow(object$X)) * df
    out$max.sc <- object$max.sc
    out$conv <- object$conv
    out$counts <- object$counts
    out$call <- object$call
    out$control <- object$control
    out$attr <- attr(object$X, "items")
    out$ancr <- attr(object$X, "anchoring")
    class(out) <- "summ.grm"
    out
}
