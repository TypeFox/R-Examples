summary.ltm <-
function (object, robust.se = FALSE, ...) {
    if (!inherits(object, "ltm"))
        stop("Use only with 'ltm' objects.\n")
    if (object$IRT.param)
        irt <- IRT.parm(object, TRUE, robust = robust.se)
    coefs <- if (object$IRT.param) irt$parms else object$coef
    Var.betas <- vcov(object, robust = robust.se)
    se <- if (!is.null(constraint <- object$constraint)) {
        p <- nrow(coefs)
        res <- matrix(NA, p, ncol(coefs))
        res[-((constraint[, 2] - 1) * p + constraint[, 1])] <- if (object$IRT.param) irt$se else sqrt(diag(Var.betas))
        c(res)
    } else {
        if (object$IRT.param) irt$se else sqrt(diag(Var.betas))
    }
    z.vals <- c(coefs)/se
    coef.tab <- cbind(value = c(coefs), std.err = se, z.vals = z.vals)
    p <- nrow(coefs)
    rownames(coef.tab) <- if (object$IRT & (object$ltst$factors == 1 & !object$ltst$quad.z1))
            paste(rep(colnames(coefs), each = p), rownames(coefs), sep = ".")
        else
            paste(rep(object$ltst$nams, each = p), rownames(coefs), sep = ".")
    out <- list(coefficients = coef.tab, Var.betas = Var.betas)
    out$logLik <- object$log.Lik
    out$AIC <- -2 * object$log.Lik + 2 * length(coefs)
    out$BIC <- -2 * object$log.Lik + log(nrow(object$X)) * length(coefs)
    out$max.sc <- object$max.sc
    out$conv <- object$conv
    out$counts <- object$counts
    out$call <- object$call
    out$ltst <- object$ltst
    out$control <- object$control
    out$attr <- attr(object$X, "items")
    out$ancr <- attr(object$X, "anchoring")
    class(out) <- "summ.ltm"
    out
}
