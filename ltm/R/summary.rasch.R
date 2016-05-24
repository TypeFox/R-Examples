summary.rasch <-
function (object, robust.se = FALSE, ...) {
    if (!inherits(object, "rasch"))
        stop("Use only with 'rasch' objects.\n")
    if (object$IRT.param)
        irt <- IRT.parm(object, TRUE, robust = robust.se)
    coefs <- if (object$IRT.param) irt$parms else c(object$coef[, 1], object$coef[1, 2])
    Var.betas <- vcov(object, robust = robust.se)
    se <- rep(NA, length(coefs))
    ind <- if (!is.null(constraint <- object$constraint)) seq(along = se)[-constraint[, 1]] else seq(along = se)
    se[ind] <- if (object$IRT.param) irt$se else sqrt(diag(Var.betas))
    z.vals <- coefs / se
    coef.tab <- cbind(value = coefs, std.err = se, z.vals = z.vals)
    p <- ncol(object$X)
    rownames(coef.tab) <- if (object$IRT) names(coefs) else c(abbreviate(names(coefs[1:p]), 5), "z") 
    out <- list(coefficients = coef.tab, Var.betas = Var.betas)
    out$logLik <- object$log.Lik
    df <- if (!is.null(constr <- object$constraint)) length(coefs) - nrow(constr) else length(coefs)
    out$AIC <- -2 * object$log.Lik + 2 * df
    out$BIC <- -2 * object$log.Lik + log(nrow(object$X)) * df
    out$max.sc <- object$max.sc
    out$conv <- object$conv
    out$counts <- object$counts
    out$call <- object$call
    out$control <- object$control
    out$IRT.param <- object$IRT.param
    out$attr <- attr(object$X, "items")
    out$ancr <- attr(object$X, "anchoring")
    class(out) <- "summ.rasch"
    out
}
