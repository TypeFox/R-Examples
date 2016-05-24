summary.tpm <-
function (object, ...) {
    if (!inherits(object, "tpm"))
        stop("Use only with 'tpm' objects.\n")
    if (object$IRT.param)
        irt <- IRT.parm(object, TRUE)
    coefs <- if (object$IRT.param) irt$parms else object$coef
    coefs[, 1] <- plogis(coefs[, 1]) * object$max.guessing
    Var.betas <- vcov(object)
    coefs <- if (object$type == "rasch") c(coefs[, 1:2], coefs[1, 3]) else c(coefs) 
    se <- if (object$IRT.param) irt$se else {
        ses <- rep(NA, length(coefs))
        ses[attr(Var.betas, "drop.ind")] <- sqrt(diag(Var.betas))
        ses
    }
    coef.tab <- cbind(value = coefs, std.err = se, z.vals = coefs / se)
    p <- ncol(object$X)
    nams <- if (object$IRT) {
         c(t(outer(colnames(irt$parms), rownames(object$coef), paste, sep = ".")))
    } else {
        as.vector(t(outer(c("c.", "beta.1", "beta.2"), as.character(1:p), paste, sep = ""))) 
    }
    rownames(coef.tab) <- if (object$type == "rasch") c(nams[seq(1, 2 * p)], "Dscrmn") else nams
    out <- list(coefficients = coef.tab, Var.betas = Var.betas)
    L <- logLik(object)
    out$logLik <- L
    df <- attr(L, "df")
    out$AIC <- AIC(object)
    out$BIC <- AIC(object, k = log(attr(L, "n")))
    out$max.sc <- object$max.sc
    out$conv <- object$conv
    out$counts <- object$counts
    out$call <- object$call
    out$control <- object$control
    out$IRT.param <- object$IRT.param
    out$type <- object$type
    out$attr <- attr(object$X, "items")
    out$ancr <- attr(object$X, "anchoring")
    class(out) <- "summ.tpm"
    out
}
