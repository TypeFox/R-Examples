summary.flexCPH <-
function (object, ...) {
    out <- list(logLik = object$logLik, AIC = AIC(object), 
        BIC = AIC(object, k = log(length(object$d))), d = object$d, knots = unique(object$knots), 
        call = object$call)
    betas <- object$coefficients$betas
    gammas <- object$coefficients$gammas
    Var <- vcov(object)
    ind <- if (length(gammas) == nrow(Var)) NULL else seq(length(gammas) + 1, nrow(Var))
    if (nind <- length(ind)) {
        sds <- if (nind == 1) sqrt(Var[ind, ind]) else sqrt(diag(Var[ind, ind]))
        out$coefTab <- cbind("value" = betas, "exp(value)" = exp(betas), "std.err" = sds, 
            "z-value" = betas / sds, "p-value" = 2 * pnorm(abs(betas / sds), lower.tail = FALSE))
    }
    class(out) <- "summary.flexCPH"
    out
}
