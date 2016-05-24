summary.weibull.frailty <-
function (object, sand.se = FALSE, ...) {
    out <- list(logLik = object$logLik, AIC = AIC(object), 
        BIC = AIC(object, k = log(length(unique(object$id)))), 
        id = object$id, d = object$y[, 2], call = object$call)
    coefs <- object$coefficients$betas
    nx <- length(coefs)
    ind <- seq_len(nx)
    betas <- coefs[ind]
    Var <- vcov(object)
    if (sand.se)
        Var2 <- vcov(object, sand.se = sand.se)
    if (nind <- length(ind)) {
        sds <- if (nind == 1) sqrt(Var[ind, ind]) else sqrt(diag(Var[ind, ind]))
        if (sand.se) {
            sds2 <- if (nind == 1) sqrt(Var2[ind, ind]) else sqrt(diag(Var2[ind, ind]))
        }
        out$coefTab <- if (sand.se) {
            cbind("value" = betas, "std.err" = sds, "sand s.e." = sds2, 
                "z-value" = betas / sds2, 
                "p-value" = 2 * pnorm(abs(betas / sds2), lower.tail = FALSE))
        } else {
            cbind("value" = betas, "std.err" = sds, "z-value" = betas / sds, 
                  "p-value" = 2 * pnorm(abs(betas / sds), lower.tail = FALSE))            
        }
    }
    out$shape <- object$coefficients$shape
    out$scale <- object$coefficients$scale
    out$var.frailty <- object$coefficients$var.frailty
    class(out) <- "summary.weibull.frailty"
    out
}
