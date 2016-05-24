residuals.rasch <-
function (object, resp.patterns = NULL, order = TRUE, ...) {
    if (!inherits(object, "rasch"))
        stop("Use only with 'rasch' objects.\n")
    if (any(is.na(object$X)))
        warning("residuals are not meaningful for patterns with missing data.\n")
    fits <- fitted(object, resp.patterns = resp.patterns)
    X <- fits[, -ncol(fits), drop = FALSE]
    Exp <- fits[, "Exp"]
    Obs <- observedFreqs(object, X)
    out <- cbind(X, Obs = Obs, Exp = Exp, Resid = (Obs - Exp) / sqrt(Exp))
    if (order)
        out <- out[order(out[, "Resid"]), ]
    out
}
