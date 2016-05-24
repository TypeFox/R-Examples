factor.scores.gpcm <-
function (object, resp.patterns = NULL, method = c("EB", "EAP", "MI"), B = 5, 
        robust.se = FALSE, prior = TRUE, return.MIvalues = FALSE, ...) {
    if (!inherits(object, "gpcm"))
        stop("Use only with 'gpcm' objects.\n")
    betas <- object$coefficients
    fits <- fitted(object, resp.patterns = resp.patterns)
    X <- fits[, -ncol(fits), drop = FALSE]
    nx <- nrow(X)
    p <- length(betas)
    method <- match.arg(method)
    if (method == "MI" && is.null(object$hessian)) {
        warning("object does not have an estimate of the Hessian; the 'EB' method is used instead.\n")
        method <- "EB"
    }
    Obs <- observedFreqs(object, X)
    res <- data.frame(X, Obs = Obs, Exp = fits[, ncol(fits)])
    names(res)[1:p] <- names(betas)
    environment(fscores.gp) <- environment()
    res <- fscores.gp(betas, X, method)
    out <- list(score.dat = res, method = method, B = B, call = object$call, resp.pats = !is.null(resp.patterns))
    class(out) <- "fscores"
    out
}
