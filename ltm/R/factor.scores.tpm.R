factor.scores.tpm <-
function (object, resp.patterns = NULL, method = c("EB", "EAP", "MI"), B = 5, 
        prior = TRUE, return.MIvalues = FALSE, ...) {
    if (!inherits(object, "tpm"))
        stop("Use only with 'tpm' objects.\n")
    thetas <- object$coef
    fits <- fitted(object, resp.patterns = resp.patterns)
    X <- fits[, -ncol(fits), drop = FALSE]
    nx <- nrow(X)
    p <- ncol(X)
    method <- match.arg(method)
    Obs <- observedFreqs(object, X)
    res <- data.frame(X, Obs = Obs, Exp = fits[, ncol(fits)])
    names(res)[1:p] <- rownames(thetas)
    environment(fscores.t) <- environment()
    res <- fscores.t(thetas, X, method)
    out <- list(score.dat = res, method = method, B = B, call = object$call, resp.pats = !is.null(resp.patterns),
                coef = coef(object))
    class(out) <- "fscores"
    out
}
