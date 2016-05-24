factor.scores.rasch <-
function (object, resp.patterns = NULL, method = c("EB", "EAP", "MI"), B = 5, 
        robust.se = FALSE, prior = TRUE, return.MIvalues = FALSE, ...) {
    if (!inherits(object, "rasch"))
        stop("Use only with 'rasch' objects.\n")
    betas <- object$coef
    fits <- fitted(object, resp.patterns = resp.patterns)
    X <- fits[, -ncol(fits), drop = FALSE]
    nx <- nrow(X)
    p <- nrow(betas)
    method <- match.arg(method)
    Obs <- observedFreqs(object, X)
    res <- data.frame(X, Obs = Obs, Exp = fits[, ncol(fits)])        
    names(res)[1:p] <- rownames(betas)
    environment(fscores.r) <- environment()
    res <- fscores.r(betas, X, method)
    out <- list(score.dat = res, method = method, B = B, call = object$call, resp.pats = !is.null(resp.patterns),
                coef = coef(object))
    class(out) <- "fscores"
    out
}
