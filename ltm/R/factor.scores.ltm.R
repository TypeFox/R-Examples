factor.scores.ltm <-
function (object, resp.patterns = NULL, method = c("EB", "EAP", "MI", "Component"), B = 5, 
        robust.se = FALSE, prior = TRUE, return.MIvalues = FALSE, ...) {
    if (!inherits(object, "ltm"))
        stop("Use only with 'ltm' objects.\n")
    betas <- object$coef
    fits <- fitted(object, resp.patterns = resp.patterns)
    X <- fits[, -ncol(fits), drop = FALSE]
    nx <- nrow(X)
    factors <- object$ltst$factors
    inter <- object$ltst$inter
    quad.z1 <- object$ltst$quad.z1
    quad.z2 <- object$ltst$quad.z2
    p <- nrow(betas)
    q. <- 1 + factors + sum(inter, quad.z1, quad.z2)
    form <- as.character(object$formula)
    form <- as.formula(paste("~", form[3]))
    method <- match.arg(method)
    if (any(inter, quad.z1, quad.z2) && method == "Component") {
        warning("In presence of nonlinear terms, Component Scores give biased factor-scores estimates. The MI method is used instead.\n")
        method <- "MI"
    }
    Obs <- observedFreqs(object, X)
    res <- data.frame(X, Obs = Obs, Exp = fits[, ncol(fits)])
    names(res)[1:p] <- rownames(betas)
    if (method == "Component") {
        res$z1 <- colSums(t(X) * betas[, 2])
        if (factors == 2)
            res$z2 <- colSums(t(X) * betas[, 3])
    }
    if (method %in% c("EB", "EAP", "MI")) {
        environment(fscores.l) <- environment()
        res <- fscores.l(betas, X, method)
    }
    out <- list(score.dat = res, method = method, B = B, call = object$call, resp.pats = !is.null(resp.patterns),
                coef = coef(object))
    class(out) <- "fscores"
    out
}
