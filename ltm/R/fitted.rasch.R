fitted.rasch <-
function (object, resp.patterns = NULL, 
    type = c("expected", "marginal-probabilities", "conditional-probabilities"), ...) {
    if (!inherits(object, "rasch"))
        stop("Use only with 'rasch' objects.\n")
    type <- match.arg(type)
    X <- if (is.null(resp.patterns)) {
        data.matrix(object$patterns$X)
    } else {
        if (!is.matrix(resp.patterns) && !is.data.frame(resp.patterns))
            stop("'resp.patterns' should be a matrix or a data.frame.\n")
        resp.patterns <- data.matrix(resp.patterns)
        resp.patterns <- apply(resp.patterns, 2, function (x) if (all(unique(x) %in% c(1, 0, NA))) x else x - 1)
        if (!is.matrix(resp.patterns))
            resp.patterns <- t(resp.patterns)
        p <- ncol(object$X)
        if (ncol(resp.patterns) != p)
            stop("the number of items in ", deparse(substitute(object)), " and the number of columns of 'resp.patterns' do not much.\n")
        check.items <- vector("logical", p)
        for (i in 1:p)
            check.items[i] <- all(unique(resp.patterns[, i]) %in% c(unique(object$patterns$X[, i]), NA))
        if (!all(check.items)) {
            its <- paste((1:p)[!check.items], collapse = ", ")
            stop("the number of levels in 'resp.patterns' does not much for item(s): ", its, "\n")
        }
        resp.patterns
    }
    colnames(X) <- colnames(object$X)
    mX <- 1 - X
    if (any(na.ind <- is.na(X)))
        X[na.ind] <- mX[na.ind] <- 0
    pr <- if (type == "expected" || type == "marginal-probabilities") {
        probs(object$GH$Z %*% t(object$coef))
    } else {
        Z <- cbind(1, factor.scores(object, resp.patterns = resp.patterns)$score.dat$z1)
        probs(Z %*% t(object$coef))
    }
    p.xz <- exp(X %*% t(log(pr)) + mX %*% t(log(1 - pr)))
    X[na.ind] <- NA
    out <- switch(type,
        "expected" = cbind(X, Exp = nrow(object$X) * colSums(object$GH$GHw * t(p.xz))),
        "marginal-probabilities" = cbind(X, "Marg-Probs" = colSums(object$GH$GHw * t(p.xz))),
        "conditional-probabilities" = pr)
    rownames(out) <- if (!is.null(resp.patterns) && !is.null(nams <- rownames(resp.patterns))) nams else NULL
    out
}
