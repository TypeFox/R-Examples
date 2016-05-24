fitted.grm <-
function (object, resp.patterns = NULL, 
        type = c("expected", "marginal-probabilities", "conditional-probabilities"), ...) {
    if (!inherits(object, "grm"))
        stop("Use only with 'grm' objects.\n")
    type <- match.arg(type)
    betas <- object$coef
    p <- length(betas)
    X <- if (is.null(resp.patterns)) {
        object$patterns$X
    } else {
        if (!is.matrix(resp.patterns) && !is.data.frame(resp.patterns))
            stop("'resp.patterns' should be a matrix or a data.frame.\n")
        resp.patterns <- data.matrix(resp.patterns)
        resp.patterns <- apply(resp.patterns, 2, function (x) { y <- x[!is.na(x)]; if (any(y == 0)) x + 1 else x })
        if (!is.matrix(resp.patterns))
            resp.patterns <- t(resp.patterns)
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
    if (type == "expected" || type == "marginal-probabilities") {
        colnames(X) <- names(betas)
        X <- data.matrix(X)
        cpr <- cprobs(betas, object$GH$Z)
        cpr <- lapply(cpr, function (x) log(rbind(x[1, ], diff(x))))
        log.p.xz <- matrix(0, nrow(X), object$control$GHk)
        for (j in 1:p) {
            log.pr <- cpr[[j]]
            xj <- X[, j]
            na.ind <- is.na(xj)
            log.pr <- log.pr[xj, , drop = FALSE]
            if (any(na.ind))
                log.pr[na.ind, ] <- 0
            log.p.xz <- log.p.xz + log.pr
        }
        p.xz <- exp(log.p.xz)
        out <- switch(type,
            "expected" = cbind(X, Exp = nrow(object$X) * colSums(object$GH$GHw * t(p.xz))),
            "marginal-probabilities" = cbind(X, "Marg-Probs" = colSums(object$GH$GHw * t(p.xz))))
        rownames(out) <- if (!is.null(resp.patterns) && !is.null(nams <- rownames(resp.patterns))) nams else NULL
        out
    } else {
        Z <- factor.scores(object, resp.patterns = resp.patterns)$score.dat$z1
        res <- vector(mode = "list", length = p)
        out <- lapply(cprobs(betas, Z), function (x) {
            res <- t(x)
            nc <- ncol(res)
            res <- cbind(res[, 1, drop = FALSE], res[, 2:nc, drop = FALSE] - res[, -nc, drop = FALSE])
            rownames(res) <- if (!is.null(resp.patterns) && !is.null(nams <- rownames(resp.patterns))) nams else NULL
            res
        })
        for (i in seq_along(out)) {
            if (is.factor(object$X[[i]]))
                colnames(out[[i]]) <- levels(object$X[[i]])
        }
        out
    }
}
