fitted.gpcm <-
function (object, resp.patterns = NULL, 
        type = c("expected", "marginal-probabilities", "conditional-probabilities"), ...) {
    if (!inherits(object, "gpcm"))
        stop("Use only with 'gpcm' objects.\n")
    type <- match.arg(type)
    betas <- object$coefficients
    p <- length(betas)
    X <- if (is.null(resp.patterns)) {
        object$patterns$X
    } else {
        if (!is.matrix(resp.patterns) && !is.data.frame(resp.patterns))
            stop("'resp.patterns' should be a matrix or a data.frame.\n")
        if (ncol(resp.patterns) != p)
            stop("the number of items in ", deparse(substitute(object)), " and the number of columns of 'resp.patterns' do not much.\n")
        rp <- resp.patterns
        if (!is.data.frame(rp))
            rp <- as.data.frame(rp)
        for (i in 1:p) {
            if (is.factor(rp[[i]])) {
                if (!all(levels(rp[[i]]) %in% levels(object$X[[i]])))
                stop("the levels in the ", i, "th column of 'resp.patterns' does not much with the levels of the ", 
                    i, " item in the original data set.\n")
            } else {
                rp[[i]] <- factor(rp[[i]], levels = sort(unique(object$patterns$X[, i])))
            }
        }
        rp <- sapply(rp, unclass)
        if (!is.matrix(rp))
            rp <- t(rp)
        rp
    }
    if (type == "expected" || type == "marginal-probabilities") {
        colnames(X) <- names(betas)
        log.crf <- crf.GPCM(betas, object$GH$Z, object$IRT.param, log = TRUE)
        log.p.xz <- matrix(0, nrow(X), object$control$GHk)
        for (j in 1:p) {
            log.pr <- log.crf[[j]]
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
        names(Z) <- if (!is.null(resp.patterns) && !is.null(nams <- rownames(resp.patterns))) nams else NULL
        res <- vector(mode = "list", length = p)
        out <- lapply(crf.GPCM(betas, Z, object$IRT.param), t)
        for (i in seq_along(out)) {
            if (is.factor(object$X[[i]]))
                colnames(out[[i]]) <- levels(object$X[[i]])
        }
        out
    }
}
