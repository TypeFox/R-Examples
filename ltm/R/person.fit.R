person.fit <-
function (object, alternative = c("less", "greater", "two.sided"), 
    resp.patterns = NULL, FUN = NULL, simulate.p.value = FALSE, B = 1000) {
    if (!class(object) %in% c("ltm", "rasch", "tpm"))
        stop("Use only with 'ltm', 'rasch' or 'tpm' objects.\n")
    if (inherits(object, "ltm") && any(object$ltst$factors > 1, unlist(object$ltst[2:4])))
        stop("currently only the two-parameter logistic model is supported.\n")
    alternative <- match.arg(alternative)
    prsFit <- if (is.null(FUN)) {
        function (X, z, betas) {
            Z <- cbind(1, z)
            pr <- if (inherits(object, "tpm")) {
                cs.mat <- matrix(plogis(betas[, 1]), npats, p, TRUE)
                cs.mat + (1 - cs.mat) * probs(Z %*% t(betas[, 2:3]))
            } else {
                probs(Z %*% t(betas))
            }
            qr <- 1 - pr
            mu <- rowSums(pr * log(pr) + qr * log(qr))
            sigma <- sqrt(rowSums(pr * qr * log(pr / qr)^2))
            l0 <- rowSums(dbinom(X, 1, pr, log = TRUE), na.rm = TRUE)
            cbind(L0 = l0, Lz = (l0 - mu) / sigma)
        }
    } else {
        if (length(formals(FUN)) != 3)
            stop("'FUN' should be a function with 3 arguments; see the help file for more info.\n")
        FUN
    }
    betas <- object$coefficients
    p <- nrow(betas)
    parms <- if (inherits(object, "tpm")) cbind(betas[, 2:3], plogis(betas[, 1])) else betas
    X <- if (is.null(resp.patterns)) {
        object$patterns$X
    } else {
        out <- fitted(object, resp.patterns = resp.patterns)[, 1:p, drop = FALSE]
        if (!is.null(rownames(resp.patterns)))
            rownames(out) <- rownames(resp.patterns)
        out
    }
    npats <- nrow(X)
    fsc <- factor.scores(object, resp.patterns = X)$score.dat
    ablts <- fsc$z1
    se.ablts <- fsc$se.z1
    Tobs <- as.matrix(prsFit(X, ablts, betas))
    if (!simulate.p.value) {
        pvals <- switch(alternative,
            "less" = pnorm(Tobs),
            "greater" = pnorm(Tobs, lower.tail = FALSE),
            "two.sided" = 2 * pnorm(-abs(Tobs)))
        if (is.null(FUN))
            pvals <- pvals[, 2, drop = FALSE]
    } else {
        T.boot <- array(0, dim = c(npats, NCOL(Tobs), B))
        for (b in 1:B) {
            ablts.new <- rnorm(npats, mean = ablts, sd = se.ablts)
            X.new <- rmvlogis(npats, parms, IRT = FALSE, z.vals = ablts.new)
            T.boot[, , b] <- prsFit(X.new, ablts.new, betas)
        }
        pvals <- switch(alternative,
            "less" = apply(T.boot, 3, function (x) x <= Tobs),
            "greater" = apply(T.boot, 3, function (x) x >= Tobs),
            "two.sided" = apply(T.boot, 3, function (x, abs.Tobs) abs(x) >= abs.Tobs, abs.Tobs = abs(Tobs)))
        pvals <- (rowSums(pvals, na.rm = TRUE) + 1) / (B + 1)
        dim(pvals) <- c(npats, ncol(Tobs))
    }
    stat.nam <- if (!is.null(FUN)) "T.stat" else c("L0", "Lz")
    rnams <- if (!is.null(resp.patterns) && !is.null(nams <- rownames(resp.patterns))) nams else 1:npats
    dimnames(Tobs) <- list(rnams, stat.nam)
    dimnames(pvals) <- if (is.null(FUN) && !simulate.p.value) list(rnams, stat.nam[2]) else list(rnams, stat.nam)
    colnames(X) <- if (is.null(cnams <- colnames(object$X))) paste("It", 1:p) else cnams
    out <- list(resp.patterns = X, Tobs = Tobs, p.values = pvals, alternative = alternative,
        FUN = FUN, simulate.p.value = simulate.p.value, B = B, call = object$call)
    class(out) <- "persFit"
    out
}
