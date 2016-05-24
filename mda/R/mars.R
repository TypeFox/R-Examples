mars <-
function (x, y, w = rep(1, nrow(x)), wp, degree = 1, nk = max(21, 
    2 * ncol(x) + 1), penalty = 2, thresh = 0.001, prune = TRUE, 
    trace.mars = FALSE, forward.step = TRUE, prevfit = NULL, ...) 
{
    this.call <- match.call()
    if ((nk%%2) != 1) 
        nk <- nk - 1
    x <- as.matrix(x)
    np <- dim(x)
    n <- np[1]
    p <- np[2]
    y <- as.matrix(y)
    nclass <- ncol(y)
    if (is.null(np)) {
        np <- c(length(x), 1)
        x <- as.matrix(x)
    }
    if (forward.step) {
        interms <- 1
        lenb <- nk
        bx <- matrix(rep(0, nrow(x) * nk), nrow = n)
        res <- matrix(rep(0, nrow(x) * ncol(y)), nrow = n)
        fullin <- rep(0, nk)
        cuts <- NULL
        factor <- NULL
    }
    else {
        bx <- model.matrix.mars(prevfit, x, full = TRUE)
        interms <- ncol(bx)
        lenb <- prevfit$lenb
        o <- prevfit$all.terms
        fullin <- rep(0, ncol(bx))
        fullin[o] <- 1
        res <- prevfit$res
        factor <- prevfit$factor
        cuts <- prevfit$cuts
        if (missing(penalty)) 
            penalty <- prevfit$penalty
        degree <- prevfit$degree
        nk <- lenb
        thresh <- prevfit$thresh
    }
    if (missing(penalty) & (degree > 1)) 
        penalty <- 3
    if (!missing(wp)) {
        if (any(wp <= 0)) 
            stop("wp should all be positive")
        wp <- sqrt(wp/sum(wp))
        y <- y * outer(rep(1, n), wp)
    }
    else wp <- NULL
    tagx <- x
    storage.mode(tagx) <- "integer"
    for (j in 1:p) {
        tagx[, j] <- order(x[, j])
    }
    bestin <- rep(0, nk)
    flag <- matrix(rep(0, nk * p), nrow = nk, ncol = p)
    if (is.null(cuts)) 
        cuts <- matrix(rep(0, nk * p), nrow = nk, ncol = p)
    if (is.null(factor)) {
        dir <- matrix(rep(0, nk * p), nrow = nk, ncol = p)
    }
    else {
        dir <- factor
    }
    alpha <- rep(0, nclass)
    beta <- matrix(rep(0, nk * nclass), nrow = nk)
    bestgcv <- 0
    storage.mode(y) <- "double"
    storage.mode(x) <- "double"
    storage.mode(bx) <- "double"
    storage.mode(flag) <- "integer"
    storage.mode(cuts) <- "double"
    storage.mode(dir) <- "double"
    storage.mode(res) <- "double"
    storage.mode(beta) <- "double"
    lenscrat <- 1 + n + 2 * n * nk + 4 * nk * nk + 3 * nk + 3 * 
        nk * nclass + 3 * nclass + 28 * n + 51
    junk <- .Fortran("marss", as.integer(n), as.integer(n), as.integer(p), 
        as.integer(nclass), as.matrix(y), as.matrix(x), as.double(w), 
        as.matrix(tagx), as.integer(degree), as.integer(nk), 
        as.double(penalty), as.double(thresh), as.logical(forward.step), 
        as.integer(interms), as.logical(prune), bx = as.matrix(bx), 
        fullin = as.integer(fullin), lenb = as.integer(lenb), 
        bestgcv = as.double(bestgcv), bestin = as.integer(bestin), 
        flag = as.matrix(flag), cuts = as.matrix(cuts), dir = as.matrix(dir), 
        res = as.matrix(res), alpha = as.double(alpha), beta = as.matrix(beta), 
        double(lenscrat), integer(4 * nk), trace.mars, PACKAGE = "mda")
    lenb <- junk$lenb
    all.terms <- seq(lenb)[junk$fullin[1:lenb] == 1]
    selected.terms <- seq(lenb)[junk$bestin[1:lenb] == 1]
    coefficients <- junk$beta[seq(selected.terms), , drop = FALSE]
    residuals <- junk$res
    fitted.values <- y - residuals
    if (!is.null(wp)) {
        TT <- outer(rep(1, n), wp)
        residuals <- residuals/TT
        fitted.values <- fitted.values/TT
        coefficients <- coefficients/outer(rep(1, length(selected.terms)), 
            wp)
    }
    dir <- junk$dir[seq(lenb), , drop = FALSE]
    dimnames(dir) <- list(NULL, dimnames(x)[[2]])
    cutss <- junk$cuts[seq(lenb), , drop = FALSE]
    x <- junk$bx[, selected.terms, drop = FALSE]
    structure(list(call = this.call, all.terms = all.terms, selected.terms = selected.terms, 
        penalty = penalty, degree = degree, nk = nk, thresh = thresh, 
        gcv = junk$bestgcv, factor = dir, cuts = cutss, residuals = residuals, 
        fitted.values = fitted.values, lenb = junk$lenb, coefficients = coefficients, 
        x = x), class = "mars")
}

