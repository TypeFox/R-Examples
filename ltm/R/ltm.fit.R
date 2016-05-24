ltm.fit <-
function (X, betas, constraint, formula, control) {
    n <- nrow(X)
    p <- ncol(X)
    colnamsX <- colnames(X)
    dimnames(X) <- NULL
    pats <- apply(X, 1, paste, collapse = "")
    freqs <- table(pats)
    obs <- as.vector(freqs)
    X <- apply(cbind(names(freqs)), 1, function (x){
                nx <- nchar(x)
                out <- substring(x, 1:nx, 1:nx)
                out <- out[out != "A"]
                out[out == "N"] <- NA
                out
        })
    X <- as.numeric(t(X))
    dim(X) <- c(length(freqs), p)
    mX <- 1 - X
    if (any(na.ind <- is.na(X)))
        X[na.ind] <- mX[na.ind] <- 0
    GH <- GHpoints(formula, control$GHk)
    Z <- GH$x
    GHw <- GH$w
    q. <- ncol(Z)
    environment(EM) <- environment(loglikltm) <- environment(scoreltm) <- environment()
    res.EM <- EM(betas, constraint, control$iter.em, control$verbose)
    if (!is.null(constraint))
        res.EM[constraint[, 1:2, drop = FALSE]] <- NA
    psc <- control$parscale 
    nb <- length(c(res.EM[!is.na(res.EM)]))
    psc <- if (!is.null(psc) && length(psc) == nb) psc else rep(1, nb)
    res.qN <- optim(c(res.EM[!is.na(res.EM)]), fn = loglikltm, gr = scoreltm, method = control$method, hessian = TRUE, 
                control = list(maxit = control$iter.qN, trace = as.numeric(control$verbose), parscale = psc), constraint = constraint)
    if (all(!is.na(res.qN$hessian) & is.finite(res.qN$hessian))) {
        ev <- eigen(res.qN$hessian, TRUE, TRUE)$values
        if (!all(ev >= -1e-06 * abs(ev[1]))) 
            warning("Hessian matrix at convergence is not positive definite; unstable solution.\n")
    } else 
        warning("Hessian matrix at convergence contains infinite or missing values; unstable solution.\n")
    betas <- betas.ltm(res.qN$par, constraint, p, q.)
    colnames(betas) <- GH$colnams
    rownames(betas) <- if (!is.null(colnamsX)) colnamsX else paste("Item", 1:p)
    if (q. == 2 && sign(betas[1, 2]) == -1)
        betas[, 2] <- -betas[, 2]
    max.sc <- max(abs(scoreltm(res.qN$par, constraint)))
    X[na.ind] <- NA
    list(coefficients = betas, log.Lik = -res.qN$value, convergence = res.qN$conv, hessian = res.qN$hessian, 
            counts = res.qN$counts, patterns = list(X = X, obs = obs), GH = list(Z = Z, GHw = GHw), max.sc = max.sc)
}
