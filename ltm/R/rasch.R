rasch <-
function (data, constraint = NULL, IRT.param = TRUE, start.val = NULL, na.action = NULL, 
    control = list(), Hessian = TRUE) {
    cl <- match.call()
    if ((!is.data.frame(data) & !is.matrix(data)) || ncol(data) == 1)
        stop("'data' must be either a numeric matrix or a data.frame, with at least two columns.\n")
    X <- data.matrix(data)
    if (any(its <- apply(X, 2, function (x) { x <- x[!is.na(x)]; length(unique(x)) } ) > 2))
        stop("'data' contain more that 2 distinct values for item(s): ", paste(which(its), collapse = ", "))
    X <- apply(X, 2, function (x) if (all(unique(x) %in% c(1, 0, NA))) x else x - 1)
    if (!is.null(na.action))
        X <- na.action(X)
    colnamsX <- colnames(X)
    dimnames(X) <- NULL
    n <- nrow(X)
    p <- ncol(X)
    con <- list(iter.qN = 150, GHk = 21, method = "BFGS", verbose = getOption("verbose"))
    con[names(control)] <- control
    betas <- start.val.rasch(start.val, X)
    if (!is.null(constraint)) {
        if (!is.matrix(constraint) || (nrow(constraint) > p + 1 | ncol(constraint) != 2))
            stop("'constraint' should be a 2-column matrix with at most ", p + 1, " rows (read help file).\n")
        if (any(constraint[, 1] < 1 | constraint[, 1] > p + 1))
            stop("the 1st column of 'constraint' should between 1 and ", p + 1, " (read help file).\n")
        constraint <- constraint[order(constraint[, 1]), , drop = FALSE]
        constraint[, 1] <- round(constraint[, 1])
        betas[constraint[, 1]] <- NA
    }
    pats <- apply(X, 1, paste, collapse = "/")
    freqs <- table(pats)
    nfreqs <- length(freqs)
    obs <- as.vector(freqs)
    X <- unlist(strsplit(cbind(names(freqs)), "/"))
    X[X == "NA"] <- as.character(NA)
    X <- matrix(as.numeric(X), nfreqs, p, TRUE)
    mX <- 1 - X
    if (any(na.ind <- is.na(X)))
        X[na.ind] <- mX[na.ind] <- 0
    GH <- GHpoints(data ~ z1, con$GHk)
    Z <- GH$x
    GHw <- GH$w
    logLik.rasch <- function (betas, constraint) {
        betas <- betas.rasch(betas, constraint, p)
        pr <- probs(Z %*% t(betas))
        p.xz <- exp(X %*% t(log(pr)) + mX %*% t(log(1 - pr)))
        p.x <- rep(c(p.xz %*% GHw), obs)
        -sum(log(p.x))
    }
    score.rasch <- function (betas, constraint) {
        betas <- betas.rasch(betas, constraint, p)
        pr <- probs(Z %*% t(betas))
        p.xz <- exp(X %*% t(log(pr)) + mX %*% t(log(1 - pr)))
        p.x <- c(p.xz %*% GHw)
        p.zx <- (p.xz / p.x) * obs
        scores <- matrix(0, p, 2)
        for (i in 1:p) {
            ind <- na.ind[, i]
            Y <- outer(X[, i], pr[, i], "-")
            Y[ind, ] <- 0
            scores[i, ] <- -colSums((p.zx * Y) %*% (Z * GHw))
        }    
        if (!is.null(constraint))
            c(scores[, 1], sum(scores[, 2]))[-constraint[, 1]]
        else
            c(scores[, 1], sum(scores[, 2]))
    }
    res.qN <- optim(betas[!is.na(betas)], fn = logLik.rasch, gr = score.rasch, method = con$method, hessian = Hessian, 
                control = list(maxit = con$iter.qN, trace = as.numeric(con$verbose)), constraint = constraint)                
    if (Hessian && all(!is.na(res.qN$hessian) & is.finite(res.qN$hessian))) {
        ev <- eigen(res.qN$hessian, TRUE, TRUE)$values
        if (!all(ev >= -1e-06 * abs(ev[1]))) 
            warning("Hessian matrix at convergence is not positive definite; unstable solution.\n")
    }
    if (Hessian && any(!is.finite(res.qN$hessian))) {
        warning("Hessian matrix at convergence contains infinite or missing values; unstable solution.\n")
    }
    betas <- betas.rasch(res.qN$par, constraint, p)
    rownames(betas) <- if (!is.null(colnamsX)) colnamsX else paste("Item", 1:p)
    colnames(betas) <- c("beta.i", "beta")
    max.sc <- max(abs(score.rasch(res.qN$par, constraint)))
    X[na.ind] <- NA
    fit <- list(coefficients = betas, log.Lik = -res.qN$value, convergence = res.qN$conv, hessian = res.qN$hessian, 
                counts = res.qN$counts, patterns = list(X = X, obs = obs), GH = list(Z = Z, GHw = GHw), max.sc = max.sc, 
                constraint = constraint, IRT.param = IRT.param, X = data, control = con, na.action = na.action, call = cl)
    class(fit) <- "rasch"
    fit
}
