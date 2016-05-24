grm <-
function (data, constrained = FALSE, IRT.param = TRUE, Hessian = FALSE, start.val = NULL, na.action = NULL, 
                 control = list()) {
    cl <- match.call()
    if ((!is.data.frame(data) & !is.matrix(data)) || ncol(data) == 1)
        stop("'data' must be either a numeric matrix or a data.frame, with at least two columns.\n")
    X <- data.matrix(data)
    if (!is.null(na.action))
        X <- na.action(X)
    X <- apply(X, 2, function (x) {
        y <- x[!is.na(x)]
        if (any(y == 0)) x + 1 else x
    })
    colnamsX <- colnames(X)
    dimnames(X) <- NULL
    ncatg <- apply(X, 2, function (x) if (any(is.na(x))) length(unique(x)) - 1 else length(unique(x)))
    n <- nrow(X)
    p <- ncol(X)
    pats <- apply(X, 1, paste, collapse = "/")
    freqs <- table(pats)
    nfreqs <- length(freqs)
    obs <- as.vector(freqs)
    X <- unlist(strsplit(cbind(names(freqs)), "/"))
    X[X == "NA"] <- as.character(NA)
    X <- matrix(as.numeric(X), nfreqs, p, TRUE)
    con <- list(iter.qN = 150, GHk = 21, method = "BFGS", verbose = getOption("verbose"), digits.abbrv = 6)
    con[names(control)] <- control
    GH <- GHpoints(data ~ z1, con$GHk)
    Z <- GH$x[, 2]
    GHw <- GH$w
    ind1 <- if (constrained) c(1, cumsum(ncatg[-p] - 1) + 1) else c(1, cumsum(ncatg[-p]) + 1)
    ind2 <- if (constrained) cumsum(ncatg - 1) else cumsum(ncatg)
    betas <- start.val.grm(start.val, X, obs, constrained, ncatg)
    environment(loglikgrm) <- environment(scoregrm) <- environment()
    old <- options(warn = (-1))
    on.exit(options(old))
    res.qN <- optim(unlist(betas), fn = loglikgrm, gr = scoregrm, method = con$method, hessian = Hessian, 
                    control = list(maxit = con$iter.qN, trace = as.numeric(con$verbose)), constrained = constrained)
    betas <- betas.grm(res.qN$par, constrained, ind1, ind2, p)
    names(betas) <- if (!is.null(colnamsX)) colnamsX else paste("Item", 1:p)
    betas <- lapply(betas, function (x) { names(x) <- c(paste("beta.", seq(1, length(x) - 1), sep = ""), "beta"); x } )
    max.sc <- max(abs(scoregrm(res.qN$par, constrained)), na.rm = TRUE)
    fit <- list(coefficients = betas, log.Lik = -res.qN$value, convergence = res.qN$conv, hessian = res.qN$hessian, 
                counts = res.qN$counts, patterns = list(X = X, obs = obs), GH = list(Z = Z, GHw = GHw), max.sc = max.sc, 
                constrained = constrained, IRT.param = IRT.param, X = data, control = con, na.action = na.action,
                call = cl)
    class(fit) <- "grm"
    fit
}
