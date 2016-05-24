gpcm <-
function (data, constraint = c("gpcm", "1PL", "rasch"), IRT.param = TRUE, start.val = NULL, na.action = NULL, control = list()) {
    cl <- match.call()
    if ((!is.data.frame(data) & !is.matrix(data)) || ncol(data) == 1)
        stop("'data' must be either a numeric matrix or a data.frame, with at least two columns.\n")
    constraint <- match.arg(constraint)
    X <- if (!is.data.frame(data)) as.data.frame(data) else data
    X[] <- lapply(X, factor)
    ncatg <- as.vector(sapply(X, function (x) length(levels(x))))
    X <- sapply(X, unclass)
    if (!is.null(na.action))
        X <- na.action(X)
    colnamsX <- colnames(X)
    dimnames(X) <- NULL
    p <- ncol(X)
    pats <- apply(X, 1, paste, collapse = "/")
    freqs <- table(pats)
    nfreqs <- length(freqs)
    obs <- as.vector(freqs)
    X <- unlist(strsplit(cbind(names(freqs)), "/"))
    X[X == "NA"] <- as.character(NA)
    X <- matrix(as.numeric(X), nfreqs, p, TRUE)
    XX <- lapply(1:p, function (j) outer(X[, j], seq(1, ncatg[j] - 1), ">") * 1)
    con <- list(iter.qN = 150, GHk = 21, optimizer = "nlminb", optimMethod = "BFGS", numrDeriv = "fd",
        epsHes = 1e-06, parscale = NULL, verbose = getOption("verbose"))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% namC]) > 0) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    GH <- GHpoints(data ~ z1, con$GHk)
    Z <- GH$x[, 2]
    GHw <- GH$w
    init.thetas <- start.val.gpcm(start.val, X, obs, constraint, ncatg, IRT.param)
    environment(loglikgpcm) <- environment(scoregpcm) <- environment()
    res.qN <- if (con$optimizer == "optim") {
        if (is.null(con$parscale) || length(con$parscale) != length(init.thetas))
            con$parscale <- rep(0.5, length(init.thetas))
        optim(init.thetas, loglikgpcm, scoregpcm, constraint = constraint, method = con$optimMethod,
              control = list(maxit = con$iter.qN, parscale = con$parscale, trace = as.numeric(con$verbose)))
    } else {
        nlb <- nlminb(init.thetas, objective = loglikgpcm, gradient = scoregpcm, constraint = constraint,
                      control = list(iter.max = con$iter.qN, eval.max = con$iter.qN, trace = as.numeric(con$verbose),
                      rel.tol = sqrt(.Machine$double.eps)))
        names(nlb) <- c("par", "value", "convergence", "message", "iterations", "counts")
        nlb
    }
    Hess <- if (con$numrDeriv == "fd")
        fd.vec(res.qN$par, scoregpcm, constraint = constraint, eps = con$epsHes)
    else
        cd.vec(res.qN$par, scoregpcm, constraint = constraint, eps = con$epsHes)
    res.qN$hessian <- 0.5 * (Hess + t(Hess))
    if (res.qN$convergence != 0) {
        if (!is.null(res.qN$message))
            warning("Not successful convergence: ", res.qN$message, ".\n")
        else
            warning("Not successful convergence.\n")
    }
    if (all(!is.na(res.qN$hessian) & is.finite(res.qN$hessian))) {
        ev <- eigen(res.qN$hessian, TRUE, TRUE)$values
        if (!all(ev >= -1e-06 * abs(ev[1]))) 
            warning("Hessian matrix at convergence is not positive definite; unstable solution.\n")
    } else 
        warning("Hessian matrix at convergence contains infinite or missing values; unstable solution.\n")
    thetas <- betas.gpcm(res.qN$par, p, ncatg, constraint)
    names(thetas) <- if (!is.null(colnamsX)) colnamsX else paste("Item", 1:p)
    thetas <- lapply(thetas, function (x) { names(x) <- c(paste("Catgr.", seq(1, length(x) - 1), sep = ""), "Dscrmn"); x })
    max.sc <- max(abs(scoregpcm(res.qN$par, constraint)), na.rm = TRUE)
    fit <- list(coefficients = thetas, log.Lik = -res.qN$value, convergence = res.qN$conv, hessian = res.qN$hessian, 
                counts = res.qN$counts, patterns = list(X = X, obs = obs), GH = list(Z = Z, GHw = GHw), max.sc = max.sc, 
                constraint = constraint, IRT.param = IRT.param, X = data, control = con, na.action = na.action,
                call = cl)
    class(fit) <- "gpcm"
    fit
}
