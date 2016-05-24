tpm <-
function (data, type = c("latent.trait", "rasch"), constraint = NULL, max.guessing = 1, IRT.param = TRUE, 
                 start.val = NULL, na.action = NULL, control = list()) {
    cl <- match.call()
    if ((!is.data.frame(data) & !is.matrix(data)) || ncol(data) == 1)
        stop("'data' must be either a numeric matrix or a data.frame, with at least two columns.\n")
    X <- data.matrix(data)
    if (any(its <- apply(X, 2, function (x) { x <- x[!is.na(x)]; length(unique(x)) } ) > 2))
        stop("'data' contain more that 2 distinct values for item(s): ", paste(which(its), collapse = ", "))
    X <- apply(X, 2, function (x) if (all(unique(x) %in% c(1, 0, NA))) x else x - 1)
    if (!is.null(na.action))
        X <- na.action(X)    
    oX <- X
    colnamsX <- colnames(X)
    dimnames(X) <- NULL
    p <- ncol(X)
    type <- match.arg(type)
    if (!is.null(constraint)) {
        if (!is.matrix(constraint) || ncol(constraint) != 3)
            stop("'constraint' should be a 3-column matrix rows (read help file).\n")
        if (any(constraint[, 1] < 1 | constraint[, 1] > p))
            stop("the 1st column of 'constraint' denotes the items and it should be a number between 1 and ", p, " (read help file).\n")
        if (any(constraint[, 2] < 1 | constraint[, 2] > 3))
            stop("the 2nd column of 'constraint' should contain numbers between 1 and 3 (read help file).\n")
        constraint[, 1:2] <- round(constraint[, 1:2])
        constraint[constraint[, 2] == 1, 3] <- qlogis(constraint[constraint[, 2] == 1, 3])
    }
    if ((!is.numeric(max.guessing) | length(max.guessing) > 1) || (max.guessing < 0 | max.guessing > 1))
        stop("'max.guessing' must be a scalar between 0 and 1.\n")
    thetas <- start.val.tpm(start.val, oX, type, constraint)
    con <- list(optimizer = "optim", iter.qN = 1000, GHk = 21, method = "BFGS", verbose = getOption("verbose"), 
                eps.hessian = 1e-03, parscale = c(rep(0.5, p - sum(constraint[, 2] == 1)), 
                rep(1, length(thetas) - p + sum(constraint[, 2] == 1))))
    con[names(control)] <- control
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
    n <- nrow(X)
    k <- con$GHk
    GH <- GHpoints(data ~ z1, k)
    Z <- GH$x
    GHw <- GH$w
    environment(logLiktpm) <- environment(scoretpm) <- environment()
    res.qN <- if (con$optimizer == "optim") {
        optim(thetas, fn = logLiktpm, gr = scoretpm, method = con$method,
              control = list(maxit = con$iter.qN, trace = as.numeric(con$verbose), parscale = con$parscale), 
              type = type, constraint = constraint, max.guessing = max.guessing)
    } else {
        out <- nlminb(thetas, objective = logLiktpm, gradient = scoretpm, scale = con$parscale,
                      control = list(iter.max = con$iter.qN, eval.max = con$iter.qN, trace = as.numeric(con$verbose),
                      rel.tol = sqrt(.Machine$double.eps)), type = type, constraint = constraint, 
                      max.guessing = max.guessing)
        names(out) <- c("par", "value", "convergence", "message", "iterations", "counts")
        out
    }
    hess <- cd.vec(res.qN$par, scoretpm, type = type, constraint = constraint, max.guessing = max.guessing, 
                    eps = con$eps.hessian)
    res.qN$hessian <- 0.5 * (hess + t(hess))
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
    thetas <- thetas.tpm(res.qN$par, type, constraint, p)
    if (type != "rasch" && sign(thetas[1, 3]) == -1)
        thetas[, 3] <- - thetas[, 3]
    rownames(thetas) <- if (!is.null(colnamsX)) colnamsX else paste("Item", 1:p)
    colnames(thetas) <- c("c.i", "beta.1i", "beta.2i")
    max.sc <- max(abs(scoretpm(res.qN$par, type, constraint, max.guessing)))
    X[na.ind] <- NA
    fit <- list(coefficients = thetas, log.Lik = -res.qN$value, convergence = res.qN$convergence, hessian = res.qN$hessian, 
                counts = res.qN$counts, patterns = list(X = X, obs = obs), GH = list(Z = Z, GHw = GHw), max.sc = max.sc, 
                type = type, constraint = constraint, max.guessing = max.guessing, IRT.param = IRT.param, X = data, 
                control = con, na.action = na.action, call = cl)
    class(fit) <- "tpm"
    fit
}
