confint.rma.uni <-
function (object, parm, level, fixed = FALSE, random = TRUE, 
    digits, transf, targs, verbose = FALSE, control, ...) 
{
    if (!is.element("rma.uni", class(object))) 
        stop("Argument 'object' must be an object of class \"rma.uni\".")
    x <- object
    k <- x$k
    p <- x$p
    yi <- x$yi
    vi <- x$vi
    X <- x$X
    Y <- cbind(yi)
    weights <- x$weights
    if (missing(level)) 
        level <- x$level
    if (missing(digits)) 
        digits <- x$digits
    if (missing(transf)) 
        transf <- FALSE
    if (missing(targs)) 
        targs <- NULL
    if (missing(control)) 
        control <- list()
    if (!fixed && !random) 
        stop("At least one of the arguments 'fixed' and 'random' must be TRUE.")
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    if (random) {
        if (k == 1) 
            stop("Stopped because k = 1.")
        if (x$method == "GENQ") {
            type <- "GENQ"
        }
        else {
            type <- "QP"
        }
        if (x$method == "FE") 
            stop("Model does not contain a random-effects component.")
        if (x$tau2.fix) 
            stop("Model does not contain an estimated random-effects component.")
        if (type == "GENQ" && x$method != "GENQ") 
            stop("Model must be fitted with 'method=\"GENQ\" to use this option.")
        tau2.min <- ifelse(is.null(x$control$tau2.min), min(0, 
            x$tau2), x$control$tau2.min)
        tau2.max <- ifelse(is.null(x$control$tau2.max), max(100, 
            abs(x$tau2) * 10, abs(tau2.min) * 10), x$control$tau2.max)
        con <- list(tol = .Machine$double.eps^0.25, maxiter = 1000, 
            tau2.min = tau2.min, tau2.max = tau2.max, verbose = FALSE)
        con[pmatch(names(control), names(con))] <- control
        if (verbose) 
            con$verbose <- verbose
        verbose <- con$verbose
        tau2.lb <- NA
        tau2.ub <- NA
        ci.null <- FALSE
        lb.conv <- FALSE
        ub.conv <- FALSE
        lb.sign <- ""
        ub.sign <- ""
        if (type == "QP") {
            if (!x$allvipos) 
                stop("Cannot compute confidence interval for the amount of (residual)\n  heterogeneity with non-positive sampling variances in the data.")
            crit.u <- qchisq(alpha/2, k - p, lower.tail = FALSE)
            crit.l <- qchisq(alpha/2, k - p, lower.tail = TRUE)
            QE.tau2.max <- .QE.func(con$tau2.max, Y = Y, vi = vi, 
                X = X, k = k, objective = 0, verbose = FALSE)
            QE.tau2.min <- .QE.func(con$tau2.min, Y = Y, vi = vi, 
                X = X, k = k, objective = 0, verbose = FALSE)
            if (QE.tau2.min < crit.l) {
                tau2.lb <- con$tau2.min
                tau2.ub <- con$tau2.min
                lb.sign <- "<"
                ub.sign <- "<"
                lb.conv <- TRUE
                ub.conv <- TRUE
                if (con$tau2.min <= 0) 
                  ci.null <- TRUE
            }
            else {
                if (QE.tau2.max > crit.l) {
                  tau2.ub <- con$tau2.max
                  ub.sign <- ">"
                  ub.conv <- TRUE
                }
                else {
                  res <- try(uniroot(.QE.func, interval = c(con$tau2.min, 
                    con$tau2.max), tol = con$tol, maxiter = con$maxiter, 
                    Y = Y, vi = vi, X = X, k = k, objective = crit.l, 
                    verbose = verbose, digits = digits)$root, 
                    silent = TRUE)
                  if (!inherits(res, "try-error")) {
                    tau2.ub <- res
                    ub.conv <- TRUE
                  }
                }
            }
            if (QE.tau2.max > crit.u) {
                tau2.lb <- con$tau2.max
                tau2.ub <- con$tau2.max
                lb.sign <- ">"
                ub.sign <- ">"
                lb.conv <- TRUE
                ub.conv <- TRUE
            }
            else {
                if (QE.tau2.min < crit.u) {
                  tau2.lb <- con$tau2.min
                  lb.conv <- TRUE
                  if (con$tau2.min > 0) 
                    lb.sign <- "<"
                }
                else {
                  res <- try(uniroot(.QE.func, interval = c(con$tau2.min, 
                    con$tau2.max), tol = con$tol, maxiter = con$maxiter, 
                    Y = Y, vi = vi, X = X, k = k, objective = crit.u, 
                    verbose = verbose, digits = digits)$root, 
                    silent = TRUE)
                  if (!inherits(res, "try-error")) {
                    tau2.lb <- res
                    lb.conv <- TRUE
                  }
                }
            }
        }
        if (type == "GENQ") {
            if (!requireNamespace("CompQuadForm", quietly = TRUE)) 
                stop("Please install the 'CompQuadForm' package when method='QGEN'.")
            A <- diag(weights, nrow = k, ncol = k)
            stXAX <- .invcalc(X = X, W = A, k = k)
            P <- A - A %*% X %*% stXAX %*% t(X) %*% A
            Q <- crossprod(Y, P) %*% Y
            GENQ.tau2.max <- .GENQ.func(con$tau2.max, P = P, 
                vi = vi, Q = Q, alpha = 0, k = k, p = p, getlower = TRUE, 
                verbose = FALSE)
            GENQ.tau2.min <- .GENQ.func(con$tau2.min, P = P, 
                vi = vi, Q = Q, alpha = 0, k = k, p = p, getlower = TRUE, 
                verbose = FALSE)
            if (GENQ.tau2.min > 1 - alpha/2) {
                tau2.lb <- con$tau2.min
                tau2.ub <- con$tau2.min
                lb.sign <- "<"
                ub.sign <- "<"
                lb.conv <- TRUE
                ub.conv <- TRUE
                if (con$tau2.min <= 0) 
                  ci.null <- TRUE
            }
            else {
                if (GENQ.tau2.max < 1 - alpha/2) {
                  tau2.ub <- con$tau2.max
                  ub.sign <- ">"
                  ub.conv <- TRUE
                }
                else {
                  res <- try(uniroot(.GENQ.func, c(con$tau2.min, 
                    con$tau2.max), P = P, vi = vi, Q = Q, alpha = alpha/2, 
                    k = k, p = p, getlower = FALSE, verbose = verbose, 
                    digits = digits)$root, silent = TRUE)
                  if (!inherits(res, "try-error")) {
                    tau2.ub <- res
                    ub.conv <- TRUE
                  }
                }
            }
            if (GENQ.tau2.max < alpha/2) {
                tau2.lb <- con$tau2.max
                tau2.ub <- con$tau2.max
                lb.sign <- ">"
                ub.sign <- ">"
                lb.conv <- TRUE
                ub.conv <- TRUE
            }
            else {
                if (GENQ.tau2.min > alpha/2) {
                  tau2.lb <- con$tau2.min
                  lb.conv <- TRUE
                  if (con$tau2.min > 0) 
                    lb.sign <- "<"
                }
                else {
                  res <- try(uniroot(.GENQ.func, c(con$tau2.min, 
                    con$tau2.max), P = P, vi = vi, Q = Q, alpha = alpha/2, 
                    k = k, p = p, getlower = TRUE, verbose = verbose, 
                    digits = digits)$root, silent = TRUE)
                  if (!inherits(res, "try-error")) {
                    tau2.lb <- res
                    lb.conv <- TRUE
                  }
                }
            }
        }
        if (type == "PL") {
            if (con$tau2.min > x$tau2) 
                stop("Lower bound of interval to be searched must be <= actual value of component.")
            if (con$tau2.max < x$tau2) 
                stop("Upper bound of interval to be searched must be >= actual value of component.")
            objective <- qchisq(1 - alpha, df = 1)
            res <- try(.profile.rma.uni(val = con$tau2.min, obj = x, 
                CI = TRUE, objective = objective, verbose = verbose), 
                silent = TRUE)
            if (!inherits(res, "try-error")) {
                if (res < 0) {
                  tau2.lb <- con$tau2.min
                  lb.conv <- TRUE
                  if (con$tau2.min > 0) 
                    lb.sign <- "<"
                }
                else {
                  res <- try(uniroot(.profile.rma.uni, interval = c(con$tau2.min, 
                    x$tau2), tol = con$tol, maxiter = con$maxiter, 
                    obj = x, CI = TRUE, objective = objective, 
                    verbose = verbose, check.conv = TRUE)$root, 
                    silent = TRUE)
                  if (!inherits(res, "try-error")) {
                    tau2.lb <- res
                    lb.conv <- TRUE
                  }
                }
            }
            res <- try(.profile.rma.uni(val = con$tau2.max, obj = x, 
                CI = TRUE, objective = objective, verbose = verbose), 
                silent = TRUE)
            if (!inherits(res, "try-error")) {
                if (res < 0) {
                  tau2.ub <- con$tau2.max
                  ub.conv <- TRUE
                  ub.sign <- ">"
                }
                else {
                  res <- try(uniroot(.profile.rma.uni, interval = c(x$tau2, 
                    con$tau2.max), tol = con$tol, maxiter = con$maxiter, 
                    obj = x, CI = TRUE, objective = objective, 
                    verbose = verbose, check.conv = TRUE)$root, 
                    silent = TRUE)
                  if (!inherits(res, "try-error")) {
                    tau2.ub <- res
                    ub.conv <- TRUE
                  }
                }
            }
        }
        if (!lb.conv) 
            warning("Error in iterative search for the lower bound.")
        if (!ub.conv) 
            warning("Error in iterative search for the upper bound.")
        wi <- 1/vi
        W <- diag(wi, nrow = k, ncol = k)
        stXWX <- .invcalc(X = X, W = W, k = k)
        P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        vi.avg <- (k - p)/.tr(P)
        I2.lb <- 100 * tau2.lb/(vi.avg + tau2.lb)
        I2.ub <- 100 * tau2.ub/(vi.avg + tau2.ub)
        H2.lb <- tau2.lb/vi.avg + 1
        H2.ub <- tau2.ub/vi.avg + 1
        tau2 <- c(x$tau2, tau2.lb, tau2.ub)
        tau <- sqrt(c(ifelse(x$tau2 >= 0, x$tau2, NA), ifelse(tau2.lb >= 
            0, tau2.lb, NA), ifelse(tau2.ub >= 0, tau2.ub, NA)))
        I2 <- c(x$I2, I2.lb, I2.ub)
        H2 <- c(x$H2, H2.lb, H2.ub)
        res.random <- rbind(`tau^2` = tau2, tau = tau, `I^2(%)` = I2, 
            `H^2` = H2)
        colnames(res.random) <- c("estimate", "ci.lb", "ci.ub")
    }
    if (fixed) {
        if (x$knha) {
            crit <- qt(alpha/2, df = x$dfs, lower.tail = FALSE)
        }
        else {
            crit <- qnorm(alpha/2, lower.tail = FALSE)
        }
        b <- x$b
        ci.lb <- c(x$b - crit * x$se)
        ci.ub <- c(x$b + crit * x$se)
        if (is.function(transf)) {
            if (is.null(targs)) {
                b <- sapply(b, transf)
                ci.lb <- sapply(ci.lb, transf)
                ci.ub <- sapply(ci.ub, transf)
            }
            else {
                b <- sapply(b, transf, targs)
                ci.lb <- sapply(ci.lb, transf, targs)
                ci.ub <- sapply(ci.ub, transf, targs)
            }
        }
        res.fixed <- cbind(estimate = b, ci.lb = ci.lb, ci.ub = ci.ub)
        rownames(res.fixed) <- rownames(x$b)
        colnames(res.fixed) <- c("estimate", "ci.lb", "ci.ub")
    }
    res <- list()
    if (fixed) 
        res$fixed <- res.fixed
    if (random) 
        res$random <- res.random
    res$digits <- digits
    if (random) {
        res$ci.null <- ci.null
        res$lb.sign <- lb.sign
        res$ub.sign <- ub.sign
        res$tau2.min <- con$tau2.min
    }
    class(res) <- c("confint.rma")
    return(res)
}
