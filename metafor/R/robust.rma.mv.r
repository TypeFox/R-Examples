robust.rma.mv <-
function (x, cluster, adjust = TRUE, digits, ...) 
{
    if (!is.element("rma.mv", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.mv\".")
    if (missing(cluster)) 
        stop("Need to specify cluster variable.")
    if (missing(digits)) 
        digits <- x$digits
    alpha <- ifelse(x$level > 1, (100 - x$level)/100, 1 - x$level)
    if (!is.null(x$subset)) 
        cluster <- cluster[x$subset]
    cluster <- cluster[x$not.na]
    if (anyNA(cluster)) 
        stop("No missing values allowed in 'cluster' variable.")
    if (length(cluster) != x$k) 
        stop("Length of variable specified via 'cluster' does not match length of data.")
    n <- length(unique(cluster))
    dfs <- n - x$p
    if (dfs <= 0) 
        stop(paste0("Number of clusters (", n, ") must be larger than the number of fixed effects (", 
            x$p, ")."))
    ocl <- order(cluster)
    cluster <- cluster[ocl]
    if (is.null(x$W)) {
        W <- try(chol(x$M[ocl, ocl]), silent = TRUE)
        if (inherits(W, "try-error")) 
            stop("Cannot take Choleski decomposition of marginal var-cov matrix.")
        W <- chol2inv(W)
        bread <- x$vb %*% crossprod(x$X[ocl, ], W)
    }
    else {
        A <- x$W[ocl, ]
        stXAX <- chol2inv(chol(as.matrix(t(x$X[ocl, ]) %*% A %*% 
            x$X[ocl, ])))
        bread <- stXAX %*% crossprod(x$X[ocl, ], A)
    }
    ei <- c(x$yi - x$X %*% x$b)
    ei <- ei[ocl]
    cluster <- factor(cluster, levels = unique(cluster))
    if (x$sparse) {
        meat <- bdiag(lapply(split(ei, cluster), function(e) tcrossprod(e)))
    }
    else {
        meat <- bldiag(lapply(split(ei, cluster), function(e) tcrossprod(e)))
    }
    vb <- bread %*% meat %*% t(bread)
    if (is.logical(adjust) && adjust) 
        vb <- (n/dfs) * vb
    if (is.character(adjust) && adjust == "Stata") 
        vb <- (n/(n - 1) * (x$k - 1)/(x$k - x$p)) * vb
    b <- x$b
    se <- sqrt(diag(vb))
    names(se) <- NULL
    tval <- c(b/se)
    pval <- 2 * pt(abs(tval), df = dfs, lower.tail = FALSE)
    crit <- qt(alpha/2, df = dfs, lower.tail = FALSE)
    ci.lb <- c(b - crit * se)
    ci.ub <- c(b + crit * se)
    QM <- c(t(b)[x$btt] %*% chol2inv(chol(vb[x$btt, x$btt])) %*% 
        b[x$btt])
    QM <- QM/x$m
    QMp <- pf(QM, df1 = x$m, df2 = dfs, lower.tail = FALSE)
    tcl <- table(cluster)
    res <- list(b = b, se = se, tval = tval, pval = pval, ci.lb = ci.lb, 
        ci.ub = ci.ub, vb = vb, k = x$k, k.f = x$k.f, p = x$p, 
        m = x$m, n = n, dfs = dfs, tcl = tcl, QM = QM, QMp = QMp, 
        yi.f = x$yi.f, vi.f = x$vi.f, X = x$X, X.f = x$X.f, method = x$method, 
        int.only = x$int.only, int.incl = x$int.incl, knha = TRUE, 
        btt = x$btt, intercept = x$intercept, digits = digits, 
        level = x$level, withG = x$withG, withH = x$withH, tau2s = x$tau2s, 
        gamma2s = x$gamma2s, mf.g.f = x$mf.g.f, mf.h.f = x$mf.h.f, 
        g.levels.f = x$g.levels.f, h.levels.f = x$h.levels.f, 
        sigma2 = x$sigma2, tau2 = x$tau2, gamma2 = x$gamma2, 
        slab = x$slab, slab.null = x$slab.null, not.na = x$not.na, 
        fit.stats = x$fit.stats, k.eff = x$k.eff, p.eff = x$p.eff, 
        parms = x$parms, measure = x$measure)
    class(res) <- c("robust.rma", "rma", "rma.mv")
    return(res)
}
