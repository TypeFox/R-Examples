permutest.rma.uni <-
function (x, exact = FALSE, iter = 1000, progbar = TRUE, retpermdist = FALSE, 
    digits, tol, ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    if (missing(digits)) 
        digits <- x$digits
    if (missing(tol)) 
        tol <- .Machine$double.eps^0.5
    if (x$int.only) {
        exact.iter <- 2^x$k
    }
    else {
        X <- as.data.frame(x$X)[do.call(order, as.data.frame(x$X)), 
            ]
        indices <- cumsum(c(TRUE, !duplicated(X)[-1]))
        indices <- rep(cumsum(rle(indices)$length) - (rle(indices)$length - 
            1), rle(indices)$length)
        ind.table <- table(indices)
        exact.iter <- round(prod((max(ind.table) + 1):x$k)/prod(factorial(ind.table[-which.max(ind.table)])))
    }
    if (exact || (exact.iter <= iter)) {
        exact <- TRUE
        iter <- exact.iter
    }
    if (iter == Inf) 
        stop("Too many iterations required for exact permutation test.\n")
    if (progbar) 
        cat("Running ", iter, " iterations for ", ifelse(exact, 
            "exact", "approximate"), " permutation test.\n", 
            sep = "")
    if (x$int.only) {
        zval.perm <- try(rep(NA_real_, iter), silent = TRUE)
        if (inherits(zval.perm, "try-error")) 
            stop("Number of iterations requested too large.")
        QM.perm <- try(rep(NA_real_, iter), silent = TRUE)
        if (inherits(QM.perm, "try-error")) 
            stop("Number of iterations requested too large.")
        if (progbar) 
            pbar <- txtProgressBar(min = 0, max = iter, style = 3)
        if (exact) {
            signmat <- as.matrix(expand.grid(replicate(x$k, list(c(1, 
                -1))), KEEP.OUT.ATTRS = FALSE))
            for (i in seq_len(iter)) {
                res <- try(rma(signmat[i, ] * x$yi, x$vi, weights = x$weights, 
                  method = x$method, weighted = x$weighted, intercept = TRUE, 
                  knha = x$knha, control = x$control, btt = 1), 
                  silent = FALSE)
                if (inherits(res, "try-error")) 
                  next
                zval.perm[i] <- res$zval
                QM.perm[i] <- res$QM
                if (progbar) 
                  setTxtProgressBar(pbar, i)
            }
        }
        else {
            i <- 1
            while (i <= iter) {
                signs <- sample(c(-1, 1), x$k, replace = TRUE)
                res <- try(rma(signs * x$yi, x$vi, weights = x$weights, 
                  method = x$method, weighted = x$weighted, intercept = TRUE, 
                  knha = x$knha, control = x$control, btt = 1), 
                  silent = FALSE)
                if (inherits(res, "try-error")) 
                  next
                zval.perm[i] <- res$zval
                QM.perm[i] <- res$QM
                i <- i + 1
                if (progbar) 
                  setTxtProgressBar(pbar, i)
            }
        }
        if (!exact) {
            zval.perm[1] <- x$zval
            QM.perm[1] <- x$QM
        }
        if (x$zval > 0) {
            pval <- 2 * mean(zval.perm >= x$zval - tol, na.rm = TRUE)
        }
        else {
            pval <- 2 * mean(zval.perm <= x$zval + tol, na.rm = TRUE)
        }
        pval[pval > 1] <- 1
        QMp <- mean(QM.perm >= x$QM - tol, na.rm = TRUE)
    }
    else {
        zval.perm <- suppressWarnings(try(matrix(NA_real_, nrow = iter, 
            ncol = x$p), silent = TRUE))
        if (inherits(zval.perm, "try-error")) 
            stop("Number of iterations requested too large.")
        QM.perm <- try(rep(NA_real_, iter), silent = TRUE)
        if (inherits(QM.perm, "try-error")) 
            stop("Number of iterations requested too large.")
        if (progbar) 
            pbar <- txtProgressBar(min = 0, max = iter, style = 3)
        if (exact) {
            permmat <- .genuperms(indices)
            for (i in seq_len(iter)) {
                res <- try(suppressWarnings(rma(x$yi, x$vi, weights = x$weights, 
                  mods = cbind(X[permmat[i, ], ]), method = x$method, 
                  weighted = x$weighted, intercept = FALSE, knha = x$knha, 
                  control = x$control, btt = x$btt)), silent = FALSE)
                if (inherits(res, "try-error")) 
                  next
                zval.perm[i, ] <- res$zval
                QM.perm[i] <- res$QM
                if (progbar) 
                  setTxtProgressBar(pbar, i)
            }
        }
        else {
            i <- 1
            while (i <= iter) {
                res <- try(suppressWarnings(rma(x$yi, x$vi, weights = x$weights, 
                  mods = cbind(X[sample(x$k), ]), method = x$method, 
                  weighted = x$weighted, intercept = FALSE, knha = x$knha, 
                  control = x$control, btt = x$btt)), silent = FALSE)
                if (inherits(res, "try-error")) 
                  next
                zval.perm[i, ] <- res$zval
                QM.perm[i] <- res$QM
                i <- i + 1
                if (progbar) 
                  setTxtProgressBar(pbar, i)
            }
        }
        if (!exact) {
            zval.perm[1, ] <- x$zval
            QM.perm[1] <- x$QM
        }
        pval <- rep(NA_real_, x$p)
        for (j in seq_len(x$p)) {
            if (x$zval[j] > 0) {
                pval[j] <- 2 * mean(zval.perm[, j] >= x$zval[j] - 
                  tol, na.rm = TRUE)
            }
            else {
                pval[j] <- 2 * mean(zval.perm[, j] <= x$zval[j] + 
                  tol, na.rm = TRUE)
            }
        }
        pval[pval > 1] <- 1
        QMp <- mean(QM.perm >= x$QM - tol, na.rm = TRUE)
    }
    if (progbar) 
        close(pbar)
    out <- list(pval = pval, QMp = QMp, b = x$b, se = x$se, zval = x$zval, 
        ci.lb = x$ci.lb, ci.ub = x$ci.ub, QM = x$QM, k = x$k, 
        p = x$p, btt = x$btt, m = x$m, knha = x$knha, dfs = x$dfs, 
        int.only = x$int.only, digits = digits, exact.iter = exact.iter)
    if (retpermdist) {
        out$QM.perm <- QM.perm
        out$zval.perm <- data.frame(zval.perm)
        names(out$zval.perm) <- colnames(x$X)
    }
    class(out) <- "permutest.rma.uni"
    return(out)
}
