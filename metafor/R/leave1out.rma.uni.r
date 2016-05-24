leave1out.rma.uni <-
function (x, digits, transf, targs, ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (!x$int.only) 
        stop("Method only applicable for models without moderators.")
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    if (missing(digits)) 
        digits <- x$digits
    if (missing(transf)) 
        transf <- FALSE
    if (missing(targs)) 
        targs <- NULL
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    b <- rep(NA_real_, x$k.f)
    se <- rep(NA_real_, x$k.f)
    zval <- rep(NA_real_, x$k.f)
    pval <- rep(NA_real_, x$k.f)
    ci.lb <- rep(NA_real_, x$k.f)
    ci.ub <- rep(NA_real_, x$k.f)
    QE <- rep(NA_real_, x$k.f)
    QEp <- rep(NA_real_, x$k.f)
    tau2 <- rep(NA_real_, x$k.f)
    I2 <- rep(NA_real_, x$k.f)
    H2 <- rep(NA_real_, x$k.f)
    for (i in seq_len(x$k.f)[x$not.na]) {
        res <- try(suppressWarnings(rma(x$yi.f[-i], x$vi.f[-i], 
            weights = x$weights.f[-i], method = x$method, weighted = x$weighted, 
            intercept = TRUE, knha = x$knha, control = x$control)), 
            silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        b[i] <- res$b
        se[i] <- res$se
        zval[i] <- res$zval
        pval[i] <- res$pval
        ci.lb[i] <- res$ci.lb
        ci.ub[i] <- res$ci.ub
        QE[i] <- res$QE
        QEp[i] <- res$QEp
        tau2[i] <- res$tau2
        I2[i] <- res$I2
        H2[i] <- res$H2
    }
    if (is.function(transf)) {
        if (is.null(targs)) {
            b <- sapply(b, transf)
            se <- rep(NA, x$k.f)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
        }
        else {
            b <- sapply(b, transf, targs)
            se <- rep(NA, x$k.f)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
        }
        transf <- TRUE
    }
    ci.bounds <- cbind(ci.lb, ci.ub)
    rev.order <- ifelse(ci.ub < ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    ci.bounds[rev.order, ] <- ci.bounds[rev.order, 2:1]
    ci.lb <- ci.bounds[, 1]
    ci.ub <- ci.bounds[, 2]
    if (na.act == "na.omit") {
        out <- list(estimate = b[x$not.na], se = se[x$not.na], 
            zval = zval[x$not.na], pval = pval[x$not.na], ci.lb = ci.lb[x$not.na], 
            ci.ub = ci.ub[x$not.na], Q = QE[x$not.na], Qp = QEp[x$not.na], 
            tau2 = tau2[x$not.na], I2 = I2[x$not.na], H2 = H2[x$not.na])
        out$slab <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(estimate = b, se = se, zval = zval, pval = pval, 
            ci.lb = ci.lb, ci.ub = ci.ub, Q = QE, Qp = QEp, tau2 = tau2, 
            I2 = I2, H2 = H2)
        out$slab <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    if (x$method == "FE") 
        out <- out[-c(9, 10, 11)]
    out$digits <- digits
    out$transf <- transf
    class(out) <- c("list.rma")
    return(out)
}
