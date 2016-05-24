cumul.rma.uni <-
function (x, order, digits, transf, targs, ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (!x$int.only) 
        stop("Method only applicable for models without moderators.")
    if (missing(order)) 
        order <- NULL
    if (missing(digits)) 
        digits <- x$digits
    if (missing(transf)) 
        transf <- FALSE
    if (missing(targs)) 
        targs <- NULL
    if (is.null(order)) 
        order <- seq_len(x$k.f)
    yi.f <- x$yi.f[order]
    vi.f <- x$vi.f[order]
    X.f <- cbind(x$X.f[order, ])
    weights.f <- x$weights.f[order]
    not.na <- x$not.na[order]
    slab <- x$slab[order]
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
    for (i in seq_len(x$k.f)[not.na]) {
        res <- try(suppressWarnings(rma(yi.f[seq_len(i)], vi.f[seq_len(i)], 
            weights = weights.f[seq_len(i)], method = x$method, 
            weighted = x$weighted, intercept = TRUE, knha = x$knha, 
            control = x$control)), silent = TRUE)
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
    alpha <- ifelse(x$level > 1, (100 - x$level)/100, 1 - x$level)
    crit <- qnorm(alpha/2, lower.tail = FALSE)
    b[1] <- yi.f[1]
    se[1] <- sqrt(vi.f[1])
    zval[1] <- yi.f[1]/se[1]
    pval[1] <- 2 * pnorm(abs(zval[1]), lower.tail = FALSE)
    ci.lb[1] <- yi.f[1] - crit * se[1]
    ci.ub[1] <- yi.f[1] + crit * se[1]
    QE[1] <- 0
    QEp[1] <- 1
    tau2[1] <- 0
    I2[1] <- 0
    H2[1] <- 1
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
        out <- list(estimate = b[not.na], se = se[not.na], zval = zval[not.na], 
            pvals = pval[not.na], ci.lb = ci.lb[not.na], ci.ub = ci.ub[not.na], 
            QE = QE[not.na], QEp = QEp[not.na], tau2 = tau2[not.na], 
            I2 = I2[not.na], H2 = H2[not.na])
        out$slab <- slab[not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(estimate = b, se = se, zval = zval, pvals = pval, 
            ci.lb = ci.lb, ci.ub = ci.ub, QE = QE, QEp = QEp, 
            tau2 = tau2, I2 = I2, H2 = H2)
        out$slab <- slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    if (x$method == "FE") 
        out <- out[-c(9, 10, 11)]
    out$digits <- digits
    out$transf <- transf
    out$slab.null <- x$slab.null
    out$level <- x$level
    out$measure <- x$measure
    out$knha <- x$knha
    if (x$measure == "GEN") {
        attr(out$estimate, "measure") <- "GEN"
    }
    else {
        attr(out$estimate, "measure") <- x$measure
    }
    class(out) <- c("list.rma", "cumul.rma")
    return(out)
}
