blup.rma.uni <-
function (x, level, digits, transf, targs, ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (missing(level)) 
        level <- x$level
    if (missing(digits)) 
        digits <- x$digits
    if (missing(transf)) 
        transf <- FALSE
    if (missing(targs)) 
        targs <- NULL
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    if (x$knha) {
        crit <- qt(alpha/2, df = x$dfs, lower.tail = FALSE)
    }
    else {
        crit <- qnorm(alpha/2, lower.tail = FALSE)
    }
    pred <- rep(NA_real_, x$k.f)
    vpred <- rep(NA_real_, x$k.f)
    li <- x$tau2/(x$tau2 + x$vi.f)
    for (i in seq_len(x$k.f)[x$not.na]) {
        Xi <- matrix(x$X.f[i, ], nrow = 1)
        pred[i] <- li[i] * x$yi.f[i] + (1 - li[i]) * Xi %*% x$b
        vpred[i] <- li[i] * x$vi.f[i] + (1 - li[i])^2 * Xi %*% 
            tcrossprod(x$vb, Xi)
    }
    se <- sqrt(vpred)
    pi.lb <- pred - crit * se
    pi.ub <- pred + crit * se
    if (is.function(transf)) {
        if (is.null(targs)) {
            pred <- sapply(pred, transf)
            se <- rep(NA, x$k.f)
            pi.lb <- sapply(pi.lb, transf)
            pi.ub <- sapply(pi.ub, transf)
        }
        else {
            pred <- sapply(pred, transf, targs)
            se <- rep(NA, x$k.f)
            pi.lb <- sapply(pi.lb, transf, targs)
            pi.ub <- sapply(pi.ub, transf, targs)
        }
        transf <- TRUE
    }
    pi.bounds <- cbind(pi.lb, pi.ub)
    rev.order <- ifelse(pi.ub < pi.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    pi.bounds[rev.order, ] <- pi.bounds[rev.order, 2:1]
    pi.lb <- pi.bounds[, 1]
    pi.ub <- pi.bounds[, 2]
    if (na.act == "na.omit") {
        out <- list(pred = pred[x$not.na], se = se[x$not.na], 
            pi.lb = pi.lb[x$not.na], pi.ub = pi.ub[x$not.na])
        out$slab <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(pred = pred, se = se, pi.lb = pi.lb, pi.ub = pi.ub)
        out$slab <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    out$digits <- digits
    out$transf <- transf
    class(out) <- c("list.rma")
    return(out)
}
