rstudent.rma.uni <-
function (model, digits, ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    if (missing(digits)) 
        digits <- x$digits
    tau2.del <- rep(NA_real_, x$k.f)
    delpred <- rep(NA_real_, x$k.f)
    vdelpred <- rep(NA_real_, x$k.f)
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    for (i in seq_len(x$k.f)[x$not.na]) {
        res <- try(suppressWarnings(rma(x$yi.f[-i], x$vi.f[-i], 
            weights = x$weights.f[-i], mods = cbind(x$X.f[-i, 
                ]), method = x$method, weighted = x$weighted, 
            intercept = FALSE, knha = x$knha, control = x$control)), 
            silent = TRUE)
        if (inherits(res, "try-error")) 
            next
        if (any(res$coef.na)) 
            next
        tau2.del[i] <- res$tau2
        Xi <- matrix(x$X.f[i, ], nrow = 1)
        delpred[i] <- Xi %*% res$b
        vdelpred[i] <- Xi %*% tcrossprod(res$vb, Xi)
    }
    delresid <- x$yi.f - delpred
    delresid[abs(delresid) < 100 * .Machine$double.eps] <- 0
    sedelresid <- sqrt(x$vi.f + vdelpred + tau2.del)
    standelres <- delresid/sedelresid
    if (na.act == "na.omit") {
        out <- list(resid = delresid[x$not.na], se = sedelresid[x$not.na], 
            z = standelres[x$not.na])
        out$slab <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- list(resid = delresid, se = sedelresid, z = standelres)
        out$slab <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    out$digits <- digits
    class(out) <- c("list.rma")
    return(out)
}
