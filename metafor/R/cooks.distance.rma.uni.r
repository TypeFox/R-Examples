cooks.distance.rma.uni <-
function (model, ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    x <- model
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    cook.d <- rep(NA_real_, x$k.f)
    svb <- chol2inv(chol(x$vb))
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
        dfbeta <- x$b - res$b
        cook.d[i] <- (crossprod(dfbeta, svb) %*% dfbeta)
    }
    if (na.act == "na.omit") {
        out <- cook.d[x$not.na]
        names(out) <- x$slab[x$not.na]
    }
    if (na.act == "na.exclude" || na.act == "na.pass") {
        out <- cook.d
        names(out) <- x$slab
    }
    if (na.act == "na.fail" && any(!x$not.na)) 
        stop("Missing values in results.")
    return(out)
}
