baujat.rma.uni <-
function (x, xlim, ylim, xlab, ylab, cex, grid = TRUE, ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    delpred <- rep(NA_real_, x$k.f)
    vdelpred <- rep(NA_real_, x$k.f)
    pred.full <- x$X.f %*% x$b
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
        Xi <- matrix(x$X.f[i, ], nrow = 1)
        delpred[i] <- Xi %*% res$b
        vdelpred[i] <- Xi %*% tcrossprod(res$vb, Xi)
    }
    yhati <- (delpred - pred.full)^2/vdelpred
    options(na.action = "na.pass")
    xhati <- 1/(x$tau2 + x$vi.f) * resid(x)^2
    options(na.action = na.act)
    if (missing(cex)) 
        cex <- 0.8
    if (missing(xlab)) {
        if (x$method == "FE") {
            xlab <- ifelse(x$int.only, "Contribution to Overall Heterogeneity", 
                "Contribution to Residual Heterogeneity")
        }
        else {
            xlab <- "Squared Pearson Residual"
        }
    }
    if (missing(ylab)) 
        ylab <- ifelse(x$int.only, "Influence on Overall Result", 
            "Influence on Fitted Value")
    if (missing(xlim)) 
        xlim <- range(xhati, na.rm = TRUE)
    if (missing(ylim)) 
        ylim <- range(yhati, na.rm = TRUE)
    plot(xhati, yhati, pch = 19, col = "white", xlab = xlab, 
        ylab = ylab, cex = cex, xlim = xlim, ylim = ylim, ...)
    if (grid) 
        grid()
    text(xhati, yhati, x$slab, cex = cex, ...)
    invisible(data.frame(x = xhati[x$not.na], y = yhati[x$not.na], 
        slab = x$slab[x$not.na]))
}
