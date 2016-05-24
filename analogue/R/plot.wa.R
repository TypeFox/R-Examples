plot.wa <- function(x, which = 1:2,
                    caption = c("Inferred vs Observed",
                    "Residuals vs Fitted"),
                    max.bias = TRUE, n.bias = 10,
                    sub.caption = NULL, main = "",
                    ask = prod(par("mfcol")) <
                    length(which) && dev.interactive(),
                    ...,
                    panel = if (add.smooth) panel.smooth else points,
                    add.smooth = getOption("add.smooth")) {
    if (!is.numeric(which) || any(which < 1) || any(which > 2))
        stop("'which' must be in 1:2")
    show <- rep(FALSE, 2)
    show[which] <- TRUE
    if (any(show[1:2])) {
        Est <- fitted(x)
        Obs <- x$orig.env
        Resi <- resid(x)
    }
    if (is.null(sub.caption)) {
        cal <- x$call
        if (!is.na(m.f <- match("formula", names(cal)))) {
            cal <- cal[c(1, m.f)]
            names(cal)[2] <- ""
        }
        cc <- deparse(cal, 80)
        nc <- nchar(cc[1])
        abbr <- length(cc) > 1 || nc > 75
        sub.caption <- if (abbr)
            paste(substr(cc[1], 1, min(75, nc)), "...")
        else cc[1]
    }
    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if (show[1]) {
        lims <- range(Est, Obs)
        plot(Est, Obs, type = "n", asp = 1, xlim = lims,
             ylim = lims, xlab = "Fitted values", ylab = "Observed", ...)
        abline(0, 1, col = "grey", ...)
        panel(Est, Obs, ...)
        if (one.fig)
            title(sub = sub.caption, ...)
        mtext(caption[1], 3, 0.25)
    }
    if (show[2]) {
        plot(Obs, Resi, type = "n", ylab = "Residuals",
             xlab = "Observed", ...)
        abline(h = 0, col = "grey", ...)
        abline(h = mean(Resi), col = "blue", lty = "dashed")
        if (max.bias) {
            groups <- cut(Obs, breaks = n.bias)
            bias <- aggregate(as.vector(Resi), list(group = groups),
                mean)$x
            interv <- lapply(strsplit(sapply(levels(groups),
                                             function(x) substr(x, 2,
                                                                nchar(x) - 1),
                                             USE.NAMES = FALSE), ","),
                             as.numeric)
            interv <- matrix(unlist(interv), ncol = 2, byrow = TRUE)
            arrows(interv[, 1], bias, interv[, 2], bias,
                   length = ifelse(one.fig, 0.05, 0.01),
                   angle = 90, code = 3, col = "blue")
        }
        panel(Obs, Resi, ...)
        if (one.fig)
            title(sub = sub.caption, ...)
        mtext(caption[2], 3, 0.25)
    }
    if (!one.fig && par("oma")[3] >= 1)
        mtext(sub.caption, outer = TRUE, cex = 1.25)
    invisible()
}
