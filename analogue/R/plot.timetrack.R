`plot.timetrack` <- function(x, choices = 1:2,
                             display = c("wa","lc"),
                             order,
                             type = c("p", "n"),
                             ptype = c("l", "p", "o", "b", "n"),
                             pch = c(1,2),
                             col = c("black","red"),
                             lty = "solid", lwd = 1,
                             xlim = NULL, ylim = NULL, ...) {
    ptype <- match.arg(ptype)
    type <- match.arg(type)
    display <- match.arg(display)
    scrs <- scores(x$ord, choices = choices, scaling = x$scaling,
                   display = display, ...)
    pass <- fitted(x, type = "passive", choices = choices)
    if (is.null(xlim)) {
        xlim <- range(scrs[,1], pass[,1])
    }
    if (is.null(ylim)) {
        ylim <- range(scrs[,2], pass[,2])
    }
    plt <- plot(x$ord, choices = choices, scaling = x$scaling,
                type = "n", display = display, ...,
                ylim = ylim, xlim = xlim)
    if (isTRUE(all.equal(type, "p"))) {
        points(scrs, pch = pch[1], col = col[1], ...)
    }
    if(!missing(order)) {
        if(length(order) != NROW(pass))
            stop("'length(order)' not equal to number of passive samples.")
        pass[order, ]
    }
    if(ptype %in% c("l", "o", "b")) {
        lines(pass, pch = pch[2], col = col[2],
              lty = lty, lwd = lwd, type = ptype, ...)
    } else if (isTRUE(all.equal(ptype, "p"))) {
        points(pass, pch = pch[2], col = col[2], ...)
    }
    invisible(x$ord)
}
