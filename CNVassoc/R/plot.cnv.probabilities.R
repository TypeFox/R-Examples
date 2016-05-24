plot.cnv.probabilities <-
function (x, my.colors = c("black", "red", "blue"), case.control, ylab = "CNV probability", xlab,
    ...)
{
    probs <- attr(x, "probabilities")
    if (missing(case.control)) {
        names(probs) <- attr(x, "num.copies")
        if (missing(xlab))
            xlab <- "copy number"
        mp <- barplot(colMeans(probs), ylab = ylab, xlab = xlab,
            col = my.colors, ...)
        axis(1, mp, attr(x, "num.copies"), tick = FALSE)
    }
    else {
        tt <- unique(case.control)
        if (length(tt) != 2)
          stop("case.control should be a dicotomous variable")
        cc <- length(cl <- unique(case.control))
        tt <- unlist(by(probs, case.control, colMeans))
        tt <- matrix(tt, ncol = cc, )
        if (missing(xlab))
            xlab <- "case/control status"
        mp <- barplot(tt, ylab = ylab, xlab = xlab, col = my.colors,
            ...)
        old.xpd <- par()$pxd
        par(xpd = NA)
        legend(mean(mp), 1.1, legend = attr(x, "num.copies"),
            fill = my.colors, xjust = 0.5, yjust = 0.5, horiz = TRUE,
            title = "copy number")
        axis(1, mp, cl, tick = FALSE)
        par(xpd = old.xpd)
    }
    invisible()
}
