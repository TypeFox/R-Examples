pairs2 <- function (x, y, xlabels, ylabels, panel = points, ...,
                    row1attop = TRUE, gap = 1) {
    localAxis <- function(side, x, y, xpd, bg, col = NULL, main,
                          oma, xlab, ylab, ... ) {
        if (side%%2 == 1){
            Axis(x, side = side, xpd = NA, ...)
            mtext(xlab,side=side, line=3)
        } else {
            Axis(y, side = side, xpd = NA, ...)
            mtext(ylab,side=side, line=3)
        }
    }
    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
    localPanel <- function(..., main, oma, font.main, cex.main) panel(...)

    dots <- list(...)
    nmdots <- names(dots)
    if (!is.matrix(x)) {
        x <- as.data.frame(x)
        for (i in seq_along(names(x))) {
            if (is.factor(x[[i]]) || is.logical(x[[i]]))
                x[[i]] <- as.numeric(x[[i]])
            if (!is.numeric(unclass(x[[i]])))
                stop("non-numeric argument to 'pairs'")
        }
    } else if (!is.numeric(x))
        stop("non-numeric argument to 'pairs'")

    if (!is.matrix(y)) {
        y <- as.data.frame(y)
        for (i in seq_along(names(y))) {
            if (is.factor(y[[i]]) || is.logical(y[[i]]))
                y[[i]] <- as.numeric(y[[i]])
            if (!is.numeric(unclass(y[[i]])))
                stop("non-numeric argument to 'pairs'")
        }
    } else if (!is.numeric(y))
        stop("non-numeric argument to 'pairs'")

    panel <- match.fun(panel)

    nc.x <- ncol(x)
    nc.y <- ncol(y)
    has.xlabs <- has.ylabs <- TRUE
    if (missing(xlabels)) {
        xlabels <- colnames(x)
        if (is.null(xlabels))
            xlabels <- paste("xvar", 1:nc.x)
    } else if (is.null(xlabels)) has.xlabs <- FALSE
    if (missing(ylabels)) {
        ylabels <- colnames(y)
        if (is.null(ylabels))
            ylabels <- paste("yvar", 1:nc.x)
    }
    else if (is.null(ylabels))
        has.ylabs <- FALSE
    oma <- if ("oma" %in% nmdots)
        dots$oma
    else NULL
    main <- if ("main" %in% nmdots)
        dots$main
    else NULL
    if (is.null(oma)) {
        oma <- c(4, 4, 4, 4)
        if (!is.null(main))
            oma[3] <- 6
    }
    opar <- par(mfrow = c(nc.y, nc.x), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar))
    for (i in if (row1attop)
        1:nc.y
    else nc.y:1) for (j in 1:nc.x) {
        localPlot(x[, j], y[, i], xlab = "", ylab = "", axes = FALSE,
            type = "n", ...)
        if (i == j || i < j || i > j ) {
            box()
            if (i == 1 && (!(j%%2)))
                localAxis(1 + 2 * row1attop, x[, j], y[, i],
                          xlab=xlabels[j], ylab=ylabels[i],
                  ...)
            if (i == nc.y && (j%%2))
                localAxis(3 - 2 * row1attop, x[, j], y[, i],
                          xlab=xlabels[j], ylab=ylabels[i],
                  ...)
            if (j == 1 && (!(i%%2) ))
                localAxis(2, x[, j], y[, i],
                          xlab=xlabels[j], ylab=ylabels[i],
                          ...)
            if (j == nc.x && (i%%2))
                localAxis(4, x[, j], y[, i],
                          xlab=xlabels[j], ylab=ylabels[i],
                          ...)
            mfg <- par("mfg")

            localPanel(as.vector(x[, j]), as.vector(y[,i]), ...)
            if (any(par("mfg") != mfg))
                stop("the 'panel' function made a new plot")
        }
        else par(new = FALSE)
    }
    if (!is.null(main)) {
        font.main <- if ("font.main" %in% nmdots)
            dots$font.main
        else par("font.main")
        cex.main <- if ("cex.main" %in% nmdots)
            dots$cex.main
        else par("cex.main")
        mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
    }
    invisible(NULL)
}
