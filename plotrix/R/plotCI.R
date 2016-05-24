plotCI <- function (x, y = NULL, uiw, liw = uiw, ui = NULL, li = NULL, 
    err = "y", sfrac = 0.01, gap = 0, slty = par("lty"), add = FALSE, 
    scol = NULL, pt.bg = par("bg"), ...) 
{
    arglist <- list(...)
    if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (is.null(y)) {
        if (is.null(x)) 
            stop("both x and y NULL")
        y <- as.numeric(x)
        x <- seq(along = x)
    }
    if (missing(uiw) && (is.null(ui) || is.null(li))) 
        stop("must specify either relative limits or both lower and upper limits")
    if (!missing(uiw)) {
        if (err == "y") 
            z <- y
        else z <- x
        ui <- z + uiw
        li <- z - liw
    }
    if (is.null(arglist$xlab)) 
        arglist$xlab <- deparse(substitute(x))
    if (is.null(arglist$ylab)) 
        arglist$ylab <- deparse(substitute(y))
    if (err == "y" && is.null(arglist$ylim)) 
        arglist$ylim <- range(c(y, ui, li), na.rm = TRUE)
    if (err == "x" && is.null(arglist$xlim)) 
        arglist$xlim <- range(c(x, ui, li), na.rm = TRUE)
    if (missing(scol)) {
        if (!is.null(arglist$col)) 
            scol <- arglist$col
        else scol <- par("col")
    }
    plotpoints <- TRUE
    if (!is.null(arglist$pch) && is.na(arglist$pch)) {
        arglist$pch <- 1
        plotpoints <- FALSE
    }
    if (!add) 
        do.call("plot", c(list(x, y, type = "n"), clean.args(arglist, 
            plot)))
    if (gap == TRUE) 
        gap <- 0.01
    ul <- c(li, ui)
    pin <- par("pin")
    usr <- par("usr")
    x.to.in <- pin[1]/diff(usr[1:2])
    y.to.in <- pin[2]/diff(usr[3:4])
    if (err == "y") {
        gap <- rep(gap, length(x)) * diff(par("usr")[3:4])
        smidge <- par("fin")[1] * sfrac
        nz <- abs(li - pmax(y - gap, li)) * y.to.in > 0.001
        scols <- rep(scol, length.out = length(x))[nz]
        arrow.args <- c(list(lty = slty, angle = 90, length = smidge, 
            code = 1, col = scols), clean.args(arglist, arrows, 
            exclude.other = c("col", "lty", "axes")))
        do.call("arrows", c(list(x[nz], li[nz], x[nz], pmax(y - 
            gap, li)[nz]), arrow.args))
        nz <- abs(ui - pmin(y + gap, ui)) * y.to.in > 0.001
        scols <- rep(scol, length.out = length(x))[nz]
        arrow.args$col <- scols
        do.call("arrows", c(list(x[nz], ui[nz], x[nz], pmin(y + 
            gap, ui)[nz]), arrow.args))
    }
    else if (err == "x") {
        gap <- rep(gap, length(x)) * diff(par("usr")[1:2])
        smidge <- par("fin")[2] * sfrac
        nz <- abs(li - pmax(x - gap, li)) * x.to.in > 0.001
        scols <- rep(scol, length.out = length(x))[nz]
        arrow.args <- c(list(lty = slty, angle = 90, length = smidge, 
            code = 1, col = scols), clean.args(arglist, arrows, 
            exclude.other = c("col", "lty", "axes")))
        do.call("arrows", c(list(li[nz], y[nz], pmax(x - gap, 
            li)[nz], y[nz]), arrow.args))
        nz <- abs(ui - pmin(x + gap, ui)) * x.to.in > 0.001
        scols <- rep(scol, length.out = length(x))[nz]
        arrow.args$col <- scols
        do.call("arrows", c(list(ui[nz], y[nz], pmin(x + gap, 
            ui)[nz], y[nz]), arrow.args))
    }
    if (plotpoints) 
        do.call("points", c(list(x, y, bg = pt.bg), clean.args(arglist, 
            points, exclude.other = c("xlab", "ylab", "xlim", 
                "ylim", "axes"))))
    invisible(list(x = x, y = y))
}
