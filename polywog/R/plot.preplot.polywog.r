##' @S3method plot preplot.polywog
plot.preplot.polywog <- function(x, auto.set.par = TRUE,
                                 FUN3D = c("contour", "filled.contour",
                                 "wireframe", "persp3d"),
                                 control.plot = list(),
                                 ...)
{
    xvars <- attr(x, "xvars")
    xcol <- attr(x, "xcol")
    whichFactors <- sapply(xcol, function(i) is.factor(x[, i]) ||
                           all(x[, i] %in% c(0, 1)))
    if (all(whichFactors == c(FALSE, TRUE))) {
        ## Reorder if only the second variable is categorical
        xvars <- rev(xvars)
        xcol <- rev(xcol)
    }

    ## Possibilities:
    ##   Two variables, at least one categorical: multiple plots broken up by
    ##     levels of the first
    ##   Two variables, both continuous: contour plot
    ##   Single variable, categorical: box plot (of sorts)
    ##   Single variable, continuous: scatterplot
    if (length(whichFactors) == 2 && any(whichFactors)) {
        ## Two variables, at least one categorical

        ## Set up the plot
        col <- xcol[1]
        nf <- length(unique(x[, col]))
        if (auto.set.par) {
            mfrow <- ceiling(sqrt(nf))
            mfcol <- ceiling(nf / mfrow)
            if (!exists("..op")) {
                ..op <- par(mfrow = c(mfrow, mfcol))
                on.exit(par(..op))
            } else {
                par(mfrow = c(mfrow, mfcol))
            }
        }

        ## Plot the relationship at each value of the factor/binary variable
        for (i in seq_len(nf)) {
            vali <- unique(x[, col])[i]
            xx <- x[x[, col] == vali, , drop = FALSE]
            attr(xx, "xcol") <- xcol[2]
            attr(xx, "xvars") <- xvars[2]
            control.plot$main <- paste(xvars[1], "=", vali)
            plot(xx, auto.set.par = auto.set.par, control.plot = control.plot)
        }
    } else if (length(whichFactors) == 2) {
        ## Two variables, both continuous

        ## Take user input about which 3D plotting function to use
        FUN3D <- match.arg(FUN3D)
        if (FUN3D == "wireframe" && !require("lattice")) {
            stop("'lattice' package required for FUN3D = \"wireframe\"")
        } else if (FUN3D == "persp3d" && !require("rgl")) {
            stop("'rgl' package required for FUN3D = \"persp3d\"")
        }

        ## Extract data
        var1 <- unique(x[, xcol[1]])
        var2 <- unique(x[, xcol[2]])
        fit <- matrix(x$fit, nrow = length(var1))

        ## Make plot
        if (FUN3D == "wireframe") {
            cl <- list(x = fit, row.values = var1, column.values = var2)
            cl <- c(cl, control.plot)
            if (is.null(cl$xlab))
                cl$xlab <- xvars[1]
            if (is.null(cl$ylab))
                cl$ylab <- xvars[2]
            if (is.null(cl$zlab))
                cl$zlab <- "fitted value"
            print(do.call(FUN3D, cl))
        } else if (FUN3D == "persp3d") {
            cl <- list(x = var1, y = var2, z = fit)
            cl <- c(cl, control.plot)
            if (is.null(cl$xlab))
                cl$xlab <- xvars[1]
            if (is.null(cl$ylab))
                cl$ylab <- xvars[2]
            if (is.null(cl$zlab))
                cl$zlab <- "fitted value"
            do.call(FUN3D, cl)
            if (attr(x, "interval")) {
                ## Confidence regions
                upr <- matrix(x$upr, nrow = length(var1))
                lwr <- matrix(x$lwr, nrow = length(var1))
                persp3d(x = var1, y = var2, z = upr,
                        col = "gray70", alpha = 0.7, add = TRUE)
                persp3d(x = var1, y = var2, z = lwr,
                        col = "gray70", alpha = 0.7, add = TRUE)
            }
        } else {
            cl <- list(z = fit, x = var1, y = var2)
            cl <- c(cl, control.plot)
            if (is.null(cl$xlab))
                cl$xlab <- xvars[1]
            if (is.null(cl$ylab))
                cl$ylab <- xvars[2]
            do.call(FUN3D, cl)
        }
    } else if (whichFactors[1]) {
        ## One variable, categorical

        ## Manually set up a "boxplot" with fitted values and bars for
        ## confidence levels
        boxStats <- list()
        boxStats$stats <- matrix(x$fit, nrow = 5, ncol = nrow(x), byrow =
                                 TRUE)
        if (attr(x, "interval")) {
            boxStats$stats[1, ] <- x$lwr
            boxStats$stats[5, ] <- x$upr
        }
        boxStats$n <- rep(1, nrow(x))
        boxStats$conf <- boxStats$stats[c(1, 5), ]
        boxStats$out <- numeric(0)
        boxStats$group <- numeric(0)
        boxStats$names <- as.character(x[, xcol])

        cl <- list(z = boxStats)
        cl <- c(cl, control.plot)
        if (is.null(cl$xlab))
            cl$xlab <- xvars[1]
        if (is.null(cl$ylab))
            cl$ylab <- "fitted value"
        do.call(bxp, cl)
    } else {
        ## One variable, continuous

        cl <- list(x = x[, xcol], y = x$fit, type = "l")
        cl <- c(cl, control.plot)
        if (is.null(cl$xlab))
            cl$xlab <- xvars[1]
        if (is.null(cl$ylab))
            cl$ylab <- "fitted value"
        if (is.null(cl$ylim))
            cl$ylim <- c(min(x$fit, x$upr, x$lwr), max(x$fit, x$upr, x$lwr))
        do.call(plot, cl)
        if (attr(x, "interval")) {
            lines(x[, xcol], x$lwr, lty = 2)
            lines(x[, xcol], x$upr, lty = 2)
        }
    }

    invisible(x)
}
