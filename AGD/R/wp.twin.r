#'Superposes two worm plots
#'
#'Superposes two worm plots from GAMLSS fitted objects. This is a diagnostic
#'tool for comparing two solutions.
#'
#'This function is a customized version of the \code{wp()} function found in
#'the \code{gamlss} package. Function \code{wp.twin()} allows overplotting of
#'two worm plots, each in its own color.  The points of \code{obj1} are plotted
#'first, the points of \code{obj2} are superposed.  This twin worm plot provide
#'a visual assessment of the differences between the solutions.  Extra
#'arguments can be specified (e.g. \code{xvar}) that are passed down to the
#'\code{wp()} function of \code{gamlss} if specified.  The worm plot is a
#'detrended normal QQ-plot that highlight departures from normality.
#'
#'Argument \code{xvar} takes priority over \code{xvar.column}. The \code{xvar}
#'variable is cut into \code{n.iter} intervals with an equal number
#'observations and detrended normal QQ (i.e. worm) plots for each interval are
#'plotted.  This is a way of highlighting failures of the model within
#'different ranges of the explanatory variable.
#'
#'If \code{line=TRUE} and \code{n.inter>1}, the fitted coefficients from
#'fitting cubic polynomials to the residuals (within each x-variable interval)
#'can be obtain by e.g.  \code{coeffs<-wp.twin(model1,xvar=x,n.iner=9)}.  van
#'Buuren \emph{et al.} (2001) used these residuals to identify regions
#'(intervals) of the explanatory variable within which the model does not fit
#'adequately the data (called "model violation")
#'
#'@param obj1 a GAMLSS fitted object
#'@param obj2 an optional second GAMLSS fitted object
#'@param xvar the explanatory variable against which the worm plots will be
#'plotted
#'@param xvar.column the number referring to the column of \code{obj1$mu.x} and
#'\code{obj2$mu.x}.  If \code{xvar=NULL} then the explanatory variable is set
#'to \code{xvar=obj1$mu.x[,xvar.column]} respectively
#'\code{xvar=obj2$mu.x[,xvar.column]}.  The default is \code{xvar.column=2},
#'which selects the variable following the intercept (which is typically age in
#'most applications).
#'@param n.inter the number of intervals in which the explanatory variable
#'\code{xvar} will be cut. The default is 16.
#'@param show.given whether to show the x-variable intervals in the top of the
#'graph, default is \code{show.given=FALSE}
#'@param ylim.worm for multiple plots, this values is the y-variable limit,
#'default value is \code{ylim.worm=0.5}
#'@param line whether to plot the polynomial line in the worm plot, default
#'value is \code{line=FALSE}
#'@param cex the cex plotting parameter with default \code{cex=1}
#'@param col1 the color for the points of \code{obj1}. The default
#'\code{col="black"}
#'@param col2 the color for the points of \code{obj2}. The default
#'\code{col="orange"}
#'@param warnings a logical indicating whether warnings should be produced. The
#'default \code{warnings=FALSE}
#'@param \dots for extra arguments, \code{overlap}, \code{xlim.worm} or
#'\code{pch}
#'@return For multiple plots the \code{xvar} intervals and the coefficients of
#'the fitted cubic polynomials to the residuals (within each \code{xvar}
#'interval) are returned.
#'@author Stef van Buuren, using R code of Mikis Stasinopoulos and Bob Rigby
#'@seealso \code{\link{wp}}
#'@references Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models
#'for location scale and shape (GAMLSS) in R. \emph{Journal of Statistical
#'Software}, Vol. \bold{23}, Issue 7, Dec 2007,
#'\url{http://www.jstatsoft.org/v23/i07}.
#'
#'van Buuren and Fredriks M. (2001) Worm plot: simple diagnostic device for
#'modelling growth reference curves.  \emph{Statistics in Medicine}, \bold{20},
#'1259--1277. 
#'
#'van Buuren and Fredriks M. (2007) Worm plot to diagnose fit in quantile
#'regression.  \emph{Statistical Modelling}, \bold{7}, 4, 363--376.
#'
#'@keywords smooth
#'@examples
#'
#'library(gamlss)
#'data(abdom)
#'a <- gamlss(y~cs(x,df=1),sigma.fo=~cs(x,0),family=LO,data=abdom)
#'b <- gamlss(y~cs(x,df=3),sigma.fo=~cs(x,1),family=LO,data=abdom)
#'coeff1 <- wp.twin(a,b,line=TRUE)
#'coeff1
#'rm(a,b,coeff1)
#'
wp.twin <- function(obj1, obj2=NULL, xvar=NULL, xvar.column=2,
           n.inter=16,
           show.given=FALSE, ylim.worm=0.5,
           line=FALSE, cex=1, col1="black", col2="orange",
           warnings=FALSE, ...)
{
    coplot2 <- function (formula, data, given.values, panel = points, rows,
        columns, show.given = TRUE, col = par("fg"), pch = par("pch"),
        bar.bg = c(num = gray(0.8), fac = gray(0.95)), xlab = c(x.name,
            paste("Given :", a.name)), ylab = c(y.name, paste("Given :",
            b.name)), subscripts = FALSE, axlabels = function(f) abbreviate(levels(f)),
        number = 6, overlap = 0.5, xlim, ylim,
        overplot = FALSE, ...)
    {
        deparen <- function(expr) {
            while (is.language(expr) && !is.name(expr) && deparse(expr[[1L]])[1L] ==
                "(") expr <- expr[[2L]]
            expr
        }
        bad.formula <- function() stop("invalid conditioning formula")
        bad.lengths <- function() stop("incompatible variable lengths")
        getOp <- function(call) deparse(call[[1L]], backtick = FALSE)[[1L]]
        formula <- deparen(formula)
        if (!inherits(formula, "formula"))
            bad.formula()
        y <- deparen(formula[[2L]])
        rhs <- deparen(formula[[3L]])
        if (getOp(rhs) != "|")
            bad.formula()
        x <- deparen(rhs[[2L]])
        rhs <- deparen(rhs[[3L]])
        if (is.language(rhs) && !is.name(rhs) && getOp(rhs) %in%
            c("*", "+")) {
            have.b <- TRUE
            a <- deparen(rhs[[2L]])
            b <- deparen(rhs[[3L]])
        }
        else {
            have.b <- FALSE
            a <- rhs
        }
        if (missing(data))
            data <- parent.frame()
        x.name <- deparse(x)
        x <- eval(x, data, parent.frame())
        nobs <- length(x)
        y.name <- deparse(y)
        y <- eval(y, data, parent.frame())
        if (length(y) != nobs)
            bad.lengths()
        a.name <- deparse(a)
        a <- eval(a, data, parent.frame())
        if (length(a) != nobs)
            bad.lengths()
        if (is.character(a))
            a <- as.factor(a)
        a.is.fac <- is.factor(a)
        if (have.b) {
            b.name <- deparse(b)
            b <- eval(b, data, parent.frame())
            if (length(b) != nobs)
                bad.lengths()
            if (is.character(b))
                b <- as.factor(b)
            b.is.fac <- is.factor(b)
            missingrows <- which(is.na(x) | is.na(y) | is.na(a) |
                is.na(b))
        }
        else {
            missingrows <- which(is.na(x) | is.na(y) | is.na(a))
            b <- NULL
            b.name <- ""
        }
        number <- as.integer(number)
        if (length(number) == 0L || any(number < 1))
            stop("'number' must be integer >= 1")
        if (any(overlap >= 1))
            stop("'overlap' must be < 1 (and typically >= 0).")
        bad.givens <- function() stop("invalid 'given.values'")
        if (missing(given.values)) {
            a.intervals <- if (a.is.fac) {
                i <- seq_along(a.levels <- levels(a))
                a <- as.numeric(a)
                cbind(i - 0.5, i + 0.5)
            }
            else co.intervals(unclass(a), number = number[1L], overlap = overlap[1L])
            b.intervals <- if (have.b) {
                if (b.is.fac) {
                    i <- seq_along(b.levels <- levels(b))
                    b <- as.numeric(b)
                    cbind(i - 0.5, i + 0.5)
                }
                else {
                    if (length(number) == 1L)
                      number <- rep.int(number, 2)
                    if (length(overlap) == 1L)
                      overlap <- rep.int(overlap, 2)
                    co.intervals(unclass(b), number = number[2L],
                      overlap = overlap[2L])
                }
            }
        }
        else {
            if (!is.list(given.values))
                given.values <- list(given.values)
            if (length(given.values) != (if (have.b)
                2L
            else 1L))
                bad.givens()
            a.intervals <- given.values[[1L]]
            if (a.is.fac) {
                a.levels <- levels(a)
                if (is.character(a.intervals))
                    a.intervals <- match(a.intervals, a.levels)
                a.intervals <- cbind(a.intervals - 0.5, a.intervals +
                    0.5)
                a <- as.numeric(a)
            }
            else if (is.numeric(a)) {
                if (!is.numeric(a.intervals))
                    bad.givens()
                if (!is.matrix(a.intervals) || ncol(a.intervals) !=
                    2)
                    a.intervals <- cbind(a.intervals - 0.5, a.intervals +
                      0.5)
            }
            if (have.b) {
                b.intervals <- given.values[[2L]]
                if (b.is.fac) {
                    b.levels <- levels(b)
                    if (is.character(b.intervals))
                      b.intervals <- match(b.intervals, b.levels)
                    b.intervals <- cbind(b.intervals - 0.5, b.intervals +
                      0.5)
                    b <- as.numeric(b)
                }
                else if (is.numeric(b)) {
                    if (!is.numeric(b.intervals))
                      bad.givens()
                    if (!is.matrix(b.intervals) || ncol(b.intervals) !=
                      2)
                      b.intervals <- cbind(b.intervals - 0.5, b.intervals +
                        0.5)
                }
            }
        }
        if (any(is.na(a.intervals)) || (have.b && any(is.na(b.intervals))))
            bad.givens()
        if (have.b) {
            rows <- nrow(b.intervals)
            columns <- nrow(a.intervals)
            nplots <- rows * columns
            if (length(show.given) < 2L)
                show.given <- rep.int(show.given, 2L)
        }
        else {
            nplots <- nrow(a.intervals)
            if (missing(rows)) {
                if (missing(columns)) {
                    rows <- ceiling(round(sqrt(nplots)))
                    columns <- ceiling(nplots/rows)
                }
                else rows <- ceiling(nplots/columns)
            }
            else if (missing(columns))
                columns <- ceiling(nplots/rows)
            if (rows * columns < nplots)
                stop("rows * columns too small")
        }
        total.columns <- columns
        total.rows <- rows
        f.col <- f.row <- 1
        if (show.given[1L]) {
            total.rows <- rows + 1
            f.row <- rows/total.rows
        }
        if (have.b && show.given[2L]) {
            total.columns <- columns + 1
            f.col <- columns/total.columns
        }
        mar <- if (have.b)
            rep.int(0, 4)
        else c(0.5, 0, 0.5, 0)
        oma <- c(5, 6, 5, 4)
        if (have.b) {
            oma[2L] <- 5
            if (!b.is.fac)
                oma[4L] <- 5
        }
        if (a.is.fac && show.given[1L])
            oma[3L] <- oma[3L] - 1
        opar <- par(mfrow = c(total.rows, total.columns), oma = oma,
            mar = mar, xaxs = "r", yaxs = "r")
        on.exit(par(opar))
        if (!overplot) plot.new()
        if (missing(xlim))
            xlim <- range(as.numeric(x), finite = TRUE)
        if (missing(ylim))
            ylim <- range(as.numeric(y), finite = TRUE)
        pch <- rep(pch, length.out = nobs)
        col <- rep(col, length.out = nobs)
        do.panel <- function(index, subscripts = FALSE, id) {
            Paxis <- function(side, x) {
                if (nlevels(x)) {
                    lab <- axlabels(x)
                    axis(side, labels = lab, at = seq(lab), xpd = NA)
                }
                else Axis(x, side = side, xpd = NA)
            }
            istart <- (total.rows - rows) + 1
            i <- total.rows - ((index - 1)%/%columns)
            j <- (index - 1)%%columns + 1
            par(mfg = c(i, j, total.rows, total.columns))
            if (!overplot) plot.new()
            plot.window(xlim, ylim)
            if (any(is.na(id)))
                id[is.na(id)] <- FALSE
            if (any(id)) {
                grid(lty = "solid")
                if (subscripts)
                    panel(x[id], y[id], subscripts = id, col = col[id],
                      pch = pch[id], ...)
                else panel(x[id], y[id], col = col[id], pch = pch[id],
                    ...)
            }
            if ((i == total.rows) && (j%%2 == 0))
                Paxis(1, x)
            else if ((i == istart || index + columns > nplots) &&
                (j%%2 == 1))
                Paxis(3, x)
            if ((j == 1) && ((total.rows - i)%%2 == 0))
                Paxis(2, y)
            else if ((j == columns || index == nplots) && ((total.rows -
                i)%%2 == 1))
                Paxis(4, y)
            box()
        }
        if (have.b) {
            count <- 1
            for (i in 1L:rows) {
                for (j in 1L:columns) {
                    id <- ((a.intervals[j, 1] <= a) & (a <= a.intervals[j,
                      2]) & (b.intervals[i, 1] <= b) & (b <= b.intervals[i,
                      2]))
                    do.panel(count, subscripts, id)
                    count <- count + 1
                }
            }
        }
        else {
            for (i in 1L:nplots) {
                id <- ((a.intervals[i, 1] <= a) & (a <= a.intervals[i,
                    2]))
                do.panel(i, subscripts, id)
            }
        }
        mtext(xlab[1L], side = 1, at = 0.5 * f.col, outer = TRUE,
            line = 3.5, xpd = NA, font = par("font.lab"), cex = par("cex.lab"))
        mtext(ylab[1L], side = 2, at = 0.5 * f.row, outer = TRUE,
            line = 3.5, xpd = NA, font = par("font.lab"), cex = par("cex.lab"))
        if (length(xlab) == 1L)
            xlab <- c(xlab, paste("Given :", a.name))
        if (show.given[1L]) {
            par(fig = c(0, f.col, f.row, 1), mar = mar + c(3 + (!a.is.fac),
                0, 0, 0), new = TRUE)
            if (!overplot) plot.new()
            nint <- nrow(a.intervals)
            a.range <- range(a.intervals, finite = TRUE)
            plot.window(a.range + c(0.03, -0.03) * diff(a.range),
                0.5 + c(0, nint))
            rect(a.intervals[, 1], 1L:nint - 0.3, a.intervals[, 2],
                1L:nint + 0.3, col = bar.bg[if (a.is.fac)
                    "fac"
                else "num"])
            if (a.is.fac) {
                text(apply(a.intervals, 1L, mean), 1L:nint, a.levels)
            }
            else {
                Axis(a, side = 3, xpd = NA)
                axis(1, labels = FALSE)
            }
            box()
            mtext(xlab[2L], 3, line = 3 - a.is.fac, at = mean(par("usr")[1L:2]),
                xpd = NA, font = par("font.lab"), cex = par("cex.lab"))
        }
        else {
            mtext(xlab[2L], 3, line = 3.25, outer = TRUE, at = 0.5 *
                f.col, xpd = NA, font = par("font.lab"), cex = par("cex.lab"))
        }
        if (have.b) {
            if (length(ylab) == 1L)
                ylab <- c(ylab, paste("Given :", b.name))
            if (show.given[2L]) {
                par(fig = c(f.col, 1, 0, f.row), mar = mar + c(0,
                    3 + (!b.is.fac), 0, 0), new = TRUE)
                if (!overplot) plot.new()
                nint <- nrow(b.intervals)
                b.range <- range(b.intervals, finite = TRUE)
                plot.window(0.5 + c(0, nint), b.range + c(0.03, -0.03) *
                    diff(b.range))
                rect(1L:nint - 0.3, b.intervals[, 1], 1L:nint + 0.3,
                    b.intervals[, 2], col = bar.bg[if (b.is.fac)
                      "fac"
                    else "num"])
                if (b.is.fac) {
                    text(1L:nint, apply(b.intervals, 1L, mean), b.levels,
                      srt = 90)
                }
                else {
                    Axis(b, side = 4, xpd = NA)
                    axis(2, labels = FALSE)
                }
                box()
                mtext(ylab[2L], 4, line = 3 - b.is.fac, at = mean(par("usr")[3:4]),
                    xpd = NA, font = par("font.lab"), cex = par("cex.lab"))
            }
            else {
                mtext(ylab[2L], 4, line = 3.25, at = 0.5 * f.row,
                    outer = TRUE, xpd = NA, font = par("font.lab"),
                    cex = par("cex.lab"))
            }
        }
        if (length(missingrows)) {
            cat("\n", gettext("Missing rows"), ": ", missingrows,
                "\n", sep = "")
            invisible(missingrows)
        }
    }

    wp2 <- function (object, xvar = NULL, n.inter = 4, xcut.points = NULL,
        overlap = 0, xlim.all = 4, xlim.worm = 3.5, show.given = TRUE,
        line = TRUE, ylim.all = 12 * sqrt(1/length(fitted(object))),
        ylim.worm = 12 * sqrt(n.inter/length(fitted(object))), cex = 1,
        pch = 21, overplot = FALSE, mline=3, color=col("col"), ...)
    {
        panel.fun <- function(x, y, col = par("col"), pch = par("pch"),
            cex = par("cex"), col.smooth = "red", span = 2/3, iter = 3,
            ...) {
            qq <- as.data.frame(qqnorm(y, plot = FALSE))
            qq$y <- qq$y - qq$x
            grid(nx = NA, ny = NA, lwd = 2)
            points(qq$x, qq$y, pch = pch, col = col, bg = col, cex = cex)
            abline(0, 0, lty = 2, col = 1)
            abline(0, 1e+05, lty = 2, col = 1)
            yuplim <- 10 * sqrt(1/length(qq$y))
            level <- 0.95
            lz <- -xlim.worm
            hz <- xlim.worm
            dz <- 0.25
            z <- seq(lz, hz, dz)
            p <- pnorm(z)
            se <- (1/dnorm(z)) * (sqrt(p * (1 - p)/length(qq$y)))
            low <- qnorm((1 - level)/2) * se
            high <- -low
            {
                no.points <- length(qq$y)
                total.points <<- total.points + no.points
                no.mis <- sum(abs(qq$y) > ylim.worm)
                warning(paste("number of missing points from plot=", no.mis,
                    " out of ", no.points, "\n",sep=""))
                if (any(abs(qq$y) > ylim.worm))
                    warning("Some points are missed out ", "\n",
                      "increase the y limits using ylim.worm")
            }
            if (any(abs(qq$x) > xlim.worm)) {
                warning("Some points are missed out ", "\n", "increase the x limits using xlim.worm")
            }
            lines(z, low, lty = 2, lwd = 0.01)
            lines(z, high, lty = 2, lwd = 0.01)
            if (line == TRUE) {
                fit <- lm(qq$y ~ qq$x + I(qq$x^2) + I(qq$x^3))
                s <- spline(qq$x, fitted(fit))
                flags <- s$x > -2.5 & s$x < 2.5
                lines(list(x = s$x[flags], y = s$y[flags]), col = col.smooth,
                    lwd = 0.01)
                coef1 <- coef(fit)
                assign("coef1", coef1, envir = parent.frame(n = 3))
                assign("coefall", c(coefall, coef1), envir = parent.frame(n = 3))
            }
        }
        check.overlap <- function(interval) {
            if (!is.matrix(interval)) {
                stop(paste("The interval specified is not a matrix."))
            }
            if (dim(interval)[2] != 2) {
                stop(paste("The interval specified is not a valid matrix.\nThe number of columns should be equal to 2."))
            }
            crows = dim(interval)[1]
            for (i in 1:(crows - 1)) {
                if (!(abs(interval[i, 2] - interval[i + 1, 1]) <
                    1e-04)) {
                    interval[i + 1, 1] = interval[i, 2]
                }
            }
            return(interval)
        }
        get.intervals <- function(xvar, xcut.points) {
            if (!is.vector(xcut.points)) {
                stop(paste("The interval is not a vector."))
            }
            if (any((xcut.points < min(xvar)) | any(xcut.points >
                max(xvar)))) {
                stop(paste("The specified `xcut.points' are not within the range of the x: (",
                    min(xvar), " , ", max(xvar), ")"))
            }
            extra <- (max(xvar) - min(xvar))/1e+05
            int <- c(min(xvar), xcut.points, (max(xvar) + 2 * extra))
            ii <- 1:(length(int) - 1)
            r <- 2:length(int)
            x1 <- int[ii]
            xr <- int[r] - extra
            if (any(x1 > xr)) {
                stop(paste("The interval is are not in a increasing order."))
            }
            cbind(x1, xr)
        }
        if (!is.gamlss(object))
            stop(paste("This is not an gamlss object", "\n", ""))
        if (is.null(xvar)) {
            qq <- as.data.frame(qqnorm(resid(object), plot = FALSE))
            qq$y <- qq$y - qq$x
            level <- 0.95
            lz <- -xlim.all
            hz <- xlim.all
            dz <- 0.25
            z <- seq(lz, hz, dz)
            p <- pnorm(z)
            se <- (1/dnorm(z)) * (sqrt(p * (1 - p)/length(qq$y)))
            low <- qnorm((1 - level)/2) * se
            high <- -low
            if (any(abs(qq$y) > ylim.all)) {
                warning("Some points are missed out ", "\n", "increase the y limits using ylim.all")
            }
            if (any(abs(qq$x) > xlim.all)) {
                warning("Some points are missed out ", "\n", "increase the x limits using xlim.all")
            }
            plot(qq$x, qq$y, ylab = "Deviation", xlab = "",
                xlim = c(-xlim.all, xlim.all), ylim = c(-ylim.all,
                    ylim.all), cex = cex, pch = pch, bg = "wheat",
                )
            grid(lty = "solid")
            abline(0, 0, lty = 2, col = 2)
            abline(0, 1e+05, lty = 2, col = 2)
            lines(z, low, lty = 2)
            lines(z, high, lty = 2)
            if (line == TRUE) {
                fit <- lm(qq$y ~ qq$x + I(qq$x^2) + I(qq$x^3))
                s <- spline(qq$x, fitted(fit))
                flags <- s$x > -3 & s$x < 3
                lines(list(x = s$x[flags], y = s$y[flags]), col = color,
                    lwd = 0.01)
            }
        }
        else {
            w <- object$weights
            if (all(trunc(w) == w))
                xvar <- rep(xvar, w)
            if (is.null(xcut.points)) {
                given.in <- co.intervals(xvar, number = n.inter,
                    overlap = overlap)
                if (overlap == 0)
                    given.in <- check.overlap(given.in)
            }
            else {
                given.in <- get.intervals(xvar, xcut.points)
            }
            total.points <- 0
            coefall <- coef1 <- NULL
            y <- resid(object)
            x <- resid(object)
            coplot2(y ~ x | xvar, given.values = given.in, panel = panel.fun,
                ylim = c(-ylim.worm, ylim.worm), xlim = c(-xlim.worm,
                    xlim.worm), ylab = "Deviation", xlab = "",
                show.given = show.given, bg = "wheat", pch = pch,
                cex = cex, bar.bg = c(num = "light blue"),
                overplot = overplot, col=color, ...)
            mtext(paste(object$call[c(-1,-length(object$call))],collapse=" "), side=1, line=mline, col=color)

            if (overlap == 0) {
                if (total.points != length(y))
                    warning("the total number of points in the plot is not equal \n to the number of observations in y \n")
            }
        }
        if (!is.null(xvar) & line) {
            mcoef <- matrix(coefall, ncol = 4, byrow = TRUE)
            out <- list(classes = given.in, coef = mcoef)
        }
    }

    is.valid <- function(xvar.column, obj) {
      return (xvar.column >= 1 & xvar.column <= ncol(obj$mu.x)) 
    }
    
    # ---- ENTRY POINT

    call <- match.call()
    
    if (!is.gamlss(obj1))
        stop(paste("Argument obj1 is not an gamlss object", "\n", ""))
    if ((!is.gamlss(obj2)) & (!is.null(obj2)))
        stop(paste("Argument obj2 is not an gamlss object", "\n", ""))

    expr1 <- expression(
        wp2(obj1, xvar=xvar, n.inter=n.inter,
            show.given=show.given, ylim.worm=ylim.worm,
            cex=cex, color=col1, col.smooth=col1,
            overplot=FALSE, line=line, mline=3, ...)
    )
    expr2 <- expression(
        wp2(obj2, xvar=xvar, n.inter=n.inter,
            show.given=show.given, ylim.worm=ylim.worm,
            cex=cex, color=col2, col.smooth=col2,
            overplot=TRUE, line=line, mline=4, ...)
    )

    if (is.null(call$xvar)) {
      if (n.inter > 1) {
        if (!is.valid(xvar.column, obj1)) {
          n.inter <- 1
          if (warnings) warning("Parameter n.inter changed to 1 because xvar.column is out-of-range.")
        } else {
          xvar <- obj1$mu.x[,xvar.column]
          if (min(xvar)==max(xvar)) stop("xvar is constant")
        }
      }
    }
    if (!warnings)
      suppressWarnings(eval(expr1))
    else
      eval(expr1)

    if (!is.null(call$obj2)) {
      if (is.null(call$xvar)) {
        if (n.inter > 1) {
          if (!is.valid(xvar.column, obj2)) {
            n.inter <- 1
            if (warnings) warning("Parameter n.inter changed to 1 because xvar.column is out-of-range.")
          } else {
            xvar <- obj2$mu.x[,xvar.column]
            if (min(xvar)==max(xvar)) stop("xvar is constant")
          }
        }
      }

      if (!warnings)
        suppressWarnings(eval(expr2))
      else
        eval(expr2)
    }
}

