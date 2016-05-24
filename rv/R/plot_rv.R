
plot.rv <- function (x, y=NULL, ...)
{ 
  .plot.default.rv(x, y, ...)
}

plot.rvsummary <- function (x, y=NULL, ...)
{ 
  .plot.default.rv(x, y, ...)
}

# ========================================================================
# plot.xy.rv  - 
# ========================================================================
#

.plot.xy.rv <- function (xy, type, pch = par("pch"), lty = par("lty"), col = par("col"), bg = NA, cex = 1, lwd = par("lwd"), ...)
{
  ##if (is.null(xy$rv)) {
  ##  return(.Internal(plot.xy(xy, type, pch, lty, col, bg, cex, lwd, ...)))
  ##}
  args <- c(xy, list(type=type, pch=pch, lty=lty, col=col, bg=bg, cex=cex, lwd=lwd, ...))
  typego <- list(p="points.rv", "l"="points.rv", "b"="points.rv")
  if (type=="n") return(invisible(NULL))
  if (is.null(plotroutine <- typego[[type]])) {
    stop("Plot type '", type, "' not yet implemented for rv objects!")
  }
  do.call(plotroutine, args)
  invisible(NULL)
}

# ========================================================================
# xy.coords.rv - modified version of xy.coords to accommodate rvs
# ========================================================================

.xy.coords.rv <- function (x, y = NULL, xlab = NULL, ylab = NULL, log = NULL, recycle = FALSE)
{
    if (is.null(y)) {
        ylab <- xlab
        if (is.language(x)) {
            if (inherits(x, "formula") && length(x) == 3) {
                ylab <- deparse(x[[2]])
                xlab <- deparse(x[[3]])
                y <- eval(x[[2]], environment(x), parent.frame())
                x <- eval(x[[3]], environment(x), parent.frame())
            }
            else stop("invalid first argument")
        }
        else if (is.complex(x)) {
            y <- Im(x)
            x <- Re(x)
            xlab <- paste("Re(", ylab, ")", sep = "")
            ylab <- paste("Im(", ylab, ")", sep = "")
        }
        else if (is.matrix(x) || is.data.frame(x)) {
            x <- data.matrix(x)
            if (ncol(x) == 1) {
                xlab <- "Index"
                y <- x[, 1]
                x <- 1:length(y)
            }
            else {
                colnames <- dimnames(x)[[2]]
                if (is.null(colnames)) {
                  xlab <- paste(ylab, "[,1]", sep = "")
                  ylab <- paste(ylab, "[,2]", sep = "")
                }
                else {
                  xlab <- colnames[1]
                  ylab <- colnames[2]
                }
                y <- x[, 2]
                x <- x[, 1]
            }
        }
        else if (is.list(x) && !is.rvobj(x)) { #### This is the only change ####
            xlab <- paste(ylab, "$x", sep = "")
            ylab <- paste(ylab, "$y", sep = "")
            y <- x[["y"]]
            x <- x[["x"]]
        }
        else {
            if (is.factor(x)) 
                x <- as.numeric(x)
            xlab <- "Index"
            y <- x
            x <- seq(along = x)
        }
    }
    if (inherits(x, "POSIXt")) 
        x <- as.POSIXct(x)
    if (length(x) != length(y)) {
        if (recycle) {
            if ((nx <- length(x)) < (ny <- length(y))) 
                x <- rep(x, length.out = ny)
            else y <- rep(y, length.out = nx)
        }
        else stop("'x' and 'y' lengths differ")
    }
    if (length(log) && log != "") {
        log <- strsplit(log, NULL)[[1]]
        f <- function (x) ((Pr(x < 0)>0) & !rv.any.na(x))
        if ("x" %in% log && any(ii <- f(x))) {
            n <- as.integer(sum(ii))
            warning(sprintf(ngettext(n, "%d x value <= 0 omitted from logarithmic plot", 
                "%d x values <= 0 omitted from logarithmic plot"), 
                n), domain = NA)
            x[ii] <- NA
        }
        if ("y" %in% log && any(ii <- f(y))) {
            n <- as.integer(sum(ii))
            warning(sprintf(ngettext(n, "%d y value <= 0 omitted from logarithmic plot", 
                "%d y values <= 0 omitted from logarithmic plot"), 
                n), domain = NA)
            y[ii] <- NA
        }
    }
    return(list(x = as.double(x), y = as.double(y), xlab = xlab, ylab = ylab))
}


.plot.default.rv <- function (x, y = NULL, type = "p", xlim = NULL, ylim = NULL, log = "", main = NULL, sub = NULL, xlab = NULL, ylab = NULL, ann = par("ann"), axes = TRUE, frame.plot = axes, panel.first = NULL, panel.last = NULL, asp = NA, rvlwd = rvpar("rvlwd"),  rvcol=rvpar("rvcol"), rvpoint=rvpar("rvpoint"), rvlex=rvpar("rvlex"), ...) 
{
    localAxis <- function(..., col, bg, pch, cex, lty, lwd) Axis(...)
    localBox <- function(..., col, bg, pch, cex, lty, lwd) box(...)
    localWindow <- function(..., col, bg, pch, cex, lty, lwd) plot.window(...)
    localTitle <- function(..., col, bg, pch, cex, lty, lwd) title(...)
    xlabel <- if (!missing(x)) 
        deparse(substitute(x))
    ylabel <- if (!missing(y)) 
        deparse(substitute(y))
    xy <- .xy.coords.rv(x, y, xlabel, ylabel, log)
    xlab <- if (is.null(xlab)) 
        xy$xlab
    else xlab
    ylab <- if (is.null(ylab)) 
        xy$ylab
    else ylab
    xlim <- if (is.null(xlim)) {
        range(rvfiniterange(xy$x))
    } else xlim
    ylim <- if (is.null(ylim)) {
        range(rvfiniterange(xy$y))
    } else ylim
    plot.new()
    localWindow(xlim, ylim, log, asp, ...)
    panel.first
    .plot.xy.rv(xy, type, rvlwd = rvlwd,  rvcol=rvcol, rvpoint=rvpoint, rvlex=rvlex, ...)
    panel.last
    if (axes) {
        localAxis(x, side = 1, ...)
        localAxis(y, side = 2, ...)
    }
    if (frame.plot) 
        localBox(...)
    if (ann) 
        localTitle(main = main, sub = sub, xlab = xlab, ylab = ylab, 
            ...)
    invisible()
}

