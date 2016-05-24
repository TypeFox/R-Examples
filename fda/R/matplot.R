matplot <- function(x, ...) {
    UseMethod('matplot')
}

matplot.matrix <- function(x, ...){
    fda::matplot.default(x, ...)
}

matplot.numeric <- function(x, ...){
    fda::matplot.default(x, ...)
}

matplot.default <- function(x, y, type = "p", lty = 1:5, lwd = 1,
    lend = par("lend"), pch = NULL, col = 1:6, cex = NULL, bg = NA,
    xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, ..., add = FALSE,
    verbose = getOption("verbose")){
#
  if(is.null(xlab))
      xlab <- deparse(substitute(x))
  if(is.null(ylab))
      ylab <- deparse(substitute(y))
#
  graphics::matplot(x, y, type = type, lty = lty, lwd = lwd,
    lend = lend, pch = pch, col = col, cex = cex, bg = bg,
    xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ..., add = add,
    verbose = getOption("verbose"))
#
}

matplot.Date <- function(x, y, type = "p", lty = 1:5, lwd = 1,
    lend = par("lend"), pch = NULL, col = 1:6, cex = NULL, bg = NA,
    xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, ..., add = FALSE,
    verbose = getOption("verbose")){
#
  if(is.null(xlab))
      xlab <- deparse(substitute(x))
  if(is.null(ylab))
      ylab <- deparse(substitute(y))
#
  matplot.POSIXct(x, y, type = type, lty = lty, lwd = lwd,
    lend = lend, pch = pch, col = col, cex = cex, bg = bg,
    xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ..., add = add,
    verbose = getOption("verbose"))
}

matplot.POSIXct <- function(x, y, type = "p", lty = 1:5, lwd = 1,
    lend = par("lend"), pch = NULL, col = 1:6, cex = NULL, bg = NA,
    xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, ..., add = FALSE,
    verbose = getOption("verbose")){
#
  if(is.null(xlab))
      xlab <- deparse(substitute(x))
  if(is.null(ylab))
      ylab <- deparse(substitute(y))
  if(is.null(xlim))
      xlim <- range(x, na.rm=TRUE)
  if(is.null(ylim))
      ylim <- range(y, na.rm=TRUE)
#
  if(!add){
    plot(range(x), range(y), type='n', cex=cex, bg=bg, xlab=xlab, ylab=ylab,
         xlim=xlim, ylim=ylim, ...)
    out <- matlines(x, y, type=type, cex=cex, ..., verbose=verbose)
    return(out)
  }
  graphics::matplot(x, y, type = type, lty = lty, lwd = lwd,
    lend = lend, pch = pch, col = col, cex = cex, bg = bg,
    xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ..., add = add,
    verbose = getOption("verbose"))
}
