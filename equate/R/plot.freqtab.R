#' Plotting Frequency Distributions
#' 
#' This function plots univariate and bivariate frequency tables of class
#' \dQuote{\code{\link{freqtab}}}.
#' 
#' For the points method, a scatterplot for \code{x} is added to the current
#' opened plot.
#' 
#' For the plot method, when \code{x} is univariate, i.e, having 2 columns, a
#' frequency plot is created for \code{x}. When \code{x} is bivariate, e.g.,
#' coming from a single group equating design or one form of a nonequivalent
#' groups design, a scatterplot is produced with frequency plots for the
#' marginal distributions.
#' 
#' \code{y} is used to superimpose lines, e.g., smoothed frequencies, over the
#' (marginal) frequencies of \code{x}.
#' 
#' Colors must be specified using \code{xcol} and \code{ycol}. When \code{ycol}
#' is missing, a vector of colors is created using \code{rainbow(ncol(y))}.
#' 
#' @aliases plot.freqtab points.freqtab bfreqplot ufreqplot
#' @param x univariate or bivariate score distribution of class
#' \dQuote{\code{\link{freqtab}}}.
#' @param y either an object of class \dQuote{\code{freqtab}}, where
#' frequencies will be extracted, or a vector or matrix of frequencies, to be
#' added to the plot of \code{x}. See below for details.
#' @param xcol,ycol colors used in plotting \code{x} and \code{y}.
#' @param pch plotting symbol used to plot bivariate points.
#' @param ylty line type used to plot frequencies in \code{y}.
#' @param xlab label for the x axis.
#' @param addlegend logical indicating whether or not a legend should be added.
#' @param legendtext character vector of text to be passed to the \code{legend}
#' argument of the \code{legend} function, defaulting to column names used in
#' \code{y}.
#' @param ds,dm integers for the scaling and center of the RGB density values,
#' with defaults of 50 and 100. These are used to convert the observed counts
#' in \code{x} to the [0, 255] range of RGB values.
#' @param \dots further arguments passed to or from other methods, such as
#' graphical parameters besides \code{col}, \code{type}, and \code{pch}.
#' @return The univariate option produces a single line plot of \code{type =
#' "h"}. Frequencies from \code{y} are then superimposed. The bivariate option
#' produces a scatterplot with a marginal frequency plot for each distribution.
#' @author Anthony Albano \email{tony.d.albano@@gmail.com}
#' @seealso \code{\link{plot.table}}, \code{\link{plot.equate}},
#' \code{\link{lines}}, \code{\link{points}}
#' @keywords misc
#' @examples
#' 
#' x <- freqtab(KBneat$x, scales = list(0:36, 0:12))
#' plot(x)
#' 
#' xs <- loglinear(x, degrees = c(4, 1),
#'   stepup = TRUE, showWarnings = FALSE)
#' plot(x, xs, lwd = 2)
#' 
#' @export
plot.freqtab <- function(x, y = NULL, xcol = 1,
  ycol, pch = 16, ylty = 1, xlab = "Total Test",
  addlegend = !missing(y), legendtext, ...) {
  
  nx <- margins(x)
  if (nx == 1)
    ufreqplot(x, y, xcol, ycol, ylty, xlab,
      horiz = FALSE, addlegend = addlegend,
      legendtext = legendtext, ...)
  else if (nx == 2)
    bfreqplot(x, y, xcol, ycol, pch, ylty, xlab,
      addlegend = addlegend,
      legendtext = legendtext, ...)
  else stop("'x' must be either univariate or bivariate")
}

#----------------------------------------------------------------
# Internal univariate plot

ufreqplot <- function(x, y = NULL, xcol = 1, ycol,
  ylty = 1, xlab = "Total Test", ylab = "Count",
  horiz = FALSE, addlegend = FALSE,
  legendtext, ...) {
  
  x <- as.data.frame(x)
  if (!is.null(y)) {
    if (is.freqtab(y))
      y <- cbind(y[, 2])
    else
      y <- cbind(y)
    if (missing(ycol))
      ycol <- rainbow(ncol(y))
  }
  
  if (horiz) {
    plot.default(round(range(0, x[, 2], y)),
      range(x[, 1]), type = "n", xlab = xlab,
      ylab = ylab, ...)
    segments(rep(0, nrow(x)), x[, 1], x[, 2],
      col = xcol)
    if (!is.null(y))
      matlines(y, x[, 1], col = ycol,
        lty = ylty, ...)
  }
  else {
    plot.default(range(x[, 1]),
      round(range(0, x[, 2], y)), type = "n",
      xlab = xlab, ylab = ylab, ...)
    segments(x[, 1], y0 = rep(0, nrow(x)),
      y1 = x[, 2], col = xcol)
    if (!is.null(y))
      matlines(x[, 1], y, col = ycol,
        lty = ylty, ...)
  }
  
  if (addlegend & !is.null(y)) {
    if (missing(legendtext))
      legendtext <- if (is.null(colnames(y)))
        1:ncol(y) else colnames(y)
    legend("topright", legend = legendtext,
      lty = ylty, col = ycol, bty = "n")
  }
}

#----------------------------------------------------------------
# Internal bivariate plot

bfreqplot <- function(x, y = NULL, xcol = 1,
  ycol, pch = 16, ylty = 1, xlab = "Total Test",
  ylab = "Anchor Test", addlegend = FALSE,
  legendtext, ...) {
  
  xtab <- margin(x)
  xvtab <- margin(x, 2)
  xd <- as.data.frame(x)
  
  if (!is.null(y)) {
    if (is.freqtab(y))
      y <- cbind(c(y))
    if (missing(ycol))
      ycol <- rainbow(ncol(y))
    ytab <- apply(y, 2, function(z)
      tapply(z, xd[, 1], sum))
    yvtab <- apply(y, 2, function(z)
      tapply(z, xd[, 2], sum))
  }
  else ytab <- yvtab <- NULL
  
  reset.par <- par(no.readonly = TRUE)
  nf <- layout(matrix(c(2, 4, 1, 3), 2, 2,
    byrow = TRUE), c(3, 1), c(1, 3), TRUE)
  par(mar = c(4, 4, 1, 1))
  plot(range(xtab), range(xvtab), type = "n",
    xlab = xlab, ylab = ylab, ...)
  points(x, xcol = xcol, pch = pch)
  
  par(mar = c(0, 4, 1, 1))
  ufreqplot(xtab, ytab, xcol, ycol, ylty,
    xlab = "", ylab = "", xaxt = "n", bty = "n")
  
  par(mar = c(4, 0, 1, 1))
  ufreqplot(xvtab, yvtab, xcol, ycol, ylty,
    xlab = "", ylab = "", yaxt = "n", bty = "n",
    horiz = TRUE)
  
  if (addlegend & !is.null(y)) {
    par(mar = c(0, 0, 0, 0))
    plot(0, 0, type = "n", bty = "n", xaxt = "n",
      yaxt = "n")
    if (missing(legendtext))
      legendtext <- if (is.null(colnames(y)))
        1:ncol(y) else colnames(y)
    legend("bottomleft", legend = legendtext,
      lty = ylty, col = ycol, bty = "n")
  }
  
  par(reset.par)
}

#----------------------------------------------------------------
# Points method

# @describeIn plot.freqtab \code{\link{points}} method for
# \dQuote{\code{\link{freqtab}}} objects.
#' @rdname plot.freqtab
#' @export
points.freqtab <- function(x, xcol = 1, pch = 16,
  ds = 50, dm = 100, ...) {
  
  x <- as.data.frame(x)
  if (ncol(x) != 3)
    stop("'x' must be a bivariate frequency table")
  
  index <- as.logical(x[, 3])
  xpoints <- x[index, 1]
  vpoints <- x[index, 2]
  if (sd(x[index, 3]) > 0)
    dens <- pmax(0, pmin(255,
      scale(x[index, 3]) * ds + dm))
  else dens <- rep(150, sum(index))
  rgbcol <- col2rgb(xcol)
  ptcol <- rgb(rgbcol[1], rgbcol[2], rgbcol[3],
    dens, maxColorValue = 255)
  points(xpoints, vpoints, col = ptcol,
    pch = pch, ...)	
}