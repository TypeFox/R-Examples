## Copyright (c) 2008-2016, James P. Howard, II <jh@jameshoward.us>
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
##     Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##
##     Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#' @title Waterfall Plot
#'
#' @description
#' Creates a waterfall plot with vertical or horizontal bars.
#'
#' @param height a vector of values describing the height of the bars
#' that make up the plot.  Matrices are not supported.
#'
#' @param width optional vector of bar widths.  Re-cycled to length the
#' number of bars drawn.  Specifying a single value will have no visible
#' effect unless 'xlim' is specified.
#'
#' @param space
#'    the amount of space (as a fraction of the average bar width)
#'    left before each bar.  May be given as a single number or one number
#'    per bar.  If not given explicitly, it defaults to 0.2.
#'
#' @param names.arg a vector of names to be plotted below each bar If
#' this argument is omitted, then the names are taken from the 'names'
#' attribute of 'height.'
#'
#' @param horiz a logical value.  If 'FALSE', the bars are drawn
#' vertically with the first bar to the left.  If 'TRUE', the bars are
#' drawn horizontally with the first at the bottom.
#'
#' @param density a vector giving the density of shading lines, in lines
#' per inch, for the bars or bar components. The default value of 'NULL'
#' means that no shading lines are drawn. Non-positive values of
#' 'density' also inhibit the drawing of shading lines.
#'
#' @param angle the slope of shading lines, given as an angle in degrees
#' (counter-clockwise), for the bars or bar components.
#'
#' @param col a vector of colors for the bars or bar components. By
#' default, grey is used.
#'
#' @param border the color to be used for the border of the bars. Use
#' 'border = NA' to omit borders.  If there are shading lines, 'border =
#' TRUE' means use the same colour for the border as forr the shading
#' lines.
#'
#' @param main overall title for the plot.
#'
#' @param sub subtitle for the plot.
#'
#' @param xlab a label for the x-axis.
#' @param ylab a label for the y-axis.
#' @param xlim limits for the x-axis.
#' @param ylim limits for the y-axis.
#' @param xpd logical. Should bars be allowed to outside region?
#'
#' @param axes logical.  If 'TRUE', a vertical (or horizontal, if
#' 'horiz' is true) axis is drawn.
#'
#' @param axisnames logical.  If 'TRUE', and if there are 'names.arg'
#' (see above), the other axis is drawn (with 'lty=0') and labeled.
#'
#' @param cex.axis expansion factor for numeric axis labels.
#' @param cex.names expansion factor for axis names (bar labels).
#' @param plot logical.  If 'FALSE', nothing is plotted.
#'
#' @param axis.lty the graphics parameter 'lty' applied to the axis and
#' tick marks of the categorical (default horizontal) axis.  Note that by
#' default the axis is suppressed.
#'
#' @param offset initial offset relative to the x axis.  The value
#' serves as the logical starting point for the first column and any
#' summary column.  Defaults to 0.
#'
#' @param add logical specifying if bars should be added to an already
#' existing plot; defaults to 'FALSE'.
#'
#' @param summary create a summary column.  A summary column provides a
#' final sum column showing the relative change from the offset.  If a
#' summary is requested and names.arg is set, the names.arg vector must
#' include one extra entry with the summary column's name.  Defaults to
#' FALSE.
#'
#' @param rev reverse the order of columns?  Defaults to FALSE.
#'
#' @param level.lines if FALSE, the lines connecting adjacent boxes are
#' ommitted from the display.
#'
#' @param ...  arguments to be passed to other methods.  For the default
#' method these can include further arguments (such as 'axes', 'asp' and
#' 'main') and graphical parameters (see 'par') which are passed to
#' 'plot.window()', 'title()' and 'axis'.
#'
#' @details
#'
#' This function closely mimics the \link{barplot} interface, but
#' provides a type of chart called a waterfall plot, showing how multiple
#' subvalues contribute to a total sum.
#'
#' This is a generic function, it currently only has a default method. A
#' formula interface may be added eventually.
#'
#' @return A numeric vector say 'mp', giving the coordinates of
#' \emph{all} the bar midpoints drawn, useful for adding to the graph.
#'
#' @references
#'
#' Andrew Jaquith, \emph{Security Metrics: Replacing Fear,
#' Uncertainty, and Doubt} (Boston: Addison-Wesley Professional, 2007),
#' 170-172.
#'
#' Ethan M. Rasiel, \emph{The McKinsey Way: Using the Techniques of the
#' World's Top Strategic Consultants to Help You and Your Business} (New
#' York: McGraw-Hill, 1999), 113-118.
#'
#' @examples
#' data(rasiel)
#' waterfallplot(rasiel$value, names.arg=rasiel$label)
#'
#' @importFrom graphics par
#' @importFrom graphics plot.new
#' @importFrom graphics plot.window
#' @importFrom graphics rect
#' @importFrom graphics lines
#' @importFrom graphics axis
#' @importFrom graphics title
#'
#' @export
waterfallplot <-
    function (height, width = 1, space = NULL, names.arg = NULL,
              horiz = FALSE, density = NULL, angle = 45, col = NULL, border = par("fg"),
              main = NA, sub = NA, xlab = NULL, ylab = NULL, xlim = NULL,
              ylim = NULL, xpd = TRUE, axes = TRUE, axisnames = TRUE, cex.axis = par("cex.axis"),
              cex.names = par("cex.axis"), plot = TRUE, axis.lty = 0, offset = 0,
              add = FALSE, summary = FALSE, rev = FALSE, level.lines = TRUE, ...)
{
	##  Set up the initial environment
	respaxis <- 2
	expaxis <- 1
	l <- length(height)
	if (summary == TRUE) {
		height[l + 1] = -sum(height)
		l = l + 1
	}
	if (rev == TRUE) {
		height = rev(height)
		if (summary == TRUE)
			height = -height
	}
	if (!is.null(names.arg) && length(names.arg) != l)
		stop("incorrect number of names")
	if (!is.null(names(height)))
		names.arg <- names(height)
	density.vec <- rep(density, l, length.out = l)
	angle.vec <- rep(angle, l, length.out = l)
	if (is.null(col))
		col <- "grey"
	col.vec <- rep(col, l, length.out = l)
	border.vec <- rep(border, l, length.out = l)

	##  Assemble the baseline and topline vectors.  Default the spacing
	##  between bars to 1/5 of the average bar width.  Create the same
	##  initial offset as present in barplot()
	width.vec <- rep(width, l, length.out = l)
	if (is.null(space))
		space <- 0.2
	space.vec <- rep(space, l + 1, length.out = l + 1)
	space.vec <- space.vec * mean(width.vec)
	leftline <- rightline <- topline <- baseline <- rep(offset, l)
	leftline[1] <- space.vec[1]
	for (i in 1:l) {
		topline[i] <- baseline[i] + height[i]
		baseline[i + 1] <- topline[i]
		rightline[i] <- leftline[i] + width.vec[i]
		leftline[i + 1] <- rightline[i] + space.vec[i + 1]
	}
	ticks.vec = rep(0, times = l)

	##  Wait, do we need to turns this on its side?  Barplot uses an
	##  internal function that rearranges the order of arguments.  This is
	##  a bit less hassle to do it once, though the line drawing functions
	##  will need to be checked later.
	if (horiz == TRUE) {
		oldtopline <- topline
		topline <- leftline
		leftline <- oldtopline
		oldbaseline <- baseline
		baseline <- rightline
		rightline <- oldbaseline
		respaxis <- 1
		expaxis <- 2
		for (i in 1:l) ticks.vec[i] <- (baseline[i] + topline[i])/2
	}
	else {
		for (i in 1:l) ticks.vec[i] <- (leftline[i] + rightline[i])/2
	}
	if (is.null(xlim))
		xlim <- c(min(leftline), max(rightline))
	if (is.null(ylim))
		ylim <- c(min(topline, baseline), max(topline, baseline))
	if (plot == TRUE) {
		if (add == FALSE)
			plot.new()
		plot.window(xlim = xlim, ylim = ylim, ...)
		for (i in 1:l) rect(leftline[i], baseline[i], rightline[i],
                            topline[i], density = density.vec[i], angle = angle.vec[i],
                            col = col.vec[i], border = border.vec[i], xpd = xpd, ...)
		if (level.lines == TRUE)
			if (horiz == TRUE)
				for (i in 1:(l - 1))
                    lines(c(leftline[i], leftline[i]), c(baseline[i], topline[i + 1]), ...)
			else
                for (i in 1:(l - 1))
                    lines(c(rightline[i], leftline[i + 1]), c(topline[i], baseline[i + 1]), ...)
		if (!is.null(names.arg) && axisnames)
			axis(expaxis, at = ticks.vec, labels = names.arg,
                 lty = axis.lty, cex.axis = cex.names, ...)
		if (axes == TRUE)
			axis(respaxis, cex.axis = cex.axis, ...)
		title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...)
	}
	invisible(ticks.vec)
}
