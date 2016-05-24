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

#' @title Waterfall Chart
#'
#' @description
#' Creates a waterfall chart using Lattice
#'
#' @param x a formula describing the form of conditioning plot. The
#' formula is generally of the form 'y ~ x | g1 * g2 * ...', indicating
#' that plots of 'y' (on the y axis) versus 'x' (on the x axis) should
#' be produced conditional on the variables 'g1, g2, ...'. However, the
#' conditioning variables 'g1,g2,...' may be omitted.
#'
#' @param data a data frame containing values (or more precisely,
#' anything that is a valid 'envir' argument in 'eval', e.g., a list or
#' an environment) for any variables in the formula, as well as 'groups'
#' and 'subset' if applicable. If not found in 'data', or if 'data' is
#' unspecified, the variables are looked for in the environment of the
#' formula.
#'
#' @param groups a vector expected to act as a grouping variable within
#' each panel, typically used to distinguish different groups by varying
#' graphical parameters like color and line type.  Unlike with the
#' \link{barchart} function, groups specifies where subtotals columns,
#' should appear.  There is a subtotal created for each group specified.
#' If no groups are given, a summary column is still reported.
#'
#' @param horizontal This argument is used to process the arguments to
#' these high level functions, but more importantly, it is passed as an
#' argument to the panel function, which is supposed to use it as
#' appropriate.
#'
#' @param panel This draws the actual plot after \link{bwplot} has done
#' the difficult work of processing the formula.
#'
#' @param prepanel This function returns the \link{bwplot} information
#' on the number of columns to display and where to place labels.
#'
#' @param box.ratio specifies the ratio of the width of the rectangles
#' to the interrectangle space.
#'
#' @param origin initial offset relative to the x axis.  The value
#' serves as the logical starting point for the first column and any
#' summary column.  Defaults to 0.
#'
#' @param level.lines if FALSE, the lines connecting adjacent boxes are
#' ommitted from the display.
#'
#' @param summaryname name of the summary column, usually "Total"
#'
#' @param ... further arguments.
#'
#' @details
#'
#' This function closely mimics the \link{barchart} interface, but
#' provides a type of chart called a waterfall plot, showing how
#' multiple subvalues contribute to a total sum.
#'
#' The bulk of the work is actually processed in \link{bwplot} which
#' defines where tickmarks and other information outside the plot itself
#' are placed.  Only a formula method is provided.
#'
#' Matrix and vector interfaces are not provided because mimicing the
#' behavior of \link{barchart} for those interfaces produces
#' unintellible and undefined graphic output.
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
#' data(jaquith)
#' waterfallchart(value~label, data=rasiel, groups=rasiel$subtotal)
#' waterfallchart(factor~score, data=jaquith)
#'
#' @import lattice
#'
#' @export
waterfallchart <-
    function (x, data = NULL, groups = NULL, horizontal = FALSE,
              panel = panel.waterfallchart,
              prepanel = prepanel.waterfallchart, summaryname = "Total",
              box.ratio = 2, origin = 0, level.lines = TRUE, ...)
	UseMethod("waterfallchart")

#' @export
waterfallchart.formula <- function (x, data = NULL, groups = NULL, horizontal = FALSE,
                                    panel = panel.waterfallchart,
                                    prepanel = prepanel.waterfallchart, summaryname = "Total",
                                    box.ratio = 2, origin = 0, level.lines = TRUE, ...) {
    ocall <- sys.call(sys.parent())
    ocall[[1]] <- quote(waterfallchart)
    ccall <- match.call()
    ccall$data <- data
    ccall$panel <- panel
    ccall$prepanel <- prepanel
    ccall$box.ratio <- box.ratio
    if(is.null(groups))
        ccall$groups <- rep(summaryname, nrow(data))


    ## Let lattice do the hard work
    ccall[[1]] <- quote(bwplot)
    ans <- eval.parent(ccall)
    ans$call <- ocall
    ans
}

panel.waterfallchart <- function (x, y, box.ratio = 1, box.width = box.ratio/(1 + box.ratio),
                                  horizontal = FALSE, origin = 0, reference = TRUE, groups = NULL,
                                  summaryname = NULL, col = if (is.null(groups)) plot.polygon$col else superpose.polygon$col,
                                  border = if (is.null(groups)) plot.polygon$border else superpose.polygon$border,
                                  lty = if (is.null(groups)) plot.polygon$lty else superpose.polygon$lty,
                                  lwd = if (is.null(groups)) plot.polygon$lwd else superpose.polygon$lwd,
                                  level.lines = TRUE, ...) {
    plot.polygon <- trellis.par.get("plot.polygon")
	superpose.polygon <- trellis.par.get("superpose.polygon")
	reference.line <- trellis.par.get("reference.line")

	keep <- (function(x, y, groups, subscripts, ...) {
        !is.na(x) & !is.na(y)
    })(x = x, y = y, ...)

	if (!any(keep))
		return()
	x <- as.numeric(x[keep])
	y <- as.numeric(y[keep])
	grplst <- sort(unique(groups))
	baseline <- rep(origin, (l = length(x) + length(grplst)) + 1)

	##  The following block of code reorganizes the data into the final
	##  format.  The same block is used above.
	if (horizontal) {
		data <- merge(data.frame(x, y, groups, groupsy = y)[order(groups, y), ], data.frame(x = rep(NA, length(grplst)), y = rep(NA, length(grplst)), groups = grplst, groupsy = grplst), all = TRUE)
		data <- data[order(data$groups, data$y), ]
		if (is.null(origin)) {
			origin <- current.panel.limits()$xlim[1]
			reference <- FALSE
		}

		if (reference)
			panel.abline(h = origin, col = reference.line$col, lty = reference.line$lty, lwd = reference.line$lwd)
		for (i in 1:l) {
			if (is.na(data[i, ]$y)) {
				height <- origin - (baseline[i + 1] = baseline[i])
				color <- col[2]
			}
			else {
				baseline[i + 1] <- baseline[i] + (height = data[i, ]$x)
				color <- col[1]
			}
			panel.rect(y = i, x = baseline[i], col = color, lty = lty, lwd = lwd, height = box.width, width = height, just = c("left", "centre"))
		}
		if (level.lines == TRUE)
			for (i in 2:l - 1) panel.lines(y = c(i, i + 1), x = baseline[i +
                                                                             1], col = border, lty = lty, lwd = lwd)
	}
	else {
		data <- merge(data.frame(x, y, groups, groupsx = x)[order(groups, x), ], data.frame(x = rep(NA, length(grplst)), y = rep(NA, length(grplst)), groups = grplst, groupsx = grplst), all = TRUE)
		data <- data[order(data$groups, data$x), ]
		if (is.null(origin)) {
			origin <- current.panel.limits()$ylim[1]
			reference <- FALSE
		}

		if (reference)
			panel.abline(h = origin, col = reference.line$col, lty = reference.line$lty, lwd = reference.line$lwd)
		for (i in 1:l) {
			if (is.na(data[i, ]$y)) {
				height <- origin - (baseline[i + 1] = baseline[i])
				color <- col[2]
			}
			else {
				baseline[i + 1] <- baseline[i] + (height = data[i, ]$y)
				color <- col[1]
			}
			panel.rect(x = i, y = baseline[i], col = color, lty = lty, lwd = lwd, width = box.width, height = height, just = c("centre", "bottom"))
		}
		if (level.lines == TRUE)
			for (i in 2:l - 1)
                panel.lines(x = c(i, i + 1), y = baseline[i + 1], col = border, lty = lty, lwd = lwd)
	}
}

prepanel.waterfallchart <- function (x, y, horizontal = FALSE, origin = 0, groups = NULL, summaryname = NULL, ...) {

    if (any(!is.na(x) & !is.na(y))) {
        grplst <- sort(unique(groups))
        baseline <- rep(origin, (l = length(x) + length(grplst)))

        ##  The following block of code reorganizes the data into the final
        ##  format.  The same block is used above.
        if (horizontal) {
            data <- merge(data.frame(x, y, groups, groupsy = y)[order(groups, y), ], data.frame(x = rep(NA, length(grplst)), y = rep(NA, length(grplst)), groups = grplst, groupsy = grplst), all = TRUE)
            data <- data[order(data$groups, data$y), ]
            list(ylim = levels(factor(data$groupsy, levels = data$groupsy)),
                 yat = sort(unique(as.numeric(data$groupsy, levels = data$groupsy))),
                 xlim = {
                     for (i in 1:(l - 1)) {
                         if (!is.na(data[i, ]$x)) baseline[i + 1] <- baseline[i] +
                             data[i, ]$x else baseline[i + 1] <- baseline[i]
                     }
                     range(baseline)
                 }, dx = 1, dy = 1)
        }
        else {
            data <- merge(data.frame(x, y, groups, groupsx = x)[order(groups, x), ], data.frame(x = rep(NA, length(grplst)), y = rep(NA, length(grplst)), groups = grplst, groupsx = grplst), all = TRUE)
            data <- data[order(data$groups, data$x), ]
            list(xlim = levels(factor(data$groupsx, levels = data$groupsx)),
                 xat = sort(unique(as.numeric(data$groupsx, levels = data$groupsx))),
                 ylim = {
                     for (i in 1:(l - 1)) {
                         if (!is.na(data[i, ]$y)) baseline[i + 1] <- baseline[i] +
                             data[i, ]$y else baseline[i + 1] <- baseline[i]
                     }
                     range(baseline)
                 }, dx = 1, dy = 1)
        }
    }
    else list(xlim = c(NA, NA), ylim = c(NA, NA), dx = 1, dy = 1)
}
