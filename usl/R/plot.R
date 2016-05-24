# Copyright (c) 2013, 2014 Stefan Moeding
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.


##############################################################################
#' Plot the scalability function from a USL model
#'
#' Create a line plot for the scalability functionh of a Universal
#' Scalability Law model.
#'
#' \code{plot} creates a plot of the scalability function for the model
#' represented by the argument \code{x}.
#'
#' If \code{from} is not specified then the range starts at the minimum value
#' given to define the model. An unspecified value for \code{to} will lead
#' to plot ending at the maximum value from the model. For \code{add = TRUE}
#' the defaults are taken from the limits of the previous plot.
#'
#' \code{xlab} and \code{ylab} can be used to set the axis titles. The defaults
#' are the names of the regressor and response variables used in the model.
#'
#' If the parameter \code{bounds} is set to \code{TRUE} then the plot also
#' shows dotted lines for the theoretical bounds of scalability. These are
#' the linear scalability for small loads and Amdahl's asymptote for the
#' limit of scalability as load approaches infinity. 
#' 
#' The parameters \code{sigma} or \code{kappa} are useful to do a what-if
#' analysis. Setting these parameters override the model parameters and show
#' how the system would behave with a different contention or coherency delay
#' parameter.
#'
#' @param x The USL object to plot.
#' @param from The start of the range over which the scalability function
#'   will be plotted.
#' @param to The end of the range over which the scalability function
#'   will be plotted.
#' @param xlab A title for the x axis: see \code{\link{title}}.
#' @param ylab A title for the y axis: see \code{\link{title}}.
#' @param bounds Add the bounds of scalability to the plot. 
#' @param sigma Optional parameter to be used for evaluation instead of the
#'   parameter computed for the model.
#' @param kappa Optional parameter to be used for evaluation instead of the
#'   parameter computed for the model.
#' @param ... Other graphical parameters passed to plot
#'   (see \code{\link{par}}, \code{\link{plot.function}}).
#'
#' @seealso \code{\link{usl}}, \code{\link{plot.function}}
#'
#' @examples
#' require(usl)
#'
#' data(raytracer)
#'
#' ## Plot result from USL model for demo dataset
#' plot(usl(throughput ~ processors, raytracer))
#'
#' @export
#'
setMethod(
  f = "plot",
  signature = "USL",
  definition = function(x, from = NULL, to = NULL, xlab = NULL, ylab = NULL,
                        bounds = FALSE, sigma, kappa, ...) {
    # Take range from the model if not specified
    if (missing(from)) from <- min(x@frame[, x@regr])
    if (missing(to)) to <- max(x@frame[, x@regr])

    # Set titles for axis
    if (missing(xlab)) xlab <- x@regr
    if (missing(ylab)) ylab <- x@resp

    # Use explicitly specified coefficients
    if (missing(sigma)) sigma <- coef(x)[['sigma']]
    if (missing(kappa)) kappa <- coef(x)[['kappa']]

    # Get the function to calculate scalability for the model
    .func <- scalability(x, sigma, kappa)

    # Plot the scalability function
    plot(x = .func, from = from, to = to, xlab = xlab, ylab = ylab, ...)
    
    # Add theoretical bounds of scalability to the plot
    if (bounds) {
      # Bound 1: linear scalability
      abline(a = 0, b = x@scale.factor, lty = "dotted")
      
      # Bound 2: Amdahl's asymptote
      if (sigma > 0) {
        abline(h = 1/sigma * x@scale.factor, lty = "dotted")
      }
    }
  }
)
