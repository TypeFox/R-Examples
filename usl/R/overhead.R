# Copyright (c) 2014 Stefan Moeding
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
#' Overhead method for Universal Scalability Law models
#'
#' \code{overhead} calculates the overhead in processing time for a system
#' modeled with the Universal Scalability Law.
#' It evaluates the regression function in the frame \code{newdata} (which
#' defaults to \code{model.frame(object)}). The result contains the ideal
#' processing time and the additional overhead caused by contention and
#' coherency delays.
#'
#' The calculated processing times are given as percentages of a
#' non-parallelized workload. So for a non-parallelized workload the ideal
#' processing time will always be given as \emph{100\%} while the overhead
#' for contention and coherency will always be zero.
#'
#' Doubling the capacity will cut the ideal processing time in half but
#' increase the overhead percentages. The increase of the overhead depends on
#' the values of the parameters \code{sigma} and \code{kappa} estimated by
#' \code{\link{usl}}.
#'
#' The calculation is based on \emph{A General Theory of Computational
#' Scalability Based on Rational Functions}, equation 26.
#'
#' @param object A USL model object for which the overhead will be calculated.
#' @param newdata An optional data frame in which to look for variables
#'   with which to calculate the overhead.
#'   If omitted, the fitted values are used.
#'
#' @return \code{overhead} produces a matrix of overhead percentages based on
#'   a non-parallelized workload. The column \code{ideal} contains the ideal
#'   percentage of execution time. The columns \code{contention} and
#'   \code{coherency} give the additional overhead percentage caused by
#'   the respective effects.
#'
#' @seealso \code{\link{usl}}, \code{\link{USL-class}}
#'
#' @references Neil J. Gunther. Guerrilla Capacity Planning: A Tactical
#'   Approach to Planning for Highly Scalable Applications and Services.
#'   Springer, Heidelberg, Germany, 1st edition, 2007.
#'
#' @references Neil J. Gunther. A General Theory of Computational Scalability
#'   Based on Rational Functions. Computing Research Repository, 2008.
#'   \code{http://arxiv.org/abs/0808.1431}
#'
#' @examples
#' require(usl)
#'
#' data(raytracer)
#'
#' ## Print overhead in processing time for demo dataset
#' overhead(usl(throughput ~ processors, raytracer))
#'
#' @aliases overhead
#' @export
#'
setMethod(
  f = "overhead",
  signature = "USL",
  definition = function(object, newdata) {
    # Calculate overhead for the initial data used to create
    # the model if no data frame 'newdata' is given as parameter
    if (missing(newdata)) newdata <- object@frame

    # Extract regressor variable from data frame
    x <- newdata[, object@regr, drop=TRUE]

    y.ideal      <- 1 / x
    y.contention <- coef(object)[['sigma']] * (x - 1) / x
    y.coherency  <- coef(object)[['kappa']] * (1/2) * (x - 1)

    col.names <-  c("ideal", "contention", "coherency")

    # Return the matrix
    matrix(c(y.ideal, y.contention, y.coherency),
           nrow = length(x), dimnames = list(seq(x), col.names))
  }
)
