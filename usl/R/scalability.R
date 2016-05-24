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
#' Scalability function of a USL model
#'
#' \code{scalability} is a higher order function and returns a function to
#' calculate the scalability for the specific USL model.
#'
#' The returned function can be used to calculate specific values once the
#' model for a system has been created.
#'
#' The parameters \code{sigma} or \code{kappa} are useful to do a what-if
#' analysis. Setting these parameters override the model parameters and show
#' how the system would behave with a different contention or coherency delay
#' parameter.
#'
#' @param object A USL object.
#' @param sigma Optional parameter to be used for evaluation instead of the
#'   parameter computed for the model.
#' @param kappa Optional parameter to be used for evaluation instead of the
#'   parameter computed for the model.
#'
#' @return A function with parameter \code{x} that calculates the
#'   scalability value of the specific model.
#'
#' @seealso \code{\link{usl}}, \code{\link{peak.scalability,USL-method}}
#'
#' @references Neil J. Gunther. Guerrilla Capacity Planning: A Tactical
#'   Approach to Planning for Highly Scalable Applications and Services.
#'   Springer, Heidelberg, Germany, 1st edition, 2007.
#'
#' @examples
#' require(usl)
#'
#' data(raytracer)
#'
#' ## Compute the scalability function
#' s <- scalability(usl(throughput ~ processors, raytracer))
#'
#' ## Print scalability for 32 CPUs for the demo dataset
#' print(s(32))
#'
#' ## Plot scalability for the range from 1 to 64 CPUs
#' plot(s, from=1, to=64)
#'
#' @aliases scalability
#' @export
#'
setMethod(
  f = "scalability",
  signature = "USL",
  definition = function(object, sigma, kappa) {
    if (missing(sigma)) sigma <- coef(object)[['sigma']]
    if (missing(kappa)) kappa <- coef(object)[['kappa']]

    .func <- function(x) {
      # Formula (4.31) on page 57 of GCaP:
      cap <- x / (1 + (sigma * (x-1)) + (kappa * x * (x-1)))

      # Scale it to the measurements
      return(object@scale.factor * cap)
    }

    # Return the usl function (lexically scoped)
    return(.func)
  }
)


##############################################################################
#' Peak scalability value of a USL model
#'
#' Calculate the point of peak scalability for a specific model.
#'
#' The peak scalability is the point where the throughput of the
#' system starts to go retrograde, i.e., starts to decrease with
#' increasing load.
#'
#' The parameters \code{sigma} or \code{kappa} are useful to do a what-if
#' analysis. Setting these parameters override the model parameters and show
#' how the system would behave with a different contention or coherency delay
#' parameter.
#'
#' See formula (4.33) in \emph{Guerilla Capacity Planning}.
#'
#' @param object A USL object.
#' @param sigma Optional parameter to be used for evaluation instead of the
#'   parameter computed for the model.
#' @param kappa Optional parameter to be used for evaluation instead of the
#'   parameter computed for the model.
#'
#' @return A numeric value for the point where peak scalability will be
#'   reached.
#'
#' @seealso \code{\link{usl}}, \code{\link{scalability,USL-method}}
#'
#' @references Neil J. Gunther. Guerrilla Capacity Planning: A Tactical
#'   Approach to Planning for Highly Scalable Applications and Services.
#'   Springer, Heidelberg, Germany, 1st edition, 2007.
#'
#' @examples
#' require(usl)
#'
#' data(raytracer)
#'
#' peak.scalability(usl(throughput ~ processors, raytracer))
#' ## Peak scalability will be reached just below 450 processors
#'
#' @aliases peak.scalability
#' @export
#'
setMethod(
  f = "peak.scalability",
  signature = "USL",
  definition = function(object, sigma, kappa) {
    if (missing(sigma)) sigma <- coef(object)[['sigma']]
    if (missing(kappa)) kappa <- coef(object)[['kappa']]

    return(sqrt((1 - sigma) / kappa))
  }
)
