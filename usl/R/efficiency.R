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
#' Efficiency of the system
#'
#' The efficiency of a system expressed in terms of the deviation from
#' linear scalability.
#'
#' The function returns a vector which contains the deviation from linearity
#' for every measurement of the model input.  A value of \code{1} indicates
#' linear scalability while values less than \code{1} correspond to the
#' fraction of the measurement compared to linear scalability.
#'
#' @param object A USL object.
#'
#' @return A vector of numeric values.
#'
#' @seealso \code{\link{usl}}
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
#' ## Show the efficiency
#' efficiency(usl(throughput ~ processors, raytracer))
#'
#' @aliases efficiency
#' @export
#'
setMethod(
  f = "efficiency",
  signature = "USL",
  definition = function(object) {
    return(object@efficiency)
  }
)
