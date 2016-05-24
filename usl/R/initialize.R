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
#' Generate an object from the "\code{USL}" class
#'
#' Initialize the object.
#'
#' @param .Object The object to initialize.
#' @param call The formula used to create the USL model.
#' @param frame The model frame containing the variables in the model.
#' @param regr The name of the regressor variable in the model.
#' @param resp The name of the response variable in the model.
#' @param scale.factor A numeric value for the scale of the model. This is the
#'     factor by which the model values have been reduced to get a normalized
#'     model.
#' @param sigma The contention parameter of the model.
#' @param kappa The coherency delay parameter of the model.
#'
#' @return An object of the specific type.
#'
#' @keywords internal
#'
setMethod(
  f = "initialize",
  signature = "USL",
  definition = function(.Object, call, frame, regr, resp, scale.factor, sigma, kappa) {
    .Object@call         <- call
    .Object@coefficients <- structure(c(sigma, kappa), names = .Object@coef.names)
    .Object@frame        <- frame
    .Object@regr         <- regr
    .Object@resp         <- resp
    .Object@scale.factor <- scale.factor
    .Object@efficiency   <- structure(frame[[resp]] / scale.factor / frame[[regr]],
                                      names = frame[, regr])
    .Object@df.residual  <- length(frame[[resp]]) - 2L

    # Call inspector
    validObject(.Object)

    return(.Object)
  }
)
