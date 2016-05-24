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
#' Confidence Intervals for USL model parameters
#'
#' Estimate confidence intervals for one or more parameters in a USL model.
#' The intervals are calculated from the parameter standard error using the
#' Student t distribution at the given level.
#' 
#' Bootstrapping is no longer used to estimate confidence intervals.
#'
#' @param object A USL object.
#' @param parm A specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing,
#'   all parameters are considered.
#' @param level The confidence level required.
#' @param type This parameter is no longer used and will be removed in the next
#'   version.
#'
#' @return A matrix (or vector) with columns giving lower and upper confidence
#'   limits for each parameter. These will be labelled as (1-level)/2 and
#'   1 - (1-level)/2 in \% (by default 2.5\% and 97.5\%).
#'
#' @seealso \code{\link{usl}}
#'
#' @examples
#' require(usl)
#'
#' data(specsdm91)
#'
#' ## Create USL model
#' usl.model <- usl(throughput ~ load, specsdm91)
#'
#' ## Print confidence intervals
#' confint(usl.model)
#'
#' @export
#'
setMethod(
  f = "confint",
  signature = "USL",
  definition = function(object, parm, level = 0.95, type) {
    ci.value <- NULL # vector with confidence interval values
    
    # Degree of freedom for Student t distribution
    df <- length(object@residuals) - 1L

    # Vectors to collect column and row names of result matrix
    col.name <- paste(formatC(100 * c((1-level)/2, 1-(1-level)/2)), "%")
    row.name <- NULL
    
    # Warn about old parameter usage
    if (!missing(type)) {
      warning("parameter 'type' is no longer used")
    }

    # Return confidence intervals for both parameters if 'parm' is unset
    if (missing(parm)) parm <- object@coef.names

    # Replace numeric parameters with named parameters
    if (mode(parm) == "numeric") {
      parm <- as.character(parm)

      parm <- gsub("1", "sigma", parm, ignore.case = TRUE)
      parm <- gsub("2", "kappa", parm, ignore.case = TRUE)
    }

    # Calculate confidence intervals for the given level
    for (i in object@coef.names) {
      if (i %in% parm) {
        pa <- object@coefficients[i]
        se <- object@coef.std.err[i] * qt(level, df)

        ci.value <- c(ci.value, pa - se, pa + se)
        row.name <- c(row.name, i)
      }
    }

    # Build dummy matrix if no sensible parameters were requested
    if (length(row.name) < 1) {
      row.name <- NA
      ci.value <- c(NA, NA)
    }

    # Return confidence intervals as matrix
    matrix(ci.value, nrow = length(row.name), ncol = 2,
           byrow = TRUE, dimnames = list(row.name, col.name))
  })
