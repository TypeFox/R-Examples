## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015 - 2016  Stefan Engblom
## Copyright (C) 2015 - 2016  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

##' Run \code{siminf_model} on scaled parameters
##'
##' @rdname run_outer-methods
##' @docType methods
##' @param x Scale the model \code{gdata} parameter values on the
##'     right hand side of the formula with \code{x} before calling
##'     \code{FUN} with the scaled model as argument.
##' @param y Scale the model \code{gdata} parameter values on the
##'     left hand side of the formula with \code{y} before calling
##'     \code{FUN} with the scaled model as argument.
##' @param model The siminf model to scale parameters on and run.
##' @param formula The parameters in the \code{gdata} vector matching
##'     the left hand side of the formula \code{a + b ~ c} will be
##'     scaled by \code{y}.  The parameters in the \code{gdata} vector
##'     matching the right hand side of the formula \code{a + b ~ c}
##'     will be scaled by \code{x}.
##' @param FUN A function to use on the scaled model 'gdata' parameters.
##' @param ... Optional arguments to be passed to \code{FUN}.
##' @return Array with dimension \code{c(dim(x), dim(y))}.
##' @include siminf_model.R
##' @importFrom stats terms
##' @examples
##' \dontrun{
##' ## Define a function that runs the model and
##' ## returns the prevalence at the last day.
##' run_f <- function(model) {prevalence(run(model))[length(model@@tspan)]}
##'
##' ## Define scaling parameters
##' x = seq(from = 0.95, to = 1.05, by = 0.01)
##' y = seq(from = 0.95, to = 1.05, by = 0.01)
##'
##' ## Run model
##' model <- demo_model(nodes = 10)
##' formula <- upsilon ~ beta_t1 + beta_t2 + beta_t3 + beta_t4
##' z <- run_outer(x, y, model, formula, run_f)
##'
##' ## Plot
##' filled.contour(x, y, z)
##' }
setGeneric("run_outer",
           signature = c("x", "y", "model"),
           function(x,
                    y,
                    model,
                    formula = NULL,
                    FUN     = NULL,
                    ...) standardGeneric("run_outer"))

##' @rdname run_outer-methods
##' @export
setMethod("run_outer",
          signature(x = "numeric", y = "numeric", model = "siminf_model"),
          function(x, y, model, formula, FUN, ...)
          {
              if (is.null(names(model@gdata)))
                  stop("'names(model@gdata)' is NULL")
              if (is.null(formula))
                  stop("'formula' argument is NULL")
              if (is.null(FUN))
                  stop("'FUN' argument is NULL")
              FUN <- match.fun(FUN)

              ## Determine indices to the 'gdata' parameters to scale
              ## by 'x'
              xx <- attr(stats::terms(formula, allowDotAsName = TRUE), "term.labels")
              xx <- xx[attr(stats::terms(formula, allowDotAsName = TRUE), "order") == 1]
              if (length(xx) < 1)
                  stop("Invalid parameters on the right side of the formula")
              x_i <- match(xx, names(model@gdata))
              if (any(is.na(x_i)))
                  stop("Unmatched parameters on the right side of the formula")

              ## Determine indices to the 'gdata' parameters to scale
              ## by 'y'
              yy <- attr(stats::terms(formula, allowDotAsName = TRUE), "response")
              if (yy < 1)
                  stop("Invalid parameters on the left side of the formula")
              vars <- attr(stats::terms(formula, allowDotAsName = TRUE), "variables")[-1]
              yy <- as.character(vars[yy])
              yy <- unlist(strsplit(yy, "+", fixed = TRUE))
              yy <- sub("^\\s", "", sub("\\s$", "", yy))
              y_i <- match(yy, names(model@gdata))
              if (any(is.na(y_i)))
                  stop("Unmatched parameters on the left hand side of the formula")

              outer(x, y, function(x, y, ...) {
                  run_internal <- function(x, y, x_i, y_i, model, ...) {
                      model@gdata[x_i] <- model@gdata[x_i] * x
                      model@gdata[y_i] <- model@gdata[y_i] * y
                      FUN(model, ...)
                  }

                  sapply(seq_len(length(x)), function(i, ...) {
                      run_internal(x[i], y[i], x_i, y_i, model, ...)
                  }, ...)
              }, ...)
          }
)
