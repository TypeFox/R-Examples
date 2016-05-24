# Copyright (c) 2013-2016 Stefan Moeding
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
#' Solve a USL model using a transformation to a 2nd degree polynom
#'
#' This function solves a USL model using the transformation introduced in
#' sections 5.5.1 - 5.5.3 of \emph{Guerrilla Capacity Planning}.
#'
#' @param model A data frame with two columns containing the values of the
#'   predictor variable in the first column and the values of the response
#'   variable in the second column.
#'
#' @return A list containing three elements: the scale.factor of the model,
#'   the model coefficients sigma and kappa.
#'
#' @seealso \code{\link{usl}}
#'
#' @references Neil J. Gunther. Guerrilla Capacity Planning: A Tactical
#'   Approach to Planning for Highly Scalable Applications and Services.
#'   Springer, Heidelberg, Germany, 1st edition, 2007.
#'
#' @keywords internal
#'
usl.solve.lm <- function(model) {
  # Verify that the scale factor for normalization is in the dataframe
  if (all(model[ , 1] != 1)) {
    stop(paste0("'data' must contain a row where '", names(model[1]), "' = 1"))
  }

  # Calculate scale factor: get throughput for entry where load=1
  scale.factor <- model[match(1, model[ , 1]), 2]

  # Rename columns
  names(model) <- c("load", "throughput")

  # normalize data (cf. GCaP chapter 5.4)
  model$capacity <- model$throughput / scale.factor

  # compute deviations from linearity (cf. GCaP chapter 5.5.2)
  model$x <- model$load - 1
  model$y <- (model$load / model$capacity) - 1

  # Solve quadratic model without intercept
  model.fit <- lm(y ~ I(x^2) + x - 1, data = model)

  # Calculate coefficients sigma & kappa used by the USL model
  sigma <- coef(model.fit)[[2]] - coef(model.fit)[[1]]
  kappa <- coef(model.fit)[[1]]

  return(list(scale.factor = scale.factor, sigma = sigma, kappa = kappa))
}


##############################################################################
#' Solve a USL model using non linear regression
#'
#' This function solves a USL model using non linear regression with least
#' squares. It uses the function \code{\link{nls}} with the "\code{port}"
#' algorithm to perform the calculation. All restrictions of the algorithm
#' apply.
#'
#' @param model A data frame with two columns containing the values of the
#'   predictor variable in the first column and the values of the response
#'   variable in the second column.
#'
#' @return A list containing three elements: the scale.factor of the model,
#'   the model coefficients sigma and kappa.
#'
#' @seealso \code{\link{usl}}
#' @keywords internal
#'
usl.solve.nls <- function(model) {
  names(model) <- c("x", "y")

  # Lower bound for scale.factor?
  sf.max <- max(model$y / model$x)

  model.fit <- nls(y ~ X1 * x/(1 + sigma * (x-1) + kappa * x * (x-1)),
                   data = model,
                   start = c(X1 = sf.max, sigma = 0.1, kappa = 0.01),
                   algorithm = "port",
                   lower = c(X1 = 0, sigma = 0, kappa = 0),
                   upper = c(X1 = Inf, sigma = 1, kappa = 1))

  scale.factor = coef(model.fit)[['X1']]
  sigma = coef(model.fit)[['sigma']]
  kappa = coef(model.fit)[['kappa']]

  return(list(scale.factor = scale.factor, sigma = sigma, kappa = kappa))
}


##############################################################################
#' Solve a USL model using non linear regression
#'
#' This function solves a USL model using non linear regression with least
#' squares. It uses the function \code{\link{nlxb}} from the \pkg{nlmrt}
#' package to perform the calculation.
#'
#' @param model A data frame with two columns containing the values of the
#'   predictor variable in the first column and the values of the response
#'   variable in the second column.
#'
#' @return A list containing three elements: the scale.factor of the model,
#'   the model coefficients sigma and kappa.
#'
#' @seealso \code{\link{usl}}
#'
#' @references John C. Nash. nlmrt: Functions for nonlinear least squares
#'   solutions, 2013. R package version 2013-8.10.
#'
#' @importFrom nlmrt nlxb
#' @keywords internal
#'
usl.solve.nlxb <- function(model) {
  names(model) <- c("x", "y")

  # Lower bound for scale.factor?
  sf.max <- max(model$y / model$x)

  model.fit <- nlxb(y ~ X1 * x/(1 + sigma * (x-1) + kappa * x * (x-1)),
                    data = model,
                    start = c(X1 = sf.max, sigma = 0.1, kappa = 0.01),
                    lower = c(X1 = 0, sigma = 0, kappa = 0),
                    upper = c(X1 = Inf, sigma = 1, kappa = 1))

  scale.factor = model.fit$coefficients[['X1']]
  sigma = model.fit$coefficients[['sigma']]
  kappa = model.fit$coefficients[['kappa']]

  return(list(scale.factor = scale.factor, sigma = sigma, kappa = kappa))
}


##############################################################################
#' Create a model for the Universal Scalability Law
#'
#' \code{usl} is used to create a model for the Universal Scalability Law.
#'
#' The Universal Scalability Law is used to forcast the scalability of
#' either a hardware or a software system.
#'
#' The USL model works with one independent variable (e.g. virtual users,
#' processes, threads, ...) and one dependent variable (e.g. throughput, ...).
#' Therefore the model formula must be in the simple
#' "\code{response ~ predictor}" format.
#'
#' The model produces two coefficients as result: \code{sigma} models the
#' contention and \code{kappa} the coherency delay of the system. The
#' function \code{\link{coef}} extracts the coefficients from the model
#' object.
#'
#' The argument \code{method} selects which solver is used to solve the
#' model:
#'
#' \itemize{
#'   \item "\code{default}" for the default method using a transformation
#'     into a 2nd degree polynom. It can only be used if the model frame
#'     contains a value for the normalization where the predictor equals
#'     "\code{1}" for one measurement. This is the algorithm introduced by
#'     Dr. Neil J. Gunther in the book \emph{Guerrilla Capacity Planning}.
#'   \item "\code{nls}" for a nonlinear regression model. This method
#'     estimates not only the coefficients \code{sigma} and \code{kappa} but
#'     also the \code{scale.factor} for the normalization. \code{\link{nls}}
#'     with the "\code{port}" algorithm is used internally to solve the
#'     model. So all restrictions of the "\code{port}" algorithm apply.
#'   \item "\code{nlxb}" for a nonliner regression model using the function
#'     \code{\link{nlxb}} from the \code{\link{nlmrt}} package. This method
#'     also estimates both coefficients and the normalization factor. It is
#'     expected to be more robust than the \code{nls} method.
#' }
#'
#' The "\code{nlxb}" solver is used as fallback if the "\code{default}"
#' method is selected and a predictor equal "\code{1}" is missing. A warning
#' message will be printed in this case.
#'
#' The Universal Scalability Law can be expressed with following formula.
#' \code{C(N)} predicts the relative capacity of the system for a given
#' load \code{N}:
#'
#' \deqn{C(N) = \frac{N}{1 + \sigma (N - 1) + \kappa N (N - 1)}}{C(N) = N / (1 + \sigma * (N - 1) + \kappa * N * (N - 1))}
#'
#' @param formula An object of class "\code{\link{formula}}" (or one that
#'   can be coerced to that class): a symbolic description of the model to be
#'   analyzed. The details of model specification are given under 'Details'.
#' @param data A data frame, list or environment (or object coercible by
#'   as.data.frame to a data frame) containing the variables in the model.
#'   If not found in data, the variables are taken from
#'   \code{environment(formula)}, typically the environment from which
#'   \code{usl} is called.
#' @param method Character value specifying the method to use. The possible
#'   values are described under 'Details'.
#'
#' @return An object of class USL.
#'
#' @seealso \code{\link{efficiency,USL-method}},
#'   \code{\link{scalability,USL-method}},
#'   \code{\link{peak.scalability,USL-method}},
#'   \code{\link{summary,USL-method}},
#'   \code{\link{predict,USL-method}},
#'   \code{\link{overhead,USL-method}},
#'   \code{\link{confint,USL-method}},
#'   \code{\link{coef}},
#'   \code{\link{fitted}},
#'   \code{\link{residuals}},
#'   \code{\link{df.residual}}
#'
#' @references Neil J. Gunther. Guerrilla Capacity Planning: A Tactical
#'   Approach to Planning for Highly Scalable Applications and Services.
#'   Springer, Heidelberg, Germany, 1st edition, 2007.
#'
#' @references John C. Nash. nlmrt: Functions for nonlinear least squares
#'   solutions, 2013. R package version 2013-8.10.
#'
#' @examples
#' require(usl)
#'
#' data(raytracer)
#'
#' ## Create USL model for "throughput" by "processors"
#' usl.model <- usl(throughput ~ processors, raytracer)
#'
#' ## Show summary of model parameters
#' summary(usl.model)
#'
#' ## Show complete list of efficiency parameters
#' efficiency(usl.model)
#'
#' ## Extract coefficients for model
#' coef(usl.model)
#'
#' ## Calculate point of peak scalability
#' peak.scalability(usl.model)
#'
#' ## Plot original data and scalability function
#' plot(raytracer)
#' plot(usl.model, add=TRUE)
#'
#' @export
#'
usl <- function(formula, data, method = "default") {
  ## canonicalize the arguments
  formula <- as.formula(formula)

  if (length(formula) < 3L) {
    stop("'formula' must be a 3-part formula")
  }

  if(!is.data.frame(data) && !is.environment(data)) {
    stop("'data' must be a data frame or an environment")
  }

  # Check parameter and variable names from formula
  var.names <- all.vars(formula)

  if (length(var.names) != 2L) {
    stop("'formula' must contain exactly 2 variables")
  }

  # Create model frame
  call <- match.call()
  frame <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(frame), 0)
  frame <- frame[c(1, m)]
  frame$na.action <- "na.omit"
  frame$drop.unused.levels <- TRUE
  frame[[1]] <- as.name("model.frame")
  frame <- eval(frame, parent.frame())

  # Verify there are enough values to do the calculation
  if (nrow(frame) < 6) {
    warning("'data' has only a few values; the result might not be accurate")
  }

  # Extract terms from the formula and get the names of the
  # predictor and response variables given by the user
  mt <- attr(frame, "terms")

  regr <- var.names[-attr(mt, "response")] # predictor
  resp <- var.names[attr(mt, "response")]  # response

  model.input <- data.frame(frame[regr], frame[resp])

  # Choose solver function
  sel <- switch(method, nls=2, nlxb=3, 1)
  usl.solve <- switch(sel, usl.solve.lm, usl.solve.nls, usl.solve.nlxb)
  
  # Use method 'nlxb' as fallback if scale factor is missing in data
  if ((sel == 1)  && (all(frame[regr] != 1))) {
    usl.solve <- usl.solve.nlxb
    warning(paste0("'data' has no row where '", regr, "' = 1; ",
                   "switching method from 'default' to 'nlxb'"))
  }

  # Solve the model for the model frame
  model.result <- usl.solve(model.input)

  # Create object for class USL
  .Object <- new(Class = "USL", call, frame, regr, resp,
                 model.result[['scale.factor']],
                 model.result[['sigma']], model.result[['kappa']])

  # Finish building the USL object
  nam <- row.names(frame)

  y.obs <- frame[, resp, drop = TRUE]
  y.fit <- predict(.Object)
  y.res <- y.obs - y.fit

  .Object@fitted    <- structure(y.fit, names = nam)
  .Object@residuals <- structure(y.res, names = nam)

  n <- length(y.obs) # sample size
  p <- 1             # number of regressors

  .Object@r.squared     <- 1 - (sum(y.res ^ 2) / sum((y.obs - mean(y.obs)) ^ 2))
  .Object@adj.r.squared <- 1 - (1 - .Object@r.squared) * ((n-1) / (n-p-1))


  # The following estimation of the standard errors is based on the
  # source code of the nls() function in R base.
  # See also: Nonlinear Regression and Nonlinear Least Squares,
  # Appendix to An R and S-PLUS Companion to Applied Regression, John
  # Fox, January 2002

  # residual variance
  df <- df.residual(.Object)
  resvar <- if(df <= 0) NaN else sum(y.res ^ 2) / df

  # Build gradient matrix
  grad <- gradient.usl(.Object)

  XtXinv <- solve(t(grad) %*% grad)

  # Standard error of coefficients
  .Object@coef.std.err <- sqrt(diag(XtXinv) * resvar)

  return(.Object)
}
