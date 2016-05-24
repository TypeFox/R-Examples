#' Calculate error between actual and forecast.
#'
#' \code{e} takes actual and forecast numeric vectors and returns a
#' numeric vector where forecast is subtracted from the actual.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @param actual a numeric vector of actuals.
#' @param forecast a numeric vector of forecasts.
#' @return  a numeric vector of errors.
#' @examples
#' # Examples
#' e(1:100, 100:1)
#' @export
e <- function(y_true = NULL, y_pred = NULL, actual = NULL, forecast = NULL) {
  warning("'e' is deprecated. Use 'absolute_error' instead.")
  if (is.null(actual)) {
    actual <- y_true
  } else {
    warning("Argument actual is deprecated; please use y_true instead.", call. = FALSE)
  }
  if (is.null(forecast)) {
    forecast <- y_pred
  } else {
    warning("Argument forecast is deprecated; please use y_pred instead.", call. = FALSE)
  }
  return(actual - forecast)
}
#' Calculate absolute error between actual and forecast.
#'
#' \code{ae} takes actual and forecast numeric vectors and returns a
#' numeric vector where forecast is subtracted from the actual and
#' then the absolute value is taken of those errors.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @param actual a numeric vector of actuals.
#' @param forecast a numeric vector of forecasts.
#' @return  a numeric vector of absolute errors.
#' @examples
#' # Examples
#' ae(1:100, 100:1)
#' @export
ae <- function(y_true = NULL, y_pred = NULL, actual = NULL, forecast = NULL) {
  warning("'ae' is deprecated. Use 'absolute_error' instead.")
  if (is.null(actual)) {
    actual <- y_true
  } else {
    warning("Argument actual is deprecated; please use y_true instead.", call. = FALSE)
  }
  if (is.null(forecast)) {
    forecast <- y_pred
  } else {
    warning("Argument forecast is deprecated; please use y_pred instead.", call. = FALSE)
  }
  return(abs(e(actual, forecast)))
}

#' Calculate percent error between actual and forecast.
#'
#' \code{pe} takes actual and forecast numeric vectors and returns a
#' numeric vector where forecast is subtracted from the actual and
#' then those errors are divided by the actuals.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @param actual a numeric vector of actuals.
#' @param forecast a numeric vector of forecasts.
#' @return  a numeric vector of percent errors.
#' @examples
#' # Examples
#' pe(1:100, 100:1)
#' @export
pe <- function(y_true = NULL, y_pred = NULL, actual = NULL, forecast = NULL) {
  warning("'pe' is deprecated. Use 'percent_error' instead.")
  if (is.null(actual)) {
    actual <- y_true
  } else {
    warning("Argument actual is deprecated; please use y_true instead.", call. = FALSE)
  }
  if (is.null(forecast)) {
    forecast <- y_pred
  } else {
    warning("Argument forecast is deprecated; please use y_pred instead.", call. = FALSE)
  }
  return(e(actual, forecast) / actual)
}

#' Calculate absolute percent error between actual and forecast.
#'
#' \code{pe} takes actual and forecast numeric vectors and returns a
#' numeric vector where forecast is subtracted from the actual and
#' then those errors are divided by the actuals. Lastly, the
#' absolute value of those percent errors are taken.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @param actual a numeric vector of actuals.
#' @param forecast a numeric vector of forecasts.
#' @return  a numeric vector of absolute percent errors.
#' @examples
#' # Examples
#' ape(1:100, 100:1)
#' @export
ape <- function(y_true = NULL, y_pred = NULL, actual = NULL, forecast = NULL) {
  warning("'ape' is deprecated. Use 'absolute_percent_error' instead.")
  if (is.null(actual)) {
    actual <- y_true
  } else {
    warning("Argument actual is deprecated; please use y_true instead.", call. = FALSE)
  }
  if (is.null(forecast)) {
    forecast <- y_pred
  } else {
    warning("Argument forecast is deprecated; please use y_pred instead.", call. = FALSE)
  }
  return(abs(pe(actual, forecast)))
}

#' Calculate mean absolute percent error between actual and forecast.
#'
#' \code{mape} takes actual and forecast numeric vectors and returns a
#' numeric vector where forecast is subtracted from the actual and
#' then those errors are divided by the actuals, the
#' absolute value of those percent errors are then taken, and lastly,
#' the mean of those absolute percent errors are taken.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @param actual a numeric vector of actuals.
#' @param forecast a numeric vector of forecasts.
#' @param na.rm a logical value indicating whether \code{NA} values should be stripped before the computation proceeds.
#' @param ... additional arguments to be passed to mean()
#' @return  a numeric vector of length one: the mean of of absolute percent errors.
#' @examples
#' # Examples
#' mape(1:100, 100:1)
#' @export
mape <- function(y_true = NULL, y_pred = NULL, actual = NULL, forecast = NULL, na.rm = TRUE, ...) {
  warning("'mape' is deprecated. Use 'mean_absolute_percent_error' instead.")
  if (is.null(actual)) {
    actual <- y_true
  } else {
    warning("Argument actual is deprecated; please use y_true instead.", call. = FALSE)
  }
  if (is.null(forecast)) {
    forecast <- y_pred
  } else {
    warning("Argument forecast is deprecated; please use y_pred instead.", call. = FALSE)
  }
  return(mean(ape(actual, forecast), na.rm = na.rm, ...))
}

#' Calculate root mean squared error.
#'
#' \code{rmse} takes actual and forecast numeric vectors and returns a
#' numeric vector where forecast is subtracted from the actual,
#' the errors are then squared, the average of those squared error is
#' taken, and lastly, the root taken of that mean squared error value.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @param actual a numeric vector of actuals.
#' @param forecast a numeric vector of forecasts.
#' @param na.rm a logical value indicating whether \code{NA} values should be stripped before the computation proceeds.
#' @param ... additional arguments to be passed to mean()
#' @return  a numeric vector of length one: the mean of of absolute percent errors.
#' @examples
#' # Examples
#' rmse(1:100, 100:1)
#' @export
rmse <- function(y_true = NULL, y_pred = NULL, actual = NULL, forecast = NULL, na.rm = TRUE, ...) {
  warning("'rmse' is deprecated. Use 'mean_squared_error' instead.")
  if (is.null(actual)) {
    actual <- y_true
  } else {
    warning("Argument actual is deprecated; please use y_true instead.", call. = FALSE)
  }
  if (is.null(forecast)) {
    forecast <- y_pred
  } else {
    warning("Argument forecast is deprecated; please use y_pred instead.", call. = FALSE)
  }
  return(sqrt(mean(e(actual, forecast) ^ 2, na.rm = na.rm, ...)))
}
