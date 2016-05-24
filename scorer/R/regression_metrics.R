#' Calculate absolute error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length equal to \code{max(length(y_true), length(y_pred))}.
#' @family regression metrics
#' @examples
#' absolute_error(1:10, 10:1)
#' @export
absolute_error <- function(y_true, y_pred) {
  return(absolute_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate percent error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length equal to \code{max(length(y_true), length(y_pred))}.
#' @family regression metrics
#' @examples
#' percent_error(1:10, 10:1)
#' @export
percent_error <- function(y_true, y_pred) {
  return(percent_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate log error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length equal to \code{max(length(y_true), length(y_pred))}.
#' @family regression metrics
#' @examples
#' log_error(1:10, 10:1)
#' @export
log_error <- function(y_true, y_pred) {
  return(log_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate squared error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length equal to \code{max(length(y_true), length(y_pred))}.
#' @family regression metrics
#' @examples
#' squared_error(1:10, 10:1)
#' @export
squared_error <- function(y_true, y_pred) {
  return(squared_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate squared log_error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length equal to \code{max(length(y_true), length(y_pred))}.
#' @family regression metrics
#' @examples
#' squared_log_error(1:10, 10:1)
#' @export
squared_log_error <- function(y_true, y_pred) {
  return(squared_log_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate absolute percent error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length equal to \code{max(length(y_true), length(y_pred))}.
#' @family regression metrics
#' @examples
#' absolute_percent_error(1:10, 10:1)
#' @export
absolute_percent_error <- function(y_true, y_pred) {
  return(absolute_percent_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate mean error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' mean_error(1:10, 10:1)
#' @export
mean_error <- function(y_true, y_pred) {
  return(mean_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate mean absolute error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' mean_absolute_error(1:10, 10:1)
#' @export
mean_absolute_error <- function(y_true, y_pred) {
  return(mean_absolute_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate median absolute error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' median_absolute_error(1:10, 10:1)
#' @export
median_absolute_error <- function(y_true, y_pred) {
  return(median_absolute_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate mean percent error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' mean_percent_error(1:10, 10:1)
#' @export
mean_percent_error <- function(y_true, y_pred) {
  return(mean_percent_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate median percent error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' median_percent_error(1:10, 10:1)
#' @export
median_percent_error <- function(y_true, y_pred) {
  return(median_percent_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate mean squared error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' mean_squared_error(1:10, 10:1)
#' @export
mean_squared_error <- function(y_true, y_pred) {
  return(mean_squared_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate median squared error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' median_squared_error(1:10, 10:1)
#' @export
median_squared_error <- function(y_true, y_pred) {
  return(median_squared_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate mean squared log error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' mean_squared_log_error(1:10, 10:1)
#' @export
mean_squared_log_error <- function(y_true, y_pred) {
  return(mean_squared_log_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate median squared log error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' median_squared_log_error(1:10, 10:1)
#' @export
median_squared_log_error <- function(y_true, y_pred) {
  return(median_squared_log_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate mean absolute percent error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' mean_absolute_percent_error(1:10, 10:1)
#' @export
mean_absolute_percent_error <- function(y_true, y_pred) {
  return(mean_absolute_percent_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate median absolute percent error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' median_absolute_percent_error(1:10, 10:1)
#' @export
median_absolute_percent_error <- function(y_true, y_pred) {
  return(median_absolute_percent_error_rcpp(y_true = y_true, y_pred = y_pred))
}

#' Calculate symmetric mean absolute percent error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' symmetric_mean_absolute_percent_error(1:10, 10:1)
#' @export
symmetric_mean_absolute_percent_error <- function(y_true, y_pred) {
  return(mean(abs(y_pred - y_true) / ((abs(y_true) + abs(y_pred)) / 2)))
}

#' Calculate symmetric median absolute percent error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' symmetric_median_absolute_percent_error(1:10, 10:1)
#' @export
symmetric_median_absolute_percent_error <- function(y_true, y_pred) {
  return(stats::median(abs(y_pred - y_true) / ((abs(y_true) + abs(y_pred)) / 2)))
}

#' Calculate mean absolute scaled error regression loss.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' mean_absolute_scaled_error(1:10, 10:1)
#' @export
mean_absolute_scaled_error <- function(y_true, y_pred) {
  n <- max(length(y_true), length(y_pred))
  numerator <- sum(abs(y_true - y_pred))
  denominator <- (n / (n - 1)) * sum(abs(y_true[2:n] - y_pred[1:(n-1)]))
  return(numerator / denominator)
}

#' Calculate total variance regression score function.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' total_variance_score(1:10, 10:1)
#' @export
total_variance_score <- function(y_true, y_pred) {
  return(sum((y_true - mean(y_true)) ^ 2))
}

#' Calculate explained variance regression score function.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' explained_variance_score(1:10, 10:1)
#' @export
explained_variance_score <- function(y_true, y_pred) {
  return(sum((y_pred - mean(y_true)) ^ 2))
}

#' Calculate unexplained variance regression score function.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' unexplained_variance_score(1:10, 10:1)
#' @export
unexplained_variance_score <- function(y_true, y_pred) {
  return(sum((y_true - y_pred) ^ 2))
}

#' Calculate R^2 (coefficient of determination) regression score function.
#'
#' @param y_true Ground truth (correct) target values.
#' @param y_pred Estimated target values.
#' @return  A numeric vector of length one.
#' @family regression metrics
#' @examples
#' r2_score(1:10, 10:1)
#' @export
r2_score <- function(y_true, y_pred) {
  return(explained_variance_score(y_true = y_true, y_pred = y_pred) / total_variance_score(y_true = y_true, y_pred = y_pred))
}
