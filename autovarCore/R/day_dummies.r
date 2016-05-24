#' Calculate day dummy variables
#'
#' This function returns either \code{NULL} (if \code{measurements_per_day} is 0) or a matrix of daydummy variables specified number of rows and measurements per day. In the latter case, we return a matrix of six columns.
#' @param number_of_rows the number of rows in the input data set.
#' @param measurements_per_day the number of measurements per day in the input data set.
#' @return Either \code{NULL} or a matrix with \code{number_of_rows} rows and \code{6} columns.
#' @examples
#' autovarCore:::day_dummies(16, 2)
#' @export
day_dummies <- function(number_of_rows, measurements_per_day) {
  if (measurements_per_day == 0)
    return(NULL)
  seasonal_dummy_columns(out_length = number_of_rows,
                         period = 7,
                         repetitions = measurements_per_day,
                         dummy_name_prefix = 'day_')
}
