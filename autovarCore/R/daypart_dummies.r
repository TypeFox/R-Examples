#' Calculate day-part dummy variables
#'
#' This function returns either \code{NULL} (if \code{measurements_per_day} is 0 or 1) or a matrix of dummy variables for the specified input configuration.
#' @param number_of_rows the number of rows in the input data set.
#' @param measurements_per_day the number of measurements per day in the input data set.
#' @return Either \code{NULL} or a matrix with \code{number_of_rows} rows and \code{measurements_per_day - 1} columns.
#' @examples
#' autovarCore:::daypart_dummies(10, 3)
#' @export
daypart_dummies <- function(number_of_rows, measurements_per_day) {
  if (measurements_per_day == 0 || measurements_per_day == 1)
    return(NULL)
  seasonal_dummy_columns(out_length = number_of_rows,
                         period = measurements_per_day,
                         repetitions = 1,
                         dummy_name_prefix = 'dailymeas_')
}

seasonal_dummy_columns <- function(out_length, period, repetitions, dummy_name_prefix) {
  # This function is used by daypart_dummies() and day_dummies().
  result <- NULL
  required_column_count <- min(floor((out_length - 1) / repetitions), period - 1)
  if (required_column_count < 1)
    return(NULL)
  for (column_index in 1:required_column_count)
    result <- cbind(result, seasonal_dummy_column(out_length,
                                                  period,
                                                  repetitions,
                                                  column_index - 1))

  result <- as.matrix(result)
  colnames(result) <- dummy_column_names(ncol(result), dummy_name_prefix)
  result
}

seasonal_dummy_column <- function(out_length, period, repetitions, offset) {
  result <- c(rep.int(0, times = offset * repetitions),
              rep.int(1, times = repetitions),
              rep.int(0, times = repetitions * (period - offset - 1)))
  rep.int(result,
          times = ceiling(out_length / (period * repetitions)))[1:out_length]
}

dummy_column_names <- function(number_of_columns, dummy_name_prefix) {
  result <- NULL
  for (i in 1:number_of_columns)
    result <- c(result, paste(dummy_name_prefix, i, sep = ''))
  result
}
