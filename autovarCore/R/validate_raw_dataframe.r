#' Validates the dataframe given to the autovar function
#'
#' This function returns the given data frame as a numeric matrix, using \code{as.numeric} to convert any columns in the data frame that are not numeric. A \code{stop()} error is thrown if there is not enough data  in the data frame.
#' @param raw_dataframe The raw, unimputed data frame.
#' @return A numeric matrix with converted values and names taken from the data frame.
#' @examples
#' raw_dataframe <- data.frame(id = rep(1, times = 5),
#'   tijdstip = c(1, 3, 5, 6, 7),
#'   home = c(1, 0, 0, NA, 1))
#' autovarCore:::validate_raw_dataframe(raw_dataframe)
#' @export
validate_raw_dataframe <- function(raw_dataframe) {
  assert_param_not_null(raw_dataframe)
  assert_param_class(raw_dataframe, 'data.frame')
  assert_param_nrow(raw_dataframe, minimum = 1)
  for (column_name in names(raw_dataframe)) {
    if (class(raw_dataframe[[column_name]]) != 'numeric')
      raw_dataframe[[column_name]] <- as.numeric(raw_dataframe[[column_name]])
  }
  as.matrix(raw_dataframe)
}
