#' Convert an outlier_mask to a vector of column indices
#'
#' This function returns an ordered vector of all the 1-toggled bits in the outlier_mask offset by 1.
#' @param outlier_mask An integer representing the outlier mask.
#' @return A vector of column indices corresponding to the outlier mask.
#' @examples
#' outlier_mask <- 7
#' autovarCore:::selected_columns(outlier_mask)
#' @export
selected_columns <- function(outlier_mask) {
  result <- NULL
  for (column_index in 1:31)
    if (bitwAnd(outlier_mask, bitwShiftL(1, column_index - 1)))
      result <- c(result, column_index)
  result
}
