#' Removes unwanted method from already existing methods in Test bench
#'
#' @param existing_method as Error observations for different methods
#' @param index_number as index number of unwanted method in study
#' @return it removes unwanted method from test bench and returns other method errors
#' @export
prediction_remove <- function(existing_method, index_number)
{
  existing_method[(3*(index_number - 1) + 3):(3*(index_number - 1) + 5)] <- NULL
  return(existing_method)
}
