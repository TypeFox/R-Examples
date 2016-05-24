#' Removes unwanted method from already existing methods in Test bench
#'
#' @param existing_method as Error observations for different methods
#' @param index_number as index number of unwanted method in study
#' @return it removes unwanted method from test bench and returns other method errors
#' @export
remove_method <- function(existing_method, index_number)
{
  existing_method[index_number+2] <- NULL
  return(existing_method)
}
