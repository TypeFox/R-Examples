#'Function to format a number
#'
#'This function takes a floating number and round to 2 decimal point
#'
#'
#' @param x a floating number
#' @export format_number
#' @aliases format_number
#' @return formatted number
#' @keywords internal
#' 
#formatting a nubmber to decimal places
format_number <- function(x){
  format(round(x, 2), nsmall = 2)
}