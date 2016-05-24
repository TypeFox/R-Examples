#' Add Commas to a Large Number
#'
#' Convert a number to a string, with commas inserted at every 3rd digit.
#'
#' @param numbers Vector of non-negative numbers (will be rounded to integers)
#'
#' @return Character string with numbers written like \code{"5,771,009"}.
#'
#' @export
#'
#' @examples
#' commas(c(2300, 9000, 21456, 987654890, 1256787, 345765, 1432))

commas <- function(numbers) {
    format(numbers, big.mark = ",", scientific = FALSE, trim = TRUE)
}
