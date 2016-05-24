#' Currency formatter: round to nearest penny and display pounds sign.
#'
#' The returned function will format a vector of values as currency.
#' Values are rounded to the nearest penny, and pennies are displayed if
#' any of the values has a non-zero pennies and the largest value is less
#' than \code{largest_with_penny} which by default is 100000.
#' 
#' Based heavily on the scales work by Hadley
#' 
#' @return a function with single paramater x, a numeric vector, that
#'   returns a character vector
#' @param largest_with_penny the value that all values of \code{x} must
#'   be less than in order for the cents to be displayed
#' @param x a numeric vector to format
#' @export
#' @examples
#' pounds_format()(c(100, 0.23, 1.456565, 2e3))
#' pounds_format()(c(1:10 * 10))
#' pounds(c(100, 0.23, 1.456565, 2e3))
#' pounds(c(1:10 * 10))
#' pounds(10^(1:8))
pounds_format <- function(x, largest_with_penny = 1e+05) {
    
    function(x) {
        x <- plyr::round_any(x, 0.01)
        if (max(x, na.rm = TRUE) < largest_with_penny & !all(x == floor(x), na.rm = TRUE)) {
            nsmall <- 2L
        } else {
            x <- plyr::round_any(x, 1)
            nsmall <- 0L
        }
        stringr::str_c("\u00a3", format(x, nsmall = nsmall, trim = TRUE, big.mark = ",", scientific = FALSE, digits = 1L))
    }
}

#' @export
#' @rdname pounds_format
pounds <- pounds_format() 
