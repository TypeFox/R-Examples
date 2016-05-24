#' Thousands formatter: format number with commas separating the number thousands and suffixed with a k.
#' Based heavily on the scales work by Hadley
#' @param x a numeric vector to format
#' @return a function with single paramater x, a numeric vector, that
#'   returns a character vector
#' @export
#' @examples
#' thousands_format()(c(1, 1e3, 2000, 1e6))
#' thousands_format()(c(1, 1e3, 2000, 1e6))
#' thousands(c(1, 1e3, 2000, 1e6))
thousands_format <- function() {
    
    function(x) {
        x <- round(x/1000, 0)
        stringr::str_c(scales::comma(x), "k")
    }
}

#' @export
#' @rdname thousands_format
thousands <- thousands_format() 
