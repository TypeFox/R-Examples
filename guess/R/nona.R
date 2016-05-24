#' No NAs
#'
#' Converts NAs to 0s
#' @param vec Required. Character or Numeric vector. 
#' @return Character vector.
#' @export
#' @examples
#' x <- c(NA, 1, 0); nona(x)
#' x <- c(NA, "dk", 0); nona(x)

nona <- function(vec = NULL)
{
    vec[is.na(vec)] <- 0 
    vec
}