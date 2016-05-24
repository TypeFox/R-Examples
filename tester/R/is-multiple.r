#' @title Is multiple
#' @description Tests if x is multiple of a given number
#' @param x a numeric object
#' @param of a given number
#' @export
#' @examples
#' is_multiple(5, of = 5) # TRUE
#' is_multiple(15, of = 5) # TRUE
#' is_multiple(3, of = 5) # FALSE
#' is_multiple(2*pi, of = pi) # TRUE
#' is_multiple(matrix(1:6, 2, 3), of = 2)
is_multiple <- function(x, of) 
{
  if (is_not_scalar(of))
    stop("\n'of' must be a scalar")
  
  (x %% of) == 0
}
