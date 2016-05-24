#' Next Best Divisor
#' 
#' Given integers n and d, with n>d, this function finds the "next
#' best divisor" of n which is greater than or equal to d.
#' 
#' Suprisingly useful for thinking about processor grid shapes.
#' 
#' @param n
#' The divident (number divided into).
#' @param d
#' The candidate divisor.
#' 
#' @examples
#' \dontrun{
#' library(pbdBASE, quiet = TRUE)
#' base.nbd(100, 10) # 10 divides 100, so 10 is returned
#' base.nbd(100, 11) # 11 does not, so the "next best" divisor, 20, is returned
#' }
#' 
#' @export
base.nbd <- function(n, d)
{
  .Call(R_nbd, as.integer(n), as.integer(d))
}



isint <- function(x, epsilon=1e-8)
{
  if (is.numeric(x))
  {
    if (abs(x-as.integer(x)) < epsilon)
      return( TRUE )
    else
      return( FALSE )
  }
  else
    return( FALSE )
}

