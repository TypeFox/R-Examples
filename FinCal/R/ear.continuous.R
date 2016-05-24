#' Convert stated annual rate to the effective annual rate with continuous compounding
#'
#' @param r stated annual rate
#' @seealso \code{\link{ear}}
#' @seealso \code{\link{r.norminal}}
#' @export
#' @examples
#' ear.continuous(r=0.1)
#'
#' ear.continuous(0.03)
ear.continuous <- function(r){
  return(exp(r) - 1)
}

