#' Convert a given continuous compounded rate to a norminal rate
#'
#' @param rc continuous compounded rate
#' @param m number of desired times compounded each year
#' @seealso \code{\link{r.continuous}}
#' @seealso \code{\link{ear.continuous}}
#' @export
#' @examples
#' r.norminal(0.03,1)
#'
#' r.norminal(0.03,4)
r.norminal <- function(rc,m) {
  return(m*(exp(rc/m)-1))
}

