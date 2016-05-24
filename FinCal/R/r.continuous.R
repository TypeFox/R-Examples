#' Convert a given norminal rate to a continuous compounded rate 
#'
#' @param r norminal rate
#' @param m number of times compounded each year
#' @seealso \code{\link{r.norminal}}
#' @export
#' @examples
#' r.continuous(0.03,4)
r.continuous <- function(r,m) {
  return(m*log(1+r/m))
}

