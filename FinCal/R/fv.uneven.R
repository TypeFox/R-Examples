#' Computing the future value of an uneven cash flow series
#'
#' @param r stated annual rate
#' @param cf uneven cash flow
#' @seealso \code{\link{fv.simple}}
#' @export
#' @examples
#' fv.uneven(r=0.1, cf=c(-1000, -500, 0, 4000, 3500, 2000))
fv.uneven <- function(r,cf){
  m <- length(cf)
  sum <- 0
  for(i in 1:m){
    n <- m - i
    sum <- sum + fv.simple(r,n,cf[i])
  }
  return(sum)
}

