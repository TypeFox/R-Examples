#' Computing bank discount yield (BDY) for a T-bill
#' 
#' @param d the dollar discount, which is equal to the difference between the face value of the bill and the purchase price 
#' @param f the face value (par value) of the bill
#' @param t number of days remaining until maturity
#' @seealso \code{\link{bdy2mmy}}
#' @export
#' @examples
#' bdy(d=1500,f=100000,t=120)
bdy <- function(d,f,t){
  return(360 * d / f / t)
}

