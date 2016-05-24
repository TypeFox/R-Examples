#' bond-equivalent yield (BEY), 2 x the semiannual discount rate
#' 
#' @param ear effective annual rate
#' @seealso \code{\link{ear}}
#' @export
#' @examples
#' ear2bey(ear=0.08)
ear2bey <- function(ear){
  return(((1 + ear)^0.5 - 1)*2)
}

