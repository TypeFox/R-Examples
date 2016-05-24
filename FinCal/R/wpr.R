#' Weighted mean as a portfolio return
#' 
#' @param r returns of the individual assets in the portfolio
#' @param w corresponding weights associated with each of the individual assets
#' @export
#' @examples
#' wpr(r=c(0.12, 0.07, 0.03),w=c(0.5,0.4,0.1))
wpr <- function(r,w){
  if(sum(w) != 1){
    warning("sum of weights is NOT equal to 1!")
    return(sum(r*w))
  }else{
    return(sum(r*w))
  }
}
