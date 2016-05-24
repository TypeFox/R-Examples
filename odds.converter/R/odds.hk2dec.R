#' Convert Hong Kong Odds to Decimal Odds
#'
#' @param x A vector of Hong Kong odds
#'
#' @return A vector of Decimal odds
#'
#' @examples
#' odds.hk2dec(c(1.93,0.05))
odds.hk2dec <- function (x){
        ifelse (x <= 0,NA,x+1)
}